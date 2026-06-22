use crate::error::{Error, Result};
use crate::hit::{hit_from_raw, Hit};
use crate::index::Index;
use crate::meth::Meth;
use crate::opts::Opts;
use minibwa_sys as sys;
use std::ffi::CString;

/// Per-thread scratch buffer for alignment. One per worker thread.
///
/// Created with a kalloc pool (`mb_tbuf_init(0)`) so it is valid for both
/// single-read and batch alignment paths.
///
/// # Panics
///
/// [`ThreadBuf::new`] panics if `mb_tbuf_init` returns NULL (out of memory).
pub struct ThreadBuf {
    ptr: *mut sys::mb_tbuf_t,
}

// SAFETY: the buffer owns its own kalloc arena and shares nothing; it may move
// between threads but must not be used from two threads at once (not Sync).
unsafe impl Send for ThreadBuf {}

impl ThreadBuf {
    pub fn new() -> ThreadBuf {
        // SAFETY: mb_tbuf_init(0) returns an owned buffer with a kalloc pool.
        let ptr = unsafe { sys::mb_tbuf_init(0) };
        assert!(!ptr.is_null(), "mb_tbuf_init returned NULL");
        ThreadBuf { ptr }
    }
}

impl Default for ThreadBuf {
    fn default() -> Self {
        ThreadBuf::new()
    }
}

impl Drop for ThreadBuf {
    fn drop(&mut self) {
        // SAFETY: ptr came from mb_tbuf_init and is destroyed exactly once.
        unsafe { sys::mb_tbuf_destroy(self.ptr) };
    }
}

/// An aligner binding an index and options. Cheap to create; the index is shared.
pub struct Aligner<'a> {
    idx: &'a Index,
    opts: &'a Opts,
}

impl<'a> Aligner<'a> {
    /// Create an aligner from a shared index and options reference.
    pub fn new(idx: &'a Index, opts: &'a Opts) -> Aligner<'a> {
        Aligner { idx, opts }
    }

    /// Align a read pair (paired-end). Returns `(hits_r1, hits_r2)`. Pairing,
    /// mate rescue, and proper-pair flagging use minibwa's PE path with the
    /// insert-size parameters on `Opts`. If the options enable methylation,
    /// R1 is treated as C2T and R2 as G2A automatically; this requires `idx`
    /// to have been built and loaded with `meth = true`, otherwise the
    /// alignments are silently wrong.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] if either sequence is empty or either
    /// read name contains a NUL byte.
    pub fn map_pair<N1, S1, N2, S2>(
        &self,
        buf: &mut ThreadBuf,
        name1: N1,
        seq1: S1,
        name2: N2,
        seq2: S2,
    ) -> Result<(Vec<Hit>, Vec<Hit>)>
    where
        N1: AsRef<[u8]>,
        S1: AsRef<[u8]>,
        N2: AsRef<[u8]>,
        S2: AsRef<[u8]>,
    {
        let name1 = name1.as_ref();
        let seq1 = seq1.as_ref();
        let name2 = name2.as_ref();
        let seq2 = seq2.as_ref();

        if seq1.is_empty() || seq2.is_empty() {
            return Err(Error::InvalidInput("empty query sequence".into()));
        }
        let c_name1 = CString::new(name1)
            .map_err(|_| Error::InvalidInput("read name contains NUL".into()))?;
        let c_name2 = CString::new(name2)
            .map_err(|_| Error::InvalidInput("read name contains NUL".into()))?;

        let pe_opt = self.opts.paired_copy();
        // SAFETY: seq1/seq2 are valid byte slices; the C code reads exactly qlen
        // bytes and does not require NUL termination on the sequence data.
        let seq1_ptr = seq1.as_ptr() as *const std::os::raw::c_char;
        let seq2_ptr = seq2.as_ptr() as *const std::os::raw::c_char;
        let mut seqs: [*const std::os::raw::c_char; 2] = [seq1_ptr, seq2_ptr];
        let mut names: [*const std::os::raw::c_char; 2] = [c_name1.as_ptr(), c_name2.as_ptr()];
        let lens: [i32; 2] = [seq1.len() as i32, seq2.len() as i32];
        let mut n_hit: [i32; 2] = [0, 0];

        // SAFETY: all arrays outlive the call; n_seq=2 with MB_F_PE pairs the
        // two reads. Returns a libc-allocated mb_hit_t** (2 slots).
        let raw = unsafe {
            sys::mb_map_batch(
                &pe_opt,
                self.idx.as_ptr(),
                2,
                lens.as_ptr(),
                seqs.as_mut_ptr(),
                n_hit.as_mut_ptr(),
                buf.ptr,
                names.as_mut_ptr(),
            )
        };
        if raw.is_null() {
            return Ok((Vec::new(), Vec::new()));
        }

        // Pass 1: copy each read's hits into owned Vecs (reads *.p, frees nothing).
        let mut out: [Vec<Hit>; 2] = [Vec::new(), Vec::new()];
        for (k, out_k) in out.iter_mut().enumerate() {
            // SAFETY: raw has 2 mb_hit_t* slots.
            let arr = unsafe { *raw.add(k) };
            if arr.is_null() || n_hit[k] <= 0 {
                continue;
            }
            // SAFETY: arr has n_hit[k] contiguous mb_hit_t.
            let slice = unsafe { std::slice::from_raw_parts(arr, n_hit[k] as usize) };
            out_k.reserve(slice.len());
            for h in slice {
                // SAFETY: h is a valid hit; hit_from_raw reads only, frees nothing.
                out_k.push(unsafe { hit_from_raw(h, self.idx) });
            }
        }

        // Pass 2: free libc memory — each .p, each per-read array, then the outer array.
        for (k, &count) in n_hit.iter().enumerate() {
            // SAFETY: raw has 2 slots.
            let arr = unsafe { *raw.add(k) };
            if arr.is_null() {
                continue;
            }
            if count > 0 {
                // SAFETY: arr has count contiguous mb_hit_t.
                let slice = unsafe { std::slice::from_raw_parts(arr, count as usize) };
                for h in slice {
                    if !h.p.is_null() {
                        // SAFETY: .p is libc-allocated, freed exactly once.
                        unsafe { libc::free(h.p as *mut libc::c_void) };
                    }
                }
            }
            // SAFETY: arr is libc-allocated, freed exactly once.
            unsafe { libc::free(arr as *mut libc::c_void) };
        }
        // SAFETY: raw is the libc-allocated outer array, freed exactly once.
        unsafe { libc::free(raw as *mut libc::c_void) };

        let [r1, r2] = out;
        Ok((r1, r2))
    }

    /// Align many reads in one batched call, returning per-read hit lists.
    ///
    /// This is the high-throughput entry point: minibwa batches the BWT
    /// seeding across the whole slice (prefetch-driven), which the single-read
    /// [`Aligner::map`] cannot do. `names` and `seqs` must be the same length;
    /// the result has one `Vec<Hit>` per input read, in order.
    ///
    /// Paired-end mode is forced off (these are independent reads). If the
    /// options enable methylation, every read is treated as the C2T (read-1)
    /// strand — for G2A or mixed methylation use [`Aligner::map`] per read.
    /// Methylation requires `idx` to have been built and loaded with
    /// `meth = true`, otherwise the alignments are silently wrong.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] if `names.len() != seqs.len()`, any
    /// sequence is empty, or any name contains a NUL byte.
    pub fn map_many(
        &self,
        buf: &mut ThreadBuf,
        names: &[impl AsRef<[u8]>],
        seqs: &[impl AsRef<[u8]>],
    ) -> Result<Vec<Vec<Hit>>> {
        if names.len() != seqs.len() {
            return Err(Error::InvalidInput(format!(
                "names ({}) and seqs ({}) length mismatch",
                names.len(),
                seqs.len()
            )));
        }
        let n = names.len();
        if n == 0 {
            return Ok(Vec::new());
        }
        // Force PE off for independent-read batching (Opts is Clone; set_paired consumes).
        let opt = self.opts.clone().set_paired(false);

        // Hold owned CStrings for names (must outlive the call); seqs are zero-copy.
        let c_names: Vec<CString> = names
            .iter()
            .map(|nm| {
                CString::new(nm.as_ref())
                    .map_err(|_| Error::InvalidInput("read name contains NUL".into()))
            })
            .collect::<Result<_>>()?;
        let mut name_ptrs: Vec<*const std::os::raw::c_char> =
            c_names.iter().map(|c| c.as_ptr()).collect();
        let mut seq_ptrs: Vec<*const std::os::raw::c_char> = Vec::with_capacity(n);
        let mut lens: Vec<i32> = Vec::with_capacity(n);
        for s in seqs {
            let s = s.as_ref();
            if s.is_empty() {
                return Err(Error::InvalidInput("empty query sequence".into()));
            }
            seq_ptrs.push(s.as_ptr() as *const std::os::raw::c_char);
            lens.push(s.len() as i32);
        }
        let mut n_hit: Vec<i32> = vec![0; n];

        // SAFETY: all arrays/CStrings/slices outlive the call; n_seq = n with PE off
        // aligns each read independently. Returns a libc-allocated mb_hit_t** (n slots).
        let raw = unsafe {
            sys::mb_map_batch(
                opt.as_ptr(),
                self.idx.as_ptr(),
                n as i32,
                lens.as_ptr(),
                seq_ptrs.as_mut_ptr(),
                n_hit.as_mut_ptr(),
                buf.ptr,
                name_ptrs.as_mut_ptr(),
            )
        };
        if raw.is_null() {
            return Ok(vec![Vec::new(); n]);
        }

        // Pass 1: copy every read's hits into owned Vecs (reads *.p, frees nothing).
        let mut out: Vec<Vec<Hit>> = Vec::with_capacity(n);
        for (k, &count) in n_hit.iter().enumerate() {
            // SAFETY: raw has n mb_hit_t* slots.
            let arr = unsafe { *raw.add(k) };
            if arr.is_null() || count <= 0 {
                out.push(Vec::new());
                continue;
            }
            // SAFETY: arr has count contiguous mb_hit_t.
            let slice = unsafe { std::slice::from_raw_parts(arr, count as usize) };
            let mut hits = Vec::with_capacity(slice.len());
            for h in slice {
                // SAFETY: h is a valid hit; hit_from_raw reads only, frees nothing.
                hits.push(unsafe { hit_from_raw(h, self.idx) });
            }
            out.push(hits);
        }

        // Pass 2: free libc memory — each .p, each per-read array, then the outer array.
        for (k, &count) in n_hit.iter().enumerate() {
            // SAFETY: raw has n slots.
            let arr = unsafe { *raw.add(k) };
            if arr.is_null() {
                continue;
            }
            if count > 0 {
                // SAFETY: arr has count contiguous mb_hit_t.
                let slice = unsafe { std::slice::from_raw_parts(arr, count as usize) };
                for h in slice {
                    if !h.p.is_null() {
                        // SAFETY: .p is libc-allocated, freed exactly once.
                        unsafe { libc::free(h.p as *mut libc::c_void) };
                    }
                }
            }
            // SAFETY: arr is libc-allocated, freed exactly once.
            unsafe { libc::free(arr as *mut libc::c_void) };
        }
        // SAFETY: raw is the libc-allocated outer array, freed exactly once.
        unsafe { libc::free(raw as *mut libc::c_void) };

        Ok(out)
    }

    /// Align one read, returning owned hits. `name` must be NUL-free.
    ///
    /// A non-[`Meth::None`] value requires `idx` to have been built and loaded
    /// with `meth = true`; otherwise the conversion is applied against an
    /// unconverted reference and the alignments are silently wrong.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] if `seq` is empty or `name` contains a
    /// NUL byte.
    pub fn map<N: AsRef<[u8]>, S: AsRef<[u8]>>(
        &self,
        buf: &mut ThreadBuf,
        name: N,
        seq: S,
        meth: Meth,
    ) -> Result<Vec<Hit>> {
        let name = name.as_ref();
        let seq = seq.as_ref();

        if seq.is_empty() {
            return Err(Error::InvalidInput("empty query sequence".into()));
        }
        let c_name =
            CString::new(name).map_err(|_| Error::InvalidInput("read name contains NUL".into()))?;

        // SAFETY: seq is a valid byte slice; the C code reads exactly seq.len()
        // bytes and does not require NUL termination on the sequence data.
        let seq_ptr = seq.as_ptr() as *const std::os::raw::c_char;
        let mut n_hit: i32 = 0;
        // SAFETY: all pointers valid for the duration of the call; mb_map allocates
        // the returned array and each .p with libc malloc/realloc.
        let raw = unsafe {
            sys::mb_map(
                self.opts.as_ptr(),
                self.idx.as_ptr(),
                seq.len() as i32,
                seq_ptr,
                meth.as_mt(),
                &mut n_hit,
                buf.ptr,
                c_name.as_ptr(),
            )
        };

        if raw.is_null() || n_hit <= 0 {
            if !raw.is_null() {
                // SAFETY: empty but allocated array still needs freeing.
                unsafe { libc::free(raw as *mut libc::c_void) };
            }
            return Ok(Vec::new());
        }

        let n = n_hit as usize;
        let mut out = Vec::with_capacity(n);
        // SAFETY: raw points to n_hit contiguous mb_hit_t.
        let slice = unsafe { std::slice::from_raw_parts(raw, n) };

        for h in slice {
            // SAFETY: h is a valid hit from mb_map; hit_from_raw reads but does not free.
            out.push(unsafe { hit_from_raw(h, self.idx) });
        }

        // Free libc-allocated C memory: each .p pointer, then the array itself.
        for h in slice {
            if !h.p.is_null() {
                // SAFETY: .p is libc-allocated by mb_map and freed exactly once here.
                unsafe { libc::free(h.p as *mut libc::c_void) };
            }
        }
        // SAFETY: raw is libc-allocated by mb_map and freed exactly once here.
        unsafe { libc::free(raw as *mut libc::c_void) };

        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_util::{synthetic_reference, write_fasta};
    use crate::{Index, Meth, Opts};

    fn revcomp(s: &[u8]) -> Vec<u8> {
        s.iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'T' => b'A',
                b'C' => b'G',
                b'G' => b'C',
                _ => b'N',
            })
            .collect()
    }

    fn build_test_index(tag: &str) -> (std::path::PathBuf, Vec<u8>, Index, Opts) {
        let dir = std::env::temp_dir().join(format!("minibwa_aln_{tag}_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let seq = synthetic_reference(5000, 13);
        let fa = write_fasta(&dir, "chr1", &seq);
        let prefix = dir.join("idx");
        Index::build_from_fasta(&fa, &prefix, false, 1).unwrap();
        let idx = Index::load(&prefix, false).unwrap();
        let opts = Opts::new();
        (dir, seq, idx, opts)
    }

    #[test]
    fn map_rejects_empty_seq() {
        let (dir, _seq, idx, opts) = build_test_index("empty_seq");
        let aligner = Aligner::new(&idx, &opts);
        let mut buf = ThreadBuf::new();
        let res = aligner.map(&mut buf, b"r", &[][..], Meth::None);
        assert!(matches!(res, Err(Error::InvalidInput(_))));
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn map_rejects_nul_in_name() {
        let (dir, seq, idx, opts) = build_test_index("nul_name");
        let aligner = Aligner::new(&idx, &opts);
        let mut buf = ThreadBuf::new();
        let res = aligner.map(&mut buf, b"bad\0name", &seq[0..100], Meth::None);
        assert!(matches!(res, Err(Error::InvalidInput(_))));
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn map_pair_rejects_empty_r2() {
        let (dir, seq, idx, opts) = build_test_index("pe_empty_r2");
        let aligner = Aligner::new(&idx, &opts);
        let mut buf = ThreadBuf::new();
        let res = aligner.map_pair(&mut buf, b"r1", &seq[0..100], b"r2", &[][..]);
        assert!(matches!(res, Err(Error::InvalidInput(_))));
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn map_reverse_strand() {
        let (dir, seq, idx, opts) = build_test_index("rev_strand");
        let aligner = Aligner::new(&idx, &opts);
        let mut buf = ThreadBuf::new();
        let rc = revcomp(&seq[2000..2150]);
        let hits = aligner.map(&mut buf, b"rc", &rc, Meth::None).unwrap();
        assert!(
            !hits.is_empty(),
            "reverse-complement query should hit reference"
        );
        assert_eq!(
            hits[0].strand,
            crate::Strand::Reverse,
            "expected Reverse strand for revcomp query"
        );
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn map_hit_fields_populated() {
        let (dir, seq, idx, opts) = build_test_index("hit_fields");
        let aligner = Aligner::new(&idx, &opts);
        let mut buf = ThreadBuf::new();
        let query = &seq[1000..1150];
        let hits = aligner.map(&mut buf, b"fwd", query, Meth::None).unwrap();
        assert!(!hits.is_empty(), "forward query should produce hits");
        let h = &hits[0];
        assert!(h.ref_end > h.ref_start, "ref_end should be > ref_start");
        assert_eq!(h.query_start, 0, "query_start should be 0");
        assert_eq!(
            h.query_end as usize,
            query.len(),
            "query_end should equal query length"
        );
        assert_eq!(h.strand, crate::Strand::Forward);
        assert!(h.mapq > 0, "mapq should be non-zero for a unique hit");
        assert!(h.score > 0, "alignment score should be positive");
        let first_kind = h.cigar[0].kind;
        assert!(
            first_kind == crate::hit::CigarKind::Match || first_kind == crate::hit::CigarKind::Eq,
            "first CIGAR op should be Match or Eq, got {first_kind:?}"
        );
        assert!(!h.cigar_string().is_empty());
        assert_ne!(h.cigar_string(), "*");
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn map_one_read_hits_reference() {
        let dir = std::env::temp_dir().join(format!("minibwa_aln_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let seq = synthetic_reference(5000, 13);
        let fa = write_fasta(&dir, "chr1", &seq);
        let prefix = dir.join("idx");
        Index::build_from_fasta(&fa, &prefix, false, 1).unwrap();
        let idx = Index::load(&prefix, false).unwrap();
        let opts = Opts::new();
        let aligner = Aligner::new(&idx, &opts);
        let mut buf = ThreadBuf::new();

        let query = &seq[1000..1150];
        let hits = aligner.map(&mut buf, b"q1", query, Meth::None).unwrap();
        assert!(!hits.is_empty());
        let h = &hits[0];
        assert_eq!(h.contig.as_deref(), Some("chr1"));
        assert!(h.ref_start >= 900 && h.ref_start <= 1100);
        assert!(!h.cigar.is_empty());
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn map_pair_proper_fr() {
        // Use insert-size opts that bracket the ~400 bp insert in this pair so
        // the proper-pair flag is reliably set by minibwa's PE scorer.
        let (dir, seq, idx, _) = build_test_index("pe_fr");
        let opts = Opts::new().set_pe_insert_size(400, 50, 200, 600);
        let aligner = Aligner::new(&idx, &opts);
        let mut buf = ThreadBuf::new();

        let r1 = seq[1000..1150].to_vec();
        let r2 = revcomp(&seq[1400..1550]);

        let (h1, h2) = aligner.map_pair(&mut buf, b"p1", &r1, b"p2", &r2).unwrap();
        assert!(!h1.is_empty(), "R1 should have hits");
        assert!(!h2.is_empty(), "R2 should have hits");
        assert_eq!(h1[0].contig.as_deref(), Some("chr1"));
        assert!(
            h1[0].ref_start >= 900 && h1[0].ref_start <= 1100,
            "R1 ref_start={} expected in [900,1100]",
            h1[0].ref_start
        );
        assert_eq!(h2[0].contig.as_deref(), Some("chr1"));
        assert!(
            h2[0].ref_start >= 1300 && h2[0].ref_start <= 1500,
            "R2 ref_start={} expected in [1300,1500]",
            h2[0].ref_start
        );
        // Strand: R1 forward (matches reference directly), R2 reverse (revcomp).
        assert_eq!(h1[0].strand, crate::Strand::Forward, "R1 should be Forward");
        assert_eq!(h2[0].strand, crate::Strand::Reverse, "R2 should be Reverse");
        // Both primary (first and only alignment for each read).
        assert!(h1[0].is_primary, "R1 primary hit should be primary");
        assert!(h2[0].is_primary, "R2 primary hit should be primary");
        // Proper-pair: insert-size window covers ~400 bp, so this should fire.
        assert!(h1[0].proper_pair, "R1 should be flagged as proper pair");
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn map_no_hit_returns_empty() {
        let (dir, _seq, idx, opts) = build_test_index("no_hit");
        let aligner = Aligner::new(&idx, &opts);
        let mut buf = ThreadBuf::new();

        let query = vec![b'N'; 150];
        let hits = aligner.map(&mut buf, b"nh", &query, Meth::None).unwrap();
        assert!(hits.is_empty(), "all-N query should return no hits");
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn map_many_batch() {
        let (dir, seq, idx, opts) = build_test_index("batch");
        let aligner = Aligner::new(&idx, &opts);
        let mut buf = ThreadBuf::new();

        let q0 = seq[500..650].to_vec();
        let q1 = seq[2000..2150].to_vec();
        let q2 = seq[3500..3650].to_vec();

        let results = aligner
            .map_many(&mut buf, &[b"r0" as &[u8], b"r1", b"r2"], &[&q0, &q1, &q2])
            .expect("map_many failed");

        assert_eq!(results.len(), 3);
        for (i, hits) in results.iter().enumerate() {
            assert!(!hits.is_empty(), "read r{i} produced no hits");
        }
        let expected_offsets: [i64; 3] = [500, 2000, 3500];
        for (i, (hits, &exp)) in results.iter().zip(expected_offsets.iter()).enumerate() {
            let pos = hits[0].ref_start;
            assert!(
                (pos - exp).abs() <= 100,
                "read r{i}: expected pos ~{exp}, got {pos}"
            );
        }
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn map_many_length_mismatch() {
        let (dir, seq, idx, opts) = build_test_index("mismatch");
        let aligner = Aligner::new(&idx, &opts);
        let mut buf = ThreadBuf::new();

        let q0 = seq[500..650].to_vec();
        let q1 = seq[2000..2150].to_vec();

        let result = aligner.map_many(&mut buf, &[b"r0" as &[u8]], &[&q0, &q1]);
        assert!(
            matches!(result, Err(Error::InvalidInput(_))),
            "expected InvalidInput error for length mismatch"
        );
        std::fs::remove_dir_all(&dir).ok();
    }
}
