use crate::error::{Error, Result};
use minibwa_sys as sys;
use std::ffi::{CStr, CString};
use std::path::Path;

/// A loaded minibwa index. Immutable after load; safe to share across threads.
pub struct Index {
    ptr: *mut sys::mb_idx_t,
}

// SAFETY: the index is read-only during mapping. mb_idx_t = {is_meth, l2b*, bwt*};
// mb_map copies opt into a local and only reads idx; the one mutable global on the
// map path (kom_nt4_table) is a read-only lookup table. minibwa's own pipeline shares
// a single const mb_idx_t* across worker threads. No interior mutability is exposed.
unsafe impl Send for Index {}
unsafe impl Sync for Index {}

fn c_path(p: &Path) -> Result<CString> {
    CString::new(p.to_string_lossy().as_bytes())
        .map_err(|_| Error::InvalidInput(format!("path has NUL: {}", p.display())))
}

impl Index {
    /// Build an index from a FASTA, writing `<prefix>.l2b`, `<prefix>.mbw`,
    /// and (when `meth`) `<prefix>.meth.mbw`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::IndexBuild`] if the FASTA cannot be read, is empty, the
    /// output directory does not exist, or the underlying C indexer reports an
    /// error.  Returns [`Error::InvalidInput`] if either path contains a NUL byte.
    pub fn build_from_fasta(
        fasta: impl AsRef<Path>,
        prefix: impl AsRef<Path>,
        meth: bool,
        threads: u32,
    ) -> Result<()> {
        let fasta = fasta.as_ref();
        let prefix = prefix.as_ref();

        // Pre-validate before the FFI call: minibwa's main_index() abort()s on
        // unreadable input, so we must catch the common failures in Rust first.
        let meta = std::fs::metadata(fasta).map_err(|e| Error::IndexBuild {
            fasta: fasta.to_path_buf(),
            msg: format!("cannot stat FASTA: {e}"),
        })?;
        if meta.len() == 0 {
            return Err(Error::IndexBuild {
                fasta: fasta.to_path_buf(),
                msg: "FASTA is empty".into(),
            });
        }
        if let Some(parent) = prefix.parent() {
            if !parent.as_os_str().is_empty() && !parent.is_dir() {
                return Err(Error::IndexBuild {
                    fasta: fasta.to_path_buf(),
                    msg: format!("output directory does not exist: {}", parent.display()),
                });
            }
        }

        let c_fa = c_path(fasta)?;
        let c_prefix = c_path(prefix)?;
        // SAFETY: both pointers are valid NUL-terminated paths; shim returns nonzero on error.
        let rc = unsafe {
            sys::mb_index_build(
                c_fa.as_ptr(),
                c_prefix.as_ptr(),
                meth as i32,
                threads as i32,
            )
        };
        if rc != 0 {
            return Err(Error::IndexBuild {
                fasta: fasta.to_path_buf(),
                msg: shim_last_error(),
            });
        }
        Ok(())
    }

    /// Load an index previously built at `prefix`. `meth` selects the 3-base index.
    ///
    /// # Errors
    ///
    /// Returns [`Error::IndexLoad`] if `mb_idx_load` returns NULL (missing or
    /// corrupt index files).  Returns [`Error::InvalidInput`] if the path
    /// contains a NUL byte.
    pub fn load(prefix: impl AsRef<Path>, meth: bool) -> Result<Index> {
        let prefix = prefix.as_ref();
        let c_prefix = c_path(prefix)?;
        // SAFETY: valid NUL-terminated prefix; mb_idx_load returns NULL on failure.
        let ptr = unsafe { sys::mb_idx_load(c_prefix.as_ptr(), meth as i32) };
        if ptr.is_null() {
            return Err(Error::IndexLoad {
                path: prefix.to_path_buf(),
                msg: "mb_idx_load returned NULL (missing or corrupt index files)".into(),
            });
        }
        Ok(Index { ptr })
    }

    /// Number of reference contigs.
    pub fn n_contigs(&self) -> usize {
        let mut n = 0usize;
        // Count by the C bounds-check (mb_idx_ctg_name returns NULL past the last
        // contig), independent of whether each name is valid UTF-8.
        // SAFETY: valid index pointer; mb_idx_ctg_name bounds-checks tid.
        while !unsafe { sys::mb_idx_ctg_name(self.ptr, n as i32) }.is_null() {
            n += 1;
        }
        n
    }

    /// Name of contig `tid`, or `None` if out of range.
    pub fn contig_name(&self, tid: usize) -> Option<&str> {
        // SAFETY: valid index; mb_idx_ctg_name bounds-checks and returns NULL out of range.
        let p = unsafe { sys::mb_idx_ctg_name(self.ptr, tid as i32) };
        if p.is_null() {
            return None;
        }
        // SAFETY: non-null contig name is a NUL-terminated string owned by the index.
        unsafe { CStr::from_ptr(p) }.to_str().ok()
    }

    /// Length of contig `tid`, or `None` if out of range.
    pub fn contig_len(&self, tid: usize) -> Option<i64> {
        // SAFETY: valid index; returns -1 out of range.
        let len = unsafe { sys::mb_idx_ctg_len(self.ptr, tid as i32) };
        if len < 0 {
            None
        } else {
            Some(len)
        }
    }

    /// Iterate `(name, len)` over all contigs.
    pub fn contigs(&self) -> impl Iterator<Item = (&str, i64)> {
        (0..self.n_contigs()).map(move |tid| {
            (
                self.contig_name(tid).unwrap(),
                self.contig_len(tid).unwrap(),
            )
        })
    }

    pub(crate) fn as_ptr(&self) -> *const sys::mb_idx_t {
        self.ptr
    }
}

impl Drop for Index {
    fn drop(&mut self) {
        // SAFETY: ptr came from mb_idx_load and is freed exactly once.
        unsafe { sys::mb_idx_destroy(self.ptr) };
    }
}

fn shim_last_error() -> String {
    // SAFETY: returns a thread-local NUL-terminated string (possibly empty).
    let p = unsafe { sys::mb_shim_last_error() };
    if p.is_null() {
        return String::new();
    }
    unsafe { CStr::from_ptr(p) }.to_string_lossy().into_owned()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_util::{synthetic_reference, write_fasta, write_multi_fasta};

    fn tmpdir(tag: &str) -> std::path::PathBuf {
        let d = std::env::temp_dir().join(format!("minibwa_idx_{tag}_{}", std::process::id()));
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    #[test]
    fn build_load_contigs() {
        let dir = tmpdir("blc");
        let seq = synthetic_reference(3000, 7);
        let fa = write_fasta(&dir, "chr1", &seq);
        let prefix = dir.join("idx");

        Index::build_from_fasta(&fa, &prefix, false, 1).unwrap();
        let idx = Index::load(&prefix, false).unwrap();
        assert_eq!(idx.n_contigs(), 1);
        assert_eq!(idx.contig_name(0), Some("chr1"));
        assert_eq!(idx.contig_len(0), Some(3000));
        assert_eq!(idx.contig_name(1), None);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn multi_contig_index() {
        let dir = tmpdir("multi");
        let seq1 = synthetic_reference(2000, 17);
        let seq2 = synthetic_reference(3000, 23);
        let fa = write_multi_fasta(&dir, &[("chr1", &seq1), ("chr2", &seq2)]);
        let prefix = dir.join("idx");

        Index::build_from_fasta(&fa, &prefix, false, 1).unwrap();
        let idx = Index::load(&prefix, false).unwrap();

        assert_eq!(idx.n_contigs(), 2);
        assert_eq!(idx.contig_name(0), Some("chr1"));
        assert_eq!(idx.contig_len(0), Some(2000i64));
        assert_eq!(idx.contig_name(1), Some("chr2"));
        assert_eq!(idx.contig_len(1), Some(3000i64));
        assert_eq!(idx.contig_len(2), None);
        assert_eq!(
            idx.contigs().collect::<Vec<_>>(),
            vec![("chr1", 2000i64), ("chr2", 3000i64)]
        );
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn build_empty_fasta_errors() {
        let dir = tmpdir("empty");
        let fa = dir.join("empty.fa");
        std::fs::File::create(&fa).unwrap(); // zero-byte file
        let res = Index::build_from_fasta(&fa, dir.join("idx"), false, 1);
        assert!(matches!(res, Err(crate::Error::IndexBuild { .. })));
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn build_bad_output_dir_errors() {
        let dir = tmpdir("baddir");
        let seq = synthetic_reference(2000, 7);
        let fa = write_fasta(&dir, "chr1", &seq);
        // Prefix inside a nonexistent subdirectory.
        let res = Index::build_from_fasta(&fa, dir.join("nodir/idx"), false, 1);
        assert!(matches!(res, Err(crate::Error::IndexBuild { .. })));
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn build_missing_fasta_errors_not_aborts() {
        let dir = tmpdir("missing");
        let res = Index::build_from_fasta(dir.join("does-not-exist.fa"), dir.join("idx"), false, 1);
        assert!(matches!(res, Err(crate::Error::IndexBuild { .. })));
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn load_missing_errors() {
        let dir = tmpdir("loadmiss");
        let res = Index::load(dir.join("nope"), false);
        assert!(matches!(res, Err(crate::Error::IndexLoad { .. })));
        std::fs::remove_dir_all(&dir).ok();
    }
}
