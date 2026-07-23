//! End-to-end raw-FFI smoke test: build an index from a synthetic FASTA,
//! load it, align one read, and verify a hit lands on the reference.
use std::ffi::{CStr, CString};
use std::io::Write;

fn write_fasta(dir: &std::path::Path) -> (std::path::PathBuf, Vec<u8>) {
    // 2 kb pseudo-random ACGT reference, deterministic.
    let mut seq = Vec::with_capacity(2000);
    let mut x: u32 = 0x1234_5678;
    for _ in 0..2000 {
        x = x.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
        seq.push(b"ACGT"[(x >> 16) as usize & 3]);
    }
    let fa = dir.join("ref.fa");
    let mut f = std::fs::File::create(&fa).unwrap();
    writeln!(f, ">chr1").unwrap();
    f.write_all(&seq).unwrap();
    writeln!(f).unwrap();
    (fa, seq)
}

/// Build the synthetic reference into a fresh temp dir and load it. Returns the temp dir
/// (caller removes it), the reference bases, and the loaded index (caller destroys it).
fn build_and_load_index(tag: &str) -> (std::path::PathBuf, Vec<u8>, *mut minibwa_sys::mb_idx_t) {
    let dir = std::env::temp_dir().join(format!("minibwa_sys_{}_{}", tag, std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let (fa, refseq) = write_fasta(&dir);
    let prefix = dir.join("ref");
    let c_fa = CString::new(fa.to_str().unwrap()).unwrap();
    let c_prefix = CString::new(prefix.to_str().unwrap()).unwrap();
    unsafe {
        let rc = minibwa_sys::mb_index_build(c_fa.as_ptr(), c_prefix.as_ptr(), 0, 1);
        assert_eq!(rc, 0, "index build failed");
        let idx = minibwa_sys::mb_idx_load(c_prefix.as_ptr(), 0);
        assert!(!idx.is_null(), "index load failed");
        (dir, refseq, idx)
    }
}

#[test]
fn build_load_map_one_read() {
    let (dir, refseq, idx) = build_and_load_index("smoke");

    unsafe {
        let mut opt: minibwa_sys::mb_opt_t = std::mem::zeroed();
        minibwa_sys::mb_opt_init(&mut opt);

        // Take a 100 bp substring of the reference as the query.
        let query = &refseq[500..600];
        let c_seq = CString::new(query).unwrap();
        let c_name = CString::new("q1").unwrap();
        let mut n_hit: i32 = 0;
        let hits = minibwa_sys::mb_map(
            &opt,
            idx,
            query.len() as i32,
            c_seq.as_ptr(),
            0,
            &mut n_hit,
            std::ptr::null_mut(),
            c_name.as_ptr(),
        );
        assert!(n_hit >= 1, "expected at least one hit");
        assert!(!hits.is_null());

        let h = &*hits;
        let name = CStr::from_ptr(minibwa_sys::mb_idx_ctg_name(idx, h.tid as i32));
        assert_eq!(name.to_str().unwrap(), "chr1");

        // Free each .p then the array (libc-allocated).
        for i in 0..n_hit as isize {
            let hp = (*hits.offset(i)).p;
            if !hp.is_null() {
                libc::free(hp as *mut libc::c_void);
            }
        }
        libc::free(hits as *mut libc::c_void);
        minibwa_sys::mb_idx_destroy(idx);
    }
    std::fs::remove_dir_all(&dir).ok();
}

/// The oracle entry points that size a raw allocation from a caller-supplied count must
/// treat an empty batch as an empty result and reject a negative count outright, rather
/// than sizing a `malloc` from it. `mb_lchain_dp` itself contracts on `n == 0 || a == 0`,
/// so the empty case has a defined answer: no chains, no anchors.
#[test]
fn shim_count_boundaries() {
    let (dir, _refseq, idx) = build_and_load_index("bounds");

    unsafe {
        let mut out_u = [0u64; 4];
        let mut out_sid = [0i32; 4];
        let mut out_len = [0i32; 4];
        let mut out_qpos = [0i32; 4];
        let mut out_tpos = [0i64; 4];
        let mut n_a: std::os::raw::c_int = -7; // poisoned: the callee must overwrite it

        // Empty anchor set: no chains, no anchors, nothing allocated.
        let n_u = minibwa_sys::mb_shim_lchain_dp(
            idx,
            5000,
            5000,
            500,
            25,
            5000,
            40,
            0.12,
            0,
            std::ptr::null(),
            std::ptr::null(),
            std::ptr::null(),
            std::ptr::null(),
            out_u.as_mut_ptr(),
            out_u.len() as i32,
            out_sid.as_mut_ptr(),
            out_len.as_mut_ptr(),
            out_qpos.as_mut_ptr(),
            out_tpos.as_mut_ptr(),
            out_sid.len() as i32,
            &mut n_a,
        );
        assert_eq!(n_u, 0, "empty anchor set should yield no chains");
        assert_eq!(n_a, 0, "empty anchor set should yield no anchors");

        // Negative count is rejected before any allocation is sized from it.
        let rc = minibwa_sys::mb_shim_lchain_dp(
            idx,
            5000,
            5000,
            500,
            25,
            5000,
            40,
            0.12,
            -1,
            std::ptr::null(),
            std::ptr::null(),
            std::ptr::null(),
            std::ptr::null(),
            out_u.as_mut_ptr(),
            out_u.len() as i32,
            out_sid.as_mut_ptr(),
            out_len.as_mut_ptr(),
            out_qpos.as_mut_ptr(),
            out_tpos.as_mut_ptr(),
            out_sid.len() as i32,
            &mut n_a,
        );
        assert_eq!(rc, -1, "negative anchor count should be rejected");

        // Same contract on the seeding side of the shim.
        let qseq = [0u8; 32];
        let n_anchor = minibwa_sys::mb_shim_anchor(
            idx,
            0,
            std::ptr::null(),
            std::ptr::null(),
            std::ptr::null(),
            19,
            qseq.len() as i32,
            qseq.as_ptr(),
            500,
            out_sid.as_mut_ptr(),
            out_len.as_mut_ptr(),
            out_qpos.as_mut_ptr(),
            out_tpos.as_mut_ptr(),
            out_sid.len() as i32,
        );
        assert_eq!(n_anchor, 0, "empty seed set should yield no anchors");

        let rc = minibwa_sys::mb_shim_anchor(
            idx,
            -1,
            std::ptr::null(),
            std::ptr::null(),
            std::ptr::null(),
            19,
            qseq.len() as i32,
            qseq.as_ptr(),
            500,
            out_sid.as_mut_ptr(),
            out_len.as_mut_ptr(),
            out_qpos.as_mut_ptr(),
            out_tpos.as_mut_ptr(),
            out_sid.len() as i32,
        );
        assert_eq!(rc, -1, "negative seed count should be rejected");

        minibwa_sys::mb_idx_destroy(idx);
    }
    std::fs::remove_dir_all(&dir).ok();
}
