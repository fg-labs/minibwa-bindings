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

#[test]
fn build_load_map_one_read() {
    let dir = std::env::temp_dir().join(format!("minibwa_sys_smoke_{}", std::process::id()));
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
