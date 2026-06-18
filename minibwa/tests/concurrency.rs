//! N threads share one Index; each owns its ThreadBuf and aligns reads.
use minibwa::test_util::{synthetic_reference, write_fasta};
use minibwa::{Aligner, Index, Meth, Opts, ThreadBuf};
use std::sync::Arc;

#[test]
fn shared_index_across_threads() {
    let dir = std::env::temp_dir().join(format!("minibwa_conc_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let seq = synthetic_reference(20_000, 99);
    let fa = write_fasta(&dir, "chr1", &seq);
    let prefix = dir.join("idx");
    Index::build_from_fasta(&fa, &prefix, false, 1).unwrap();
    let idx = Arc::new(Index::load(&prefix, false).unwrap());

    let mut handles = Vec::new();
    for t in 0..4u32 {
        let idx = Arc::clone(&idx);
        let seq = seq.clone();
        handles.push(std::thread::spawn(move || {
            let opts = Opts::new();
            let aligner = Aligner::new(&idx, &opts);
            let mut buf = ThreadBuf::new();
            let mut mapped = 0;
            for i in 0..50 {
                let start = ((t as usize * 137 + i * 311) % (seq.len() - 200)) as usize;
                let q = &seq[start..start + 150];
                let hits = aligner
                    .map(&mut buf, format!("t{t}_r{i}").as_bytes(), q, Meth::None)
                    .unwrap();
                if !hits.is_empty() {
                    mapped += 1;
                }
            }
            mapped
        }));
    }
    let total: i32 = handles.into_iter().map(|h| h.join().unwrap()).sum();
    assert!(total > 150, "expected most reads to map, got {total}");
    std::fs::remove_dir_all(&dir).ok();
}
