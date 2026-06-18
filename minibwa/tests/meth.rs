//! A meth index (superset) loads both ways; C2T/G2A reads align.
use minibwa::test_util::{synthetic_reference, write_fasta};
use minibwa::{Aligner, Index, Meth, Opts, ThreadBuf};

#[test]
fn meth_index_builds_and_loads_both_ways() {
    let dir = std::env::temp_dir().join(format!("minibwa_meth_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let seq = synthetic_reference(10_000, 5);
    let fa = write_fasta(&dir, "chr1", &seq);
    let prefix = dir.join("idx");

    // meth=true writes both normal and meth indexes.
    Index::build_from_fasta(&fa, &prefix, true, 1).unwrap();

    // Loads normally (uses .mbw).
    let plain = Index::load(&prefix, false).unwrap();
    assert_eq!(plain.contig_name(0), Some("chr1"));

    // Loads as meth (uses .meth.mbw).
    let meth_idx = Index::load(&prefix, true).unwrap();
    let opts = Opts::new().set_methylation(true);
    let aligner = Aligner::new(&meth_idx, &opts);
    let mut buf = ThreadBuf::new();

    // A C2T-converted copy of a reference substring should still map under Meth::C2T.
    let mut q: Vec<u8> = seq[2000..2150].to_vec();
    for b in q.iter_mut() {
        if *b == b'C' {
            *b = b'T';
        }
    }
    let hits = aligner.map(&mut buf, b"m1", &q, Meth::C2T).unwrap();
    assert!(!hits.is_empty(), "expected a methylation hit");
    std::fs::remove_dir_all(&dir).ok();
}
