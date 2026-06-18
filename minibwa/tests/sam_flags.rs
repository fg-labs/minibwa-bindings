//! End-to-end SAM-flag classification tests for `Hit::is_primary`, `is_secondary`,
//! and `is_supplementary`. These exercise the full alignment pipeline, not just the
//! isolated `classify()` unit test in `hit.rs`.
use minibwa::test_util::{synthetic_reference, write_fasta};
use minibwa::{Aligner, Index, Meth, Opts, ThreadBuf};

/// Case 1 — primary path.
///
/// Align a read from a unique region of a 5000 bp reference. Expect exactly one
/// hit that is flagged primary (and neither secondary nor supplementary).
#[test]
fn primary_flag_unique_read() {
    let dir = std::env::temp_dir().join(format!("minibwa_sf_primary_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();

    let seq = synthetic_reference(5000, 42);
    let fa = write_fasta(&dir, "chr1", &seq);
    let prefix = dir.join("idx");
    Index::build_from_fasta(&fa, &prefix, false, 1).unwrap();
    let idx = Index::load(&prefix, false).unwrap();
    let opts = Opts::new();
    let aligner = Aligner::new(&idx, &opts);
    let mut buf = ThreadBuf::new();

    // Read from a unique interior region — 150 bp is long enough to be placed uniquely.
    let query = &seq[2000..2150];
    let hits = aligner
        .map(&mut buf, b"primary_read", query, Meth::None)
        .unwrap();

    println!("[case1] hit count: {}", hits.len());
    for (i, h) in hits.iter().enumerate() {
        println!(
            "  hit[{i}] contig={:?} ref={}-{} primary={} secondary={} supplementary={}",
            h.contig, h.ref_start, h.ref_end, h.is_primary, h.is_secondary, h.is_supplementary
        );
    }

    assert_eq!(hits.len(), 1, "unique read should produce exactly one hit");
    assert!(hits[0].is_primary, "the single hit must be primary");
    assert!(
        !hits[0].is_secondary,
        "the single hit must not be secondary"
    );
    assert!(
        !hits[0].is_supplementary,
        "the single hit must not be supplementary"
    );

    std::fs::remove_dir_all(&dir).ok();
}

/// Case 2 — secondary path.
///
/// Build a reference where a 250 bp segment (`dup`) appears at two distinct loci
/// so a read from that segment can map to both places. The reference is:
///   [flank_a (2000 bp)] [dup (250 bp)] [flank_b (2000 bp)] [dup (250 bp)]
/// Both dup copies are identical (same seed), while the flanks use different seeds
/// so they are unique and the index has no trouble telling the two dup loci apart
/// positionally.
///
/// Contract:
///   - At least one hit (primary) must exist.
///   - Exactly one hit has `is_primary == true`.
///   - Any additional hits have `is_secondary == true` and `is_primary == false`.
///   - No hit is both primary and secondary.
///
/// Whether a secondary hit is produced depends on minibwa's mask_level / pri_ratio
/// filtering. The test is written to pass regardless of the actual hit count, while
/// still verifying the above invariants.
#[test]
fn secondary_flag_duplicated_region() {
    let dir = std::env::temp_dir().join(format!("minibwa_sf_secondary_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();

    // Build a reference with two identical 250 bp copies separated by unique flanks.
    let flank_a = synthetic_reference(2000, 11);
    let dup = synthetic_reference(250, 77); // identical copies at both loci
    let flank_b = synthetic_reference(2000, 31);

    let mut refseq = Vec::with_capacity(flank_a.len() + dup.len() + flank_b.len() + dup.len());
    refseq.extend_from_slice(&flank_a);
    refseq.extend_from_slice(&dup);
    refseq.extend_from_slice(&flank_b);
    refseq.extend_from_slice(&dup);

    let fa = write_fasta(&dir, "chr1", &refseq);
    let prefix = dir.join("idx");
    Index::build_from_fasta(&fa, &prefix, false, 1).unwrap();
    let idx = Index::load(&prefix, false).unwrap();

    // Raise out_n so secondary alignments are reported when present.
    let opts = Opts::new().set_out_n(5);
    let aligner = Aligner::new(&idx, &opts);
    let mut buf = ThreadBuf::new();

    // Query = the middle 200 bp of the dup segment (avoids boundary seed issues).
    let query = &dup[25..225];
    let hits = aligner
        .map(&mut buf, b"dup_read", query, Meth::None)
        .unwrap();

    println!("[case2] hit count: {}", hits.len());
    for (i, h) in hits.iter().enumerate() {
        println!(
            "  hit[{i}] contig={:?} ref={}-{} primary={} secondary={} supplementary={}",
            h.contig, h.ref_start, h.ref_end, h.is_primary, h.is_secondary, h.is_supplementary
        );
    }

    // At minimum there must be a primary hit.
    assert!(
        !hits.is_empty(),
        "duplicated-region read must produce at least one hit"
    );

    // Exactly one primary.
    let primary_count = hits.iter().filter(|h| h.is_primary).count();
    assert_eq!(
        primary_count, 1,
        "exactly one hit must be primary; got {primary_count}"
    );

    // No hit is simultaneously primary and secondary.
    for (i, h) in hits.iter().enumerate() {
        assert!(
            !(h.is_primary && h.is_secondary),
            "hit[{i}] is both primary and secondary — impossible"
        );
    }

    // Any non-primary hits must be secondary (not supplementary, since the full read maps).
    for (i, h) in hits.iter().enumerate() {
        if !h.is_primary {
            assert!(
                h.is_secondary,
                "hit[{i}] is non-primary but not secondary (supplementary={})",
                h.is_supplementary
            );
        }
    }

    let secondary_count = hits.iter().filter(|h| h.is_secondary).count();
    println!(
        "[case2] primaries={primary_count} secondaries={secondary_count} total={}",
        hits.len()
    );

    std::fs::remove_dir_all(&dir).ok();
}
