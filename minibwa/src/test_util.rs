//! Deterministic synthetic test data generators (no committed fixtures).
#![cfg(any(test, feature = "test-util"))]

use std::io::Write;
use std::path::{Path, PathBuf};

/// Deterministic pseudo-random ACGT reference of `len` bases.
pub fn synthetic_reference(len: usize, seed: u32) -> Vec<u8> {
    let mut x = seed | 1;
    let mut v = Vec::with_capacity(len);
    for _ in 0..len {
        x = x.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
        v.push(b"ACGT"[(x >> 16) as usize & 3]);
    }
    v
}

/// Write a single-contig FASTA; returns the path.
pub fn write_fasta(dir: &Path, contig: &str, seq: &[u8]) -> PathBuf {
    let fa = dir.join(format!("{contig}.fa"));
    let mut f = std::fs::File::create(&fa).unwrap();
    writeln!(f, ">{contig}").unwrap();
    f.write_all(seq).unwrap();
    writeln!(f).unwrap();
    fa
}

/// Write a multi-contig FASTA; returns the path.
///
/// Each entry in `contigs` is `(name, sequence)`.  The output file is named
/// `multi.fa` inside `dir`.
pub fn write_multi_fasta(dir: &Path, contigs: &[(&str, &[u8])]) -> PathBuf {
    let fa = dir.join("multi.fa");
    let mut f = std::fs::File::create(&fa).unwrap();
    for (name, seq) in contigs {
        writeln!(f, ">{name}").unwrap();
        f.write_all(seq).unwrap();
        writeln!(f).unwrap();
    }
    fa
}
