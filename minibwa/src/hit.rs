use crate::index::Index;
use minibwa_sys as sys;

/// Alignment strand.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
}

/// A CIGAR operation kind, mirroring minibwa's `MB_CIGAR_*`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CigarKind {
    /// Alignment match or mismatch (M).
    Match,
    /// Insertion to the reference (I).
    Ins,
    /// Deletion from the reference (D).
    Del,
    /// Skipped region from the reference (N).
    RefSkip,
    /// Soft clip — bases present in SEQ (S).
    SoftClip,
    /// Hard clip — bases absent from SEQ (H).
    HardClip,
    /// Padding (silent deletion from padded reference) (P).
    Pad,
    /// Sequence match (=).
    Eq,
    /// Sequence mismatch (X).
    Diff,
}

impl CigarKind {
    /// Convert the 4-bit op code from a minibwa CIGAR word into a [`CigarKind`].
    ///
    /// # Panics
    ///
    /// Panics if `op > 8`; this is unreachable for well-formed minibwa output.
    pub(crate) fn from_op(op: u32) -> CigarKind {
        match op {
            0 => CigarKind::Match,
            1 => CigarKind::Ins,
            2 => CigarKind::Del,
            3 => CigarKind::RefSkip,
            4 => CigarKind::SoftClip,
            5 => CigarKind::HardClip,
            6 => CigarKind::Pad,
            7 => CigarKind::Eq,
            8 => CigarKind::Diff,
            other => panic!("unknown minibwa cigar op {other}"),
        }
    }
}

/// A single CIGAR operation: a kind and a length.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CigarOp {
    /// The operation kind (M, I, D, …).
    pub kind: CigarKind,
    /// Number of bases consumed by this operation.
    pub len: u32,
}

impl std::fmt::Display for CigarOp {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let ch = match self.kind {
            CigarKind::Match => 'M',
            CigarKind::Ins => 'I',
            CigarKind::Del => 'D',
            CigarKind::RefSkip => 'N',
            CigarKind::SoftClip => 'S',
            CigarKind::HardClip => 'H',
            CigarKind::Pad => 'P',
            CigarKind::Eq => '=',
            CigarKind::Diff => 'X',
        };
        write!(f, "{}{}", self.len, ch)
    }
}

/// A fully-owned alignment record converted from minibwa's `mb_hit_t`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Hit {
    /// 0-based reference sequence index (target id).
    pub tid: i64,
    /// Reference contig name, or `None` for unmapped hits.
    pub contig: Option<String>,
    /// 0-based start position on the reference (inclusive).
    pub ref_start: i64,
    /// 0-based end position on the reference (exclusive).
    pub ref_end: i64,
    /// 0-based start position on the query (inclusive).
    pub query_start: i32,
    /// 0-based end position on the query (exclusive).
    pub query_end: i32,
    /// Alignment strand.
    pub strand: Strand,
    /// Mapping quality (0–255).
    pub mapq: u8,
    /// Alignment score.
    pub score: i32,
    /// Number of suboptimal hits.
    pub n_sub: i32,
    /// True when this hit is part of a proper pair.
    pub proper_pair: bool,
    /// True when this is the primary alignment.
    pub is_primary: bool,
    /// True when this is a secondary alignment (SAM flag 0x100).
    pub is_secondary: bool,
    /// True when this is a supplementary alignment (SAM flag 0x800).
    pub is_supplementary: bool,
    /// CIGAR over the aligned region only (no leading/trailing clips).
    pub cigar: Vec<CigarOp>,
}

impl Hit {
    /// Return the CIGAR as a SAM-style string (e.g. `"100M5I45M"`).
    ///
    /// Returns `"*"` if the CIGAR is empty (unmapped).
    pub fn cigar_string(&self) -> String {
        if self.cigar.is_empty() {
            return "*".to_owned();
        }
        self.cigar.iter().map(|op| op.to_string()).collect()
    }
}

/// Classify a hit's SAM-flag status from its `parent`/`id`/`sam_pri` fields,
/// matching minibwa's own SAM writer (`format.c`):
///   secondary     = parent != id             (SAM 0x100)
///   supplementary = parent == id && !sam_pri  (SAM 0x800)
///   primary       = parent == id && sam_pri
/// Returns `(is_primary, is_secondary, is_supplementary)`.
pub(crate) fn classify(parent: i32, id: i32, sam_pri: bool) -> (bool, bool, bool) {
    if parent != id {
        (false, true, false) // secondary
    } else if !sam_pri {
        (false, false, true) // supplementary
    } else {
        (true, false, false) // primary
    }
}

/// Convert a raw `mb_hit_t` into an owned `Hit`. Does NOT free `h.p`.
///
/// # Safety
/// `h` must be a valid `mb_hit_t` from `mb_map`/`mb_map_batch`, and `h.p`, if
/// non-null, must point to a valid `mb_extra_t` with `n_cigar` cigar words.
pub(crate) unsafe fn hit_from_raw(h: &sys::mb_hit_t, idx: &Index) -> Hit {
    let mut cigar = Vec::new();
    if !h.p.is_null() {
        // SAFETY: caller guarantees h.p is a valid mb_extra_t.
        let extra = unsafe { &*h.p };
        let n = extra.n_cigar.max(0) as usize;
        // Read exactly n_cigar words; cap includes appended cs/MD bytes.
        // SAFETY: cigar[] has n_cigar valid u32 words.
        let words = unsafe { extra.cigar.as_slice(n) };
        cigar.reserve(n);
        for &w in words {
            cigar.push(CigarOp {
                kind: CigarKind::from_op(w & 0xf),
                len: w >> 4,
            });
        }
    }
    let (is_primary, is_secondary, is_supplementary) = classify(h.parent, h.id, h.sam_pri() != 0);
    Hit {
        tid: h.tid,
        contig: idx.contig_name(h.tid as usize).map(str::to_owned),
        ref_start: h.ts,
        ref_end: h.te,
        query_start: h.qs,
        query_end: h.qe,
        strand: if h.rev() != 0 {
            Strand::Reverse
        } else {
            Strand::Forward
        },
        mapq: h.mapq.clamp(0, 255) as u8,
        score: h.score,
        n_sub: h.n_sub,
        proper_pair: h.proper_pair() != 0,
        is_primary,
        is_secondary,
        is_supplementary,
        cigar,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::meth::Meth;

    #[test]
    fn cigar_kind_from_op() {
        assert_eq!(CigarKind::from_op(0), CigarKind::Match);
        assert_eq!(CigarKind::from_op(4), CigarKind::SoftClip);
        assert_eq!(CigarKind::from_op(8), CigarKind::Diff);
    }

    #[test]
    fn meth_codes() {
        assert_eq!(Meth::None.as_mt(), 0);
        assert_eq!(Meth::C2T.as_mt(), 1);
        assert_eq!(Meth::G2A.as_mt(), 2);
    }

    #[test]
    fn classify_flags() {
        assert_eq!(classify(5, 5, true), (true, false, false)); // primary
        assert_eq!(classify(5, 5, false), (false, false, true)); // supplementary
        assert_eq!(classify(2, 5, true), (false, true, false)); // secondary (parent != id)
    }

    #[test]
    fn cigar_op_display() {
        assert_eq!(
            CigarOp {
                kind: CigarKind::Match,
                len: 100
            }
            .to_string(),
            "100M"
        );
        assert_eq!(
            CigarOp {
                kind: CigarKind::Ins,
                len: 5
            }
            .to_string(),
            "5I"
        );
        assert_eq!(
            CigarOp {
                kind: CigarKind::Del,
                len: 3
            }
            .to_string(),
            "3D"
        );
        assert_eq!(
            CigarOp {
                kind: CigarKind::Eq,
                len: 10
            }
            .to_string(),
            "10="
        );
        assert_eq!(
            CigarOp {
                kind: CigarKind::Diff,
                len: 2
            }
            .to_string(),
            "2X"
        );
        assert_eq!(
            CigarOp {
                kind: CigarKind::SoftClip,
                len: 7
            }
            .to_string(),
            "7S"
        );
    }

    #[test]
    fn cigar_string_empty_is_star() {
        let hit = Hit {
            tid: -1,
            contig: None,
            ref_start: 0,
            ref_end: 0,
            query_start: 0,
            query_end: 0,
            strand: Strand::Forward,
            mapq: 0,
            score: 0,
            n_sub: 0,
            proper_pair: false,
            is_primary: false,
            is_secondary: false,
            is_supplementary: false,
            cigar: vec![],
        };
        assert_eq!(hit.cigar_string(), "*");
    }

    #[test]
    fn cigar_string_nonempty() {
        let hit = Hit {
            tid: 0,
            contig: Some("chr1".to_owned()),
            ref_start: 100,
            ref_end: 250,
            query_start: 0,
            query_end: 150,
            strand: Strand::Forward,
            mapq: 60,
            score: 100,
            n_sub: 0,
            proper_pair: false,
            is_primary: true,
            is_secondary: false,
            is_supplementary: false,
            cigar: vec![
                CigarOp {
                    kind: CigarKind::Match,
                    len: 100,
                },
                CigarOp {
                    kind: CigarKind::Ins,
                    len: 5,
                },
                CigarOp {
                    kind: CigarKind::Match,
                    len: 45,
                },
            ],
        };
        assert_eq!(hit.cigar_string(), "100M5I45M");
    }
}
