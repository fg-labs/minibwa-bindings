// pyo3 0.22 generates Into<PyErr> conversion shims in #[pymethods]/#[pyfunction]
// that clippy::useless_conversion flags; the lint fires in macro-generated code
// that refers back to the original function signature spans.
#![allow(clippy::useless_conversion)]

use minibwa::{Aligner, Index as RsIndex, Meth, Opts as RsOpts, ThreadBuf};
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use std::sync::Arc;

fn map_err(e: minibwa::Error) -> PyErr {
    match e {
        minibwa::Error::InvalidInput(m) | minibwa::Error::InvalidOpts(m) => {
            PyValueError::new_err(m)
        }
        other => PyRuntimeError::new_err(other.to_string()),
    }
}

/// A minibwa reference index.
///
/// Build an index once with ``Index.build``, then load it with
/// ``Index.load`` before aligning reads.  The loaded index is
/// reference-counted internally so it can be shared cheaply across
/// aligners.
#[pyclass(name = "Index")]
struct PyIndex {
    inner: Arc<RsIndex>,
}

#[pymethods]
impl PyIndex {
    /// Build a minibwa index from a FASTA file.
    ///
    /// Writes ``<prefix>.mbw`` and ``<prefix>.l2b`` to disk.
    ///
    /// Args:
    ///     fasta:   Path to the reference FASTA (may be gzipped).
    ///     prefix:  Output file prefix (directory must exist).
    ///     meth:    Build a bisulfite/EM-seq methylation index (default
    ///              ``False``).
    ///     threads: Number of threads for index construction (default 1).
    ///
    /// Raises:
    ///     ValueError: On invalid input.
    ///     RuntimeError: On I/O or internal errors.
    #[staticmethod]
    #[pyo3(signature = (fasta, prefix, meth=false, threads=1))]
    fn build(fasta: &str, prefix: &str, meth: bool, threads: u32) -> PyResult<()> {
        RsIndex::build_from_fasta(fasta, prefix, meth, threads).map_err(map_err)
    }

    /// Load a previously built minibwa index from disk.
    ///
    /// Args:
    ///     prefix: File prefix used when the index was built.
    ///     meth:   Load as a methylation index (must match build setting,
    ///             default ``False``).
    ///
    /// Returns:
    ///     A ready-to-use ``Index`` object.
    ///
    /// Raises:
    ///     RuntimeError: If the index files cannot be read.
    #[staticmethod]
    #[pyo3(signature = (prefix, meth=false))]
    fn load(prefix: &str, meth: bool) -> PyResult<Self> {
        Ok(PyIndex {
            inner: Arc::new(RsIndex::load(prefix, meth).map_err(map_err)?),
        })
    }

    /// Return the number of contigs (sequences) in the index.
    fn n_contigs(&self) -> usize {
        self.inner.n_contigs()
    }

    /// Return a list of ``(name, length)`` tuples for every contig.
    ///
    /// Returns:
    ///     ``list[tuple[str, int]]`` in the order they appear in the index.
    fn contigs(&self) -> Vec<(String, i64)> {
        self.inner
            .contigs()
            .map(|(n, l)| (n.to_string(), l))
            .collect()
    }
}

/// Alignment options for minibwa.
///
/// Create with ``Opts()`` for defaults (``"adap"`` preset) or
/// ``Opts(preset="sr")`` for short-read mode.  Individual parameters can
/// be adjusted with the ``set_*`` methods after construction.
#[pyclass(name = "Opts")]
struct PyOpts {
    inner: RsOpts,
}

#[pymethods]
impl PyOpts {
    /// Create alignment options, optionally from a named preset.
    ///
    /// Args:
    ///     preset: One of ``"sr"`` (short-read), ``"adap"`` (adaptive), or
    ///             ``"lr"`` (long-read). ``None`` (default) uses the ``"adap"``
    ///             defaults from minibwa.
    ///
    /// Raises:
    ///     ValueError: If the preset name is unknown.
    #[new]
    #[pyo3(signature = (preset=None))]
    fn new(preset: Option<&str>) -> PyResult<Self> {
        let inner = match preset {
            Some(p) => RsOpts::with_preset(p).map_err(map_err)?,
            None => RsOpts::new(),
        };
        Ok(PyOpts { inner })
    }

    /// Enable or disable paired-end mode.
    ///
    /// Args:
    ///     on: ``True`` to enable PE mode (``MB_F_PE`` flag).
    fn set_paired(&mut self, on: bool) {
        // Opts is a consuming builder; clone, transform, replace.
        self.inner = self.inner.clone().set_paired(on);
    }

    /// Enable or disable bisulfite/EM-seq methylation mode.
    ///
    /// Args:
    ///     on: ``True`` to enable methylation mode (``MB_F_METH`` flag).
    fn set_methylation(&mut self, on: bool) {
        self.inner = self.inner.clone().set_methylation(on);
    }

    /// Set the minimum seed length for BWT seeding.
    ///
    /// Args:
    ///     v: Minimum seed length (default depends on preset).
    fn set_min_seed_len(&mut self, v: i32) {
        self.inner = self.inner.clone().set_min_seed_len(v);
    }

    /// Set the maximum number of secondary alignments to output.
    ///
    /// Args:
    ///     v: Max secondary alignments per read.
    fn set_out_n(&mut self, v: i32) {
        self.inner = self.inner.clone().set_out_n(v);
    }

    /// Set the Smith-Waterman match score.
    ///
    /// Args:
    ///     score: Match score (``a`` parameter, positive integer).
    fn set_match_score(&mut self, score: i32) {
        self.inner = self.inner.clone().set_match_score(score);
    }

    /// Set the Smith-Waterman mismatch penalty.
    ///
    /// Args:
    ///     penalty: Mismatch penalty (``b`` parameter, positive integer).
    fn set_mismatch_penalty(&mut self, penalty: i32) {
        self.inner = self.inner.clone().set_mismatch_penalty(penalty);
    }

    /// Set the gap-open penalty.
    ///
    /// Args:
    ///     open: Gap-open penalty (``q`` parameter, positive integer).
    fn set_gap_open(&mut self, open: i32) {
        self.inner = self.inner.clone().set_gap_open(open);
    }

    /// Set the gap-extend penalty.
    ///
    /// Args:
    ///     extend: Gap-extend penalty (``e`` parameter, positive integer).
    fn set_gap_extend(&mut self, extend: i32) {
        self.inner = self.inner.clone().set_gap_extend(extend);
    }

    /// Set paired-end insert-size parameters for proper-pair classification.
    ///
    /// Args:
    ///     avg: Expected insert-size mean.
    ///     std: Expected insert-size standard deviation.
    ///     lo:  Lower insert-size bound for proper-pair classification.
    ///     hi:  Upper insert-size bound for proper-pair classification.
    fn set_pe_insert_size(&mut self, avg: i32, std: i32, lo: i32, hi: i32) {
        self.inner = self.inner.clone().set_pe_insert_size(avg, std, lo, hi);
    }
}

/// One alignment result returned by ``map``, ``map_pair``, or ``map_many``.
///
/// Attributes:
///     contig (str | None): Reference contig name, or ``None`` for unmapped.
///     ref_start (int):     0-based start on the reference (inclusive).
///     ref_end (int):       0-based end on the reference (exclusive).
///     query_start (int):   0-based start on the query (inclusive).
///     query_end (int):     0-based end on the query (exclusive).
///     reverse (bool):      ``True`` if aligned to the reverse strand.
///                          Prefer ``strand`` for new code.
///     strand (str):        ``"+"`` (forward) or ``"-"`` (reverse).
///     mapq (int):          Mapping quality (0–255).
///     score (int):         Alignment score.
///     n_sub (int):         Number of sub-optimal hits.
///     proper_pair (bool):  ``True`` if insert size is within the expected
///                          range (only meaningful for paired-end hits).
///     is_primary (bool):   ``True`` for the primary alignment.
///     is_secondary (bool): ``True`` for secondary alignments.
///     is_supplementary (bool): ``True`` for supplementary alignments.
///     cigar (list[tuple[str, int]]): CIGAR operations as (op_char, length)
///                          pairs, e.g. ``[('M', 150)]``.
#[pyclass(name = "Hit")]
struct PyHit {
    #[pyo3(get)]
    contig: Option<String>,
    #[pyo3(get)]
    ref_start: i64,
    #[pyo3(get)]
    ref_end: i64,
    #[pyo3(get)]
    query_start: i32,
    #[pyo3(get)]
    query_end: i32,
    #[pyo3(get)]
    reverse: bool,
    #[pyo3(get)]
    mapq: u8,
    #[pyo3(get)]
    score: i32,
    #[pyo3(get)]
    n_sub: i32,
    #[pyo3(get)]
    proper_pair: bool,
    #[pyo3(get)]
    is_primary: bool,
    #[pyo3(get)]
    is_secondary: bool,
    #[pyo3(get)]
    is_supplementary: bool,
    /// CIGAR as a list of (op_char, len).
    #[pyo3(get)]
    cigar: Vec<(char, u32)>,
}

#[pymethods]
impl PyHit {
    /// Return ``"+"`` for a forward-strand hit or ``"-"`` for reverse.
    ///
    /// The ``reverse`` field is kept for backward compatibility; ``strand``
    /// is the preferred accessor.
    #[getter]
    fn strand(&self) -> &'static str {
        if self.reverse { "-" } else { "+" }
    }

    fn __repr__(&self) -> String {
        format!(
            "Hit(contig={:?}, ref_start={}, ref_end={}, mapq={}, score={}, strand={}, cigar_ops={})",
            self.contig,
            self.ref_start,
            self.ref_end,
            self.mapq,
            self.score,
            if self.reverse { "-" } else { "+" },
            self.cigar.len()
        )
    }
}

fn cigar_char(k: minibwa::CigarKind) -> char {
    use minibwa::CigarKind::*;
    match k {
        Match => 'M',
        Ins => 'I',
        Del => 'D',
        RefSkip => 'N',
        SoftClip => 'S',
        HardClip => 'H',
        Pad => 'P',
        Eq => '=',
        Diff => 'X',
    }
}

fn to_pyhit(h: minibwa::Hit) -> PyHit {
    PyHit {
        contig: h.contig,
        ref_start: h.ref_start,
        ref_end: h.ref_end,
        query_start: h.query_start,
        query_end: h.query_end,
        reverse: matches!(h.strand, minibwa::Strand::Reverse),
        mapq: h.mapq,
        score: h.score,
        n_sub: h.n_sub,
        proper_pair: h.proper_pair,
        is_primary: h.is_primary,
        is_secondary: h.is_secondary,
        is_supplementary: h.is_supplementary,
        cigar: h.cigar.iter().map(|c| (cigar_char(c.kind), c.len)).collect(),
    }
}

/// Align one read against the index.
///
/// The GIL is released during alignment so the call does not block other
/// Python threads.
///
/// Args:
///     index: A loaded ``Index``.
///     opts:  Alignment options.
///     name:  Read name (must not contain NUL bytes).
///     seq:   Read sequence (non-empty DNA string).
///     meth:  Methylation strand — ``"none"`` (default), ``"c2t"``, or
///            ``"g2a"``.  ``Meth`` enum values are accepted too.
///
/// Returns:
///     ``list[Hit]`` — may be empty if the read does not align.
///
/// Raises:
///     ValueError: If ``seq`` is empty, ``meth`` is unrecognised, or
///                 ``name`` contains a NUL byte.
#[pyfunction]
#[pyo3(signature = (index, opts, name, seq, meth="none"))]
fn map(
    py: Python<'_>,
    index: &PyIndex,
    opts: &PyOpts,
    name: &str,
    seq: &str,
    meth: &str,
) -> PyResult<Vec<PyHit>> {
    let mt = match meth {
        "none" => Meth::None,
        "c2t" => Meth::C2T,
        "g2a" => Meth::G2A,
        other => return Err(PyValueError::new_err(format!("bad meth {other:?}"))),
    };
    let idx = Arc::clone(&index.inner);
    let name = name.as_bytes().to_vec();
    let seq = seq.as_bytes().to_vec();
    // Release the GIL during alignment.
    let hits = py.allow_threads(|| {
        let aligner = Aligner::new(&idx, &opts.inner);
        let mut buf = ThreadBuf::new();
        aligner.map(&mut buf, &name, &seq, mt)
    });
    let hits = hits.map_err(map_err)?;
    Ok(hits.into_iter().map(to_pyhit).collect())
}

/// Align a read pair (paired-end) against the index.
///
/// Enables paired-end scoring, mate rescue, and proper-pair flagging via
/// minibwa's PE path.  If the options enable methylation, R1 is treated as
/// C2T and R2 as G2A automatically.  The GIL is released during alignment.
///
/// Args:
///     index: A loaded ``Index``.
///     opts:  Alignment options (insert-size parameters are read from here).
///     name1: R1 read name.
///     seq1:  R1 read sequence (non-empty).
///     name2: R2 read name.
///     seq2:  R2 read sequence (non-empty).
///
/// Returns:
///     ``(list[Hit], list[Hit])`` — hits for R1 and R2 respectively.
///
/// Raises:
///     ValueError: If either sequence is empty or either name contains a
///                 NUL byte.
#[pyfunction]
#[pyo3(signature = (index, opts, name1, seq1, name2, seq2))]
fn map_pair(
    py: Python<'_>,
    index: &PyIndex,
    opts: &PyOpts,
    name1: &str,
    seq1: &str,
    name2: &str,
    seq2: &str,
) -> PyResult<(Vec<PyHit>, Vec<PyHit>)> {
    let idx = Arc::clone(&index.inner);
    let (n1, s1, n2, s2) = (
        name1.as_bytes().to_vec(),
        seq1.as_bytes().to_vec(),
        name2.as_bytes().to_vec(),
        seq2.as_bytes().to_vec(),
    );
    let (h1, h2) = py
        .allow_threads(|| {
            let aligner = Aligner::new(&idx, &opts.inner);
            let mut buf = ThreadBuf::new();
            aligner.map_pair(&mut buf, &n1, &s1, &n2, &s2)
        })
        .map_err(map_err)?;
    Ok((
        h1.into_iter().map(to_pyhit).collect(),
        h2.into_iter().map(to_pyhit).collect(),
    ))
}

/// Align many independent reads in a single batched call.
///
/// This is the high-throughput entry point: minibwa batches BWT seeding
/// across the whole slice (prefetch-driven), which the single-read ``map``
/// cannot do.  Paired-end mode is forced off — reads are aligned
/// independently.  If the options enable methylation, every read is treated
/// as the C2T (read-1) strand; use ``map`` per read for G2A.
///
/// ``names`` and ``seqs`` must be the same length.
///
/// Args:
///     index: A loaded ``Index``.
///     opts:  Alignment options.
///     names: Read names (one per read).
///     seqs:  Read sequences (one per read, same order as ``names``).
///
/// Returns:
///     ``list[list[Hit]]`` — one inner list per input read, in order.
///
/// Raises:
///     ValueError: If ``names`` and ``seqs`` differ in length, any sequence
///                 is empty, or any name contains a NUL byte.
#[pyfunction]
#[pyo3(signature = (index, opts, names, seqs))]
fn map_many(
    py: Python<'_>,
    index: &PyIndex,
    opts: &PyOpts,
    names: Vec<String>,
    seqs: Vec<String>,
) -> PyResult<Vec<Vec<PyHit>>> {
    let idx = Arc::clone(&index.inner);
    // Copy to owned bytes before releasing the GIL.
    let names: Vec<Vec<u8>> = names.into_iter().map(|s| s.into_bytes()).collect();
    let seqs: Vec<Vec<u8>> = seqs.into_iter().map(|s| s.into_bytes()).collect();
    let result = py.allow_threads(|| {
        let aligner = Aligner::new(&idx, &opts.inner);
        let mut buf = ThreadBuf::new();
        aligner.map_many(&mut buf, &names, &seqs)
    });
    let result = result.map_err(map_err)?;
    Ok(result
        .into_iter()
        .map(|hits| hits.into_iter().map(to_pyhit).collect())
        .collect())
}

#[pymodule]
fn _minibwa(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyIndex>()?;
    m.add_class::<PyOpts>()?;
    m.add_class::<PyHit>()?;
    m.add_function(wrap_pyfunction!(map, m)?)?;
    m.add_function(wrap_pyfunction!(map_pair, m)?)?;
    m.add_function(wrap_pyfunction!(map_many, m)?)?;
    Ok(())
}
