"""Type stubs for the compiled minibwa._minibwa extension module."""

import os

class Index:
    """A minibwa reference index."""

    @staticmethod
    def build(
        fasta: str | os.PathLike[str],
        prefix: str | os.PathLike[str],
        meth: bool = False,
        threads: int = 1,
    ) -> None:
        """Build a minibwa index from a FASTA file."""
        ...

    @staticmethod
    def load(prefix: str | os.PathLike[str], meth: bool = False) -> "Index":
        """Load a previously built minibwa index from disk."""
        ...

    def n_contigs(self) -> int:
        """Return the number of contigs (sequences) in the index."""
        ...

    def contigs(self) -> list[tuple[str, int]]:
        """Return a list of (name, length) tuples for every contig."""
        ...

class Opts:
    """Alignment options for minibwa."""

    def __init__(self, preset: str | None = None) -> None:
        """Create alignment options, optionally from a named preset."""
        ...

    def set_paired(self, on: bool) -> None:
        """Enable or disable paired-end mode."""
        ...

    def set_methylation(self, on: bool) -> None:
        """Enable or disable bisulfite/EM-seq methylation mode."""
        ...

    def set_min_seed_len(self, v: int) -> None:
        """Set the minimum seed length for BWT seeding."""
        ...

    def set_out_n(self, v: int) -> None:
        """Set the maximum number of secondary alignments to output."""
        ...

    def set_match_score(self, score: int) -> None:
        """Set the Smith-Waterman match score."""
        ...

    def set_mismatch_penalty(self, penalty: int) -> None:
        """Set the Smith-Waterman mismatch penalty."""
        ...

    def set_gap_open(self, penalty: int) -> None:
        """Set the gap-open penalty."""
        ...

    def set_gap_extend(self, extend: int) -> None:
        """Set the gap-extend penalty."""
        ...

    def set_pe_insert_size(self, avg: int, std: int, lo: int, hi: int) -> None:
        """Set paired-end insert-size parameters for proper-pair classification."""
        ...

class Hit:
    """One alignment result."""

    # Read-only getters (pyo3 `#[pyo3(get)]`), declared as properties so a type
    # checker rejects assignment, matching the runtime.
    @property
    def contig(self) -> str | None: ...
    @property
    def ref_start(self) -> int: ...
    @property
    def ref_end(self) -> int: ...
    @property
    def query_start(self) -> int: ...
    @property
    def query_end(self) -> int: ...
    @property
    def reverse(self) -> bool: ...
    @property
    def mapq(self) -> int: ...
    @property
    def score(self) -> int: ...
    @property
    def n_sub(self) -> int: ...
    @property
    def proper_pair(self) -> bool: ...
    @property
    def is_primary(self) -> bool: ...
    @property
    def is_secondary(self) -> bool: ...
    @property
    def is_supplementary(self) -> bool: ...
    @property
    def cigar(self) -> list[tuple[str, int]]: ...
    @property
    def strand(self) -> str:
        """Return '+' for forward-strand or '-' for reverse."""
        ...

    @property
    def cigar_string(self) -> str:
        """Return the CIGAR as a SAM-style string (e.g. '150M'), or '*' if empty."""
        ...

    def __repr__(self) -> str: ...

def map(  # noqa: A001
    index: Index,
    opts: Opts,
    name: str,
    seq: str,
    meth: str = "none",
) -> list[Hit]:
    """Align one read against the index."""
    ...

def map_pair(
    index: Index,
    opts: Opts,
    name1: str,
    seq1: str,
    name2: str,
    seq2: str,
) -> tuple[list[Hit], list[Hit]]:
    """Align a read pair (paired-end) against the index."""
    ...

def map_many(
    index: Index,
    opts: Opts,
    names: list[str],
    seqs: list[str],
) -> list[list[Hit]]:
    """Align many independent reads in a single batched call."""
    ...
