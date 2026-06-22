"""Tests for minibwa alignment functions."""

import os
import pathlib
import tempfile

import minibwa
import pytest


def _synthetic_reference(n: int, seed: int = 7) -> str:
    x = seed | 1
    out = []
    for _ in range(n):
        x = (x * 1664525 + 1013904223) & 0xFFFFFFFF
        out.append("ACGT"[(x >> 16) & 3])
    return "".join(out)


def _revcomp(s: str) -> str:
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return s.translate(comp)[::-1]


def _build_index(tmpdir: str, ref: str) -> tuple[minibwa.Index, minibwa.Opts]:
    """Build and load a test index; return (idx, opts)."""
    fa = os.path.join(tmpdir, "ref.fa")
    with open(fa, "w") as fh:
        fh.write(">chr1\n")
        fh.write(ref + "\n")
    prefix = os.path.join(tmpdir, "idx")
    minibwa.Index.build(fa, prefix, meth=False, threads=1)
    idx = minibwa.Index.load(prefix, meth=False)
    opts = minibwa.Opts()
    return idx, opts


def test_build_load_map() -> None:
    ref = _synthetic_reference(5000)
    with tempfile.TemporaryDirectory() as d:
        fa = os.path.join(d, "ref.fa")
        with open(fa, "w") as fh:
            fh.write(">chr1\n")
            fh.write(ref + "\n")
        prefix = os.path.join(d, "idx")
        minibwa.Index.build(fa, prefix, meth=False, threads=1)
        idx = minibwa.Index.load(prefix, meth=False)
        assert idx.n_contigs() == 1
        assert idx.contigs()[0][0] == "chr1"

        opts = minibwa.Opts()
        query = ref[1000:1150]
        hits = minibwa.map(idx, opts, "q1", query)
        assert len(hits) >= 1
        h = hits[0]
        assert h.contig == "chr1"
        assert 900 <= h.ref_start <= 1100
        assert len(h.cigar) >= 1
        assert h.cigar[0][0] in "M=X"
        assert repr(hits[0]).startswith("Hit(")


def test_build_load_accepts_pathlike() -> None:
    # Bioinformatics callers pass pathlib.Path; build/load must accept os.PathLike,
    # not just str.
    ref = _synthetic_reference(5000)
    with tempfile.TemporaryDirectory() as d:
        fa = pathlib.Path(d) / "ref.fa"
        fa.write_text(">chr1\n" + ref + "\n")
        prefix = pathlib.Path(d) / "idx"
        minibwa.Index.build(fa, prefix, meth=False, threads=1)
        idx = minibwa.Index.load(prefix, meth=False)
        assert idx.n_contigs() == 1


def test_map_pair() -> None:
    ref = _synthetic_reference(5000)
    with tempfile.TemporaryDirectory() as d:
        fa = os.path.join(d, "ref.fa")
        with open(fa, "w") as fh:
            fh.write(">chr1\n")
            fh.write(ref + "\n")
        prefix = os.path.join(d, "idx")
        minibwa.Index.build(fa, prefix, meth=False, threads=1)
        idx = minibwa.Index.load(prefix, meth=False)

        opts = minibwa.Opts()
        r1 = ref[1000:1150]
        r2 = _revcomp(ref[1400:1550])
        result = minibwa.map_pair(idx, opts, "p1", r1, "p2", r2)
        assert isinstance(result, tuple) and len(result) == 2
        hits1, hits2 = result
        assert len(hits1) >= 1
        assert len(hits2) >= 1
        assert hits1[0].contig == "chr1"
        assert 900 <= hits1[0].ref_start <= 1100
        assert hits2[0].contig == "chr1"
        assert 1300 <= hits2[0].ref_start <= 1500


def test_map_many() -> None:
    ref = _synthetic_reference(5000)
    with tempfile.TemporaryDirectory() as d:
        idx, opts = _build_index(d, ref)

        # Three distinct 150 bp queries from well-separated regions.
        q0 = ref[500:650]
        q1 = ref[2000:2150]
        q2 = ref[3500:3650]

        result = minibwa.map_many(idx, opts, ["r0", "r1", "r2"], [q0, q1, q2])

        assert len(result) == 3
        for i, hits in enumerate(result):
            assert len(hits) >= 1, f"read r{i} produced no hits"
        expected_offsets = [500, 2000, 3500]
        for i, (hits, exp) in enumerate(zip(result, expected_offsets)):
            pos = hits[0].ref_start
            assert abs(pos - exp) <= 100, f"read r{i}: expected pos ~{exp}, got {pos}"
            assert hits[0].contig == "chr1"


def test_meth_enum() -> None:
    ref = _synthetic_reference(5000)
    with tempfile.TemporaryDirectory() as d:
        idx, opts = _build_index(d, ref)
        query = ref[1000:1150]

        # Meth.C2T is a str subclass — passes through to the Rust &str matcher.
        assert minibwa.Meth.C2T == "c2t"  # type: ignore[comparison-overlap]
        assert minibwa.Meth.G2A == "g2a"  # type: ignore[comparison-overlap]
        assert minibwa.Meth.NONE == "none"  # type: ignore[comparison-overlap]

        # Using the enum value as the meth argument must not raise.
        hits = minibwa.map(idx, opts, "r", query, meth=minibwa.Meth.C2T)
        assert isinstance(hits, list)


def test_hit_strand_and_fields() -> None:
    ref = _synthetic_reference(5000)
    with tempfile.TemporaryDirectory() as d:
        idx, opts = _build_index(d, ref)
        query = ref[1000:1150]

        hits = minibwa.map(idx, opts, "r", query)
        assert len(hits) >= 1
        h = hits[0]

        # Strand accessors.
        assert h.strand == "+"
        assert h.reverse is False

        # Coordinate fields.
        assert h.ref_end > h.ref_start
        assert h.query_start == 0
        assert h.query_end == len(query)

        # Quality / score.
        assert h.mapq > 0
        assert h.score > 0

        # Alignment flags.
        assert h.is_primary
        assert not h.is_secondary
        assert not h.is_supplementary

        # cigar_string mirrors the cigar tuples as a SAM-style string. The
        # empty-CIGAR "*" branch isn't reachable from Python (unmapped reads
        # come back as an empty list, never a Hit); the Rust crate covers it
        # directly (hit.rs::cigar_string_empty_is_star).
        assert h.cigar_string == "".join(f"{n}{op}" for op, n in h.cigar)
        assert h.cigar_string != "*"


def test_error_paths() -> None:
    ref = _synthetic_reference(5000)
    with tempfile.TemporaryDirectory() as d:
        idx, opts = _build_index(d, ref)
        query = ref[1000:1150]

        # Empty sequence → ValueError.
        with pytest.raises(ValueError):
            minibwa.map(idx, opts, "r", "")

        # Unknown meth string → ValueError.
        with pytest.raises(ValueError):
            minibwa.map(idx, opts, "r", query, meth="bogus")

        # Non-existent index prefix → RuntimeError.
        with pytest.raises(RuntimeError):
            minibwa.Index.load("/no/such/path")

        # Unknown preset → ValueError.
        with pytest.raises(ValueError):
            minibwa.Opts(preset="nope")


def test_error_paths_names_and_batches() -> None:
    ref = _synthetic_reference(5000)
    with tempfile.TemporaryDirectory() as d:
        idx, opts = _build_index(d, ref)
        query = ref[1000:1150]

        # NUL byte in a read name → ValueError.
        with pytest.raises(ValueError):
            minibwa.map(idx, opts, "bad\x00name", query)

        # map_pair with an empty mate → ValueError.
        with pytest.raises(ValueError):
            minibwa.map_pair(idx, opts, "p1", query, "p2", "")

        # map_many with mismatched names/seqs lengths → ValueError.
        with pytest.raises(ValueError):
            minibwa.map_many(idx, opts, ["r0", "r1"], [query])

        # map_many with an empty sequence → ValueError.
        with pytest.raises(ValueError):
            minibwa.map_many(idx, opts, ["r0"], [""])


def test_opts_setters() -> None:
    ref = _synthetic_reference(5000)
    with tempfile.TemporaryDirectory() as d:
        idx, _opts = _build_index(d, ref)

        # Chain-free style: mutate in place, no exception expected.
        o = minibwa.Opts()
        o.set_min_seed_len(19)
        o.set_out_n(5)
        o.set_match_score(1)
        o.set_mismatch_penalty(4)
        o.set_gap_open(6)
        o.set_gap_extend(1)
        o.set_pe_insert_size(400, 50, 200, 600)

        # Opts must still be usable for alignment after mutation.
        query = ref[1000:1150]
        hits = minibwa.map(idx, o, "r", query)
        assert isinstance(hits, list)


def test_opts_setters_keyword_form() -> None:
    # The documented keyword names (mirrored in the .pyi stub) must match the
    # actual parameter names, so callers and type checkers agree.
    o = minibwa.Opts()
    o.set_min_seed_len(v=19)
    o.set_out_n(v=5)
    o.set_match_score(score=1)
    o.set_mismatch_penalty(penalty=4)
    o.set_gap_open(penalty=6)
    o.set_gap_extend(extend=1)
    o.set_pe_insert_size(avg=400, std=50, lo=200, hi=600)
