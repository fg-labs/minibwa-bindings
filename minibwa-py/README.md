# minibwa

Python bindings for [minibwa](https://github.com/lh3/minibwa), a minimal BWA-style short- and long-read aligner. Built on the Rust `minibwa` crate via PyO3; the alignment runs in native code with the GIL released.

## Install

```bash
pip install minibwa
```

## Quick start

```python
import minibwa

# Build an index from a FASTA (writes <prefix>.l2b and <prefix>.mbw), then load it.
minibwa.Index.build("ref.fa", "ref")
idx = minibwa.Index.load("ref")

opts = minibwa.Opts()
for hit in minibwa.map(idx, opts, "read1", "ACGTACGT..."):
    print(hit.contig, hit.ref_start, hit.ref_end, hit.strand, hit.cigar_string)
```

- Paired-end: `minibwa.map_pair(idx, opts, n1, s1, n2, s2)`
- Batched (higher throughput): `minibwa.map_many(idx, opts, names, seqs)`
- Methylation (bisulfite): pass `meth=minibwa.Meth.C2T` / `Meth.G2A` (build/load the index with `meth=True`)

See the [project repository](https://github.com/fg-labs/minibwa-bindings) for the full API and documentation.

## License

MIT. The package vendors and links minibwa (MIT) and libsais (Apache-2.0); see the repository's `THIRD-PARTY.md` for full notices.
