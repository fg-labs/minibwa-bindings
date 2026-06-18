# minibwa

[![CI](https://github.com/fg-labs/minibwa-bindings/actions/workflows/check.yml/badge.svg)](https://github.com/fg-labs/minibwa-bindings/actions/workflows/check.yml)
[![crates.io](https://img.shields.io/crates/v/minibwa.svg)](https://crates.io/crates/minibwa)
[![docs.rs](https://docs.rs/minibwa/badge.svg)](https://docs.rs/minibwa)
[![PyPI](https://img.shields.io/pypi/v/minibwa.svg)](https://pypi.org/project/minibwa)

Rust and Python bindings for [minibwa](https://github.com/lh3/minibwa).

- `minibwa-sys` — low-level FFI over vendored minibwa C (compiled via `build.rs`).
- `minibwa` — safe Rust wrapper exposing structured alignment `Hit`s.
- `minibwa-py` — Python bindings (pyo3 + maturin).

## Quick start (Rust)

```rust
use minibwa::{Aligner, Index, Meth, Opts, ThreadBuf};

Index::build_from_fasta("ref.fa", "ref", false, 4)?;
let idx = Index::load("ref", false)?;
let opts = Opts::new();
let aligner = Aligner::new(&idx, &opts);
let mut buf = ThreadBuf::new();
for hit in aligner.map(&mut buf, b"read1", b"ACGT...", Meth::None)? {
    println!("{:?} {}..{}", hit.contig, hit.ref_start, hit.ref_end);
}
# Ok::<(), minibwa::Error>(())
```

## Quick start (Python)

```python
import minibwa

# Build a BWA index from a FASTA file (one-time)
minibwa.Index.build("ref.fa", "ref")

# Load the index
idx = minibwa.Index.load("ref")

# Configure alignment options
opts = minibwa.Opts()

# Map a single read
hits = minibwa.map(idx, opts, "read1", "ACGT...")
for h in hits:
    print(h)
```

## Licensing

MIT. The default build uses only Apache-2.0 `libsais` for indexing. The GPL
`bwtgen` path is available behind the opt-in `gpl-bwtgen` feature, which changes
the effective license to GPL.

## Updating the vendored source

`scripts/refresh-minibwa.sh <clean-commit> [local-source-path]`
