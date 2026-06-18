# minibwa — agent notes

- Three crates: `minibwa-sys` (FFI + vendored C + shim), `minibwa` (safe), `minibwa-py` (pyo3, standalone).
- minibwa is pure C99: `cc::Build` must not enable C++; link `z`/`pthread`/`m` only.
- `mb_map` results are libc-allocated — free with `libc::free`; read `n_cigar` (not `cap`) cigar words.
- `ThreadBuf` uses `mb_tbuf_init(0)` (kalloc pool) so it works for batch too.
- Index build wraps `main_index` via the shim; it can `abort()` on bad input, so `Index::build_from_fasta` pre-validates paths.
- Default build is permissive (libsais only); GPL `bwtgen` is behind the `gpl-bwtgen` feature.
