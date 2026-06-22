# minibwa — agent notes

- Three crates: `minibwa-sys` (FFI + vendored C + shim), `minibwa` (safe), `minibwa-py` (pyo3, standalone).
- minibwa is pure C99: `cc::Build` must not enable C++; link `z`/`pthread`/`m` only.
- `mb_map` results are libc-allocated — free with `libc::free`; read `n_cigar` (not `cap`) cigar words.
- `ThreadBuf` uses `mb_tbuf_init(0)` (kalloc pool) so it works for batch too.
- Index build wraps `main_index` via the shim; it can `abort()` on bad input, so `Index::build_from_fasta` pre-validates paths.
- Default build is permissive (libsais only); GPL `bwtgen` is behind the `gpl-bwtgen` feature.
- Commands: `cargo ci-fmt`/`ci-lint`/`ci-test` (aliases in `.cargo/config.toml`) for the workspace; Python in `minibwa-py/` via pixi (`pixi run test` rebuilds then pytests; also `lint`/`typecheck`/`format-check`).
- `minibwa-py` is excluded from the workspace, so workspace `ci-fmt`/`ci-lint` skip its Rust — it's gated separately in `python.yml`. When editing `minibwa-py/src/lib.rs`, run `cargo ci-fmt`/`ci-lint` from `minibwa-py/` (rust-toolchain.toml pins 1.95.0, so local rustfmt matches CI).
- pyo3 ↔ stub: `_minibwa.pyi` is hand-written; pyo3 arg names are the Python kwargs, so keep names/types synced with the stub (read-only `#[pyo3(get)]` fields → `@property`). Take paths as `PathBuf` (accepts `str`/`os.PathLike`), not `&str`.
- Methylation: `map(meth=…)` and the index must agree — a conversion against a non-meth index runs but is silently wrong (single-read `mt` is unguarded).
