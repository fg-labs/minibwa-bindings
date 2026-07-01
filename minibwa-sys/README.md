# minibwa-sys

[![crates.io](https://img.shields.io/crates/v/minibwa-sys.svg)](https://crates.io/crates/minibwa-sys)
[![docs.rs](https://docs.rs/minibwa-sys/badge.svg)](https://docs.rs/minibwa-sys)

Low-level, unsafe FFI bindings to [minibwa](https://github.com/lh3/minibwa), Heng Li's lightweight bwa.

This is the `-sys` crate: it vendors the minibwa C source, compiles it via `build.rs`, and exposes the raw `bindgen`-generated bindings. Most users should depend on the safe [`minibwa`](https://crates.io/crates/minibwa) wrapper instead, which builds structured alignment `Hit`s on top of these bindings.

## What it does

- Vendors the minibwa C source under `vendor/` and a small `shim/` layer.
- Compiles the C at build time (`build.rs`), then runs `bindgen` to generate the Rust FFI declarations.
- On docs.rs, where system zlib is unavailable, it runs `bindgen` only and skips C compilation.

## Licensing

MIT. The default build uses only the Apache-2.0 `libsais` for indexing. The GPL `bwtgen` indexing path is available behind the opt-in `gpl-bwtgen` feature, which changes the effective license to GPL.

## Updating the vendored source

From the workspace root:

```sh
scripts/refresh-minibwa.sh <clean-commit> [local-source-path]
```
