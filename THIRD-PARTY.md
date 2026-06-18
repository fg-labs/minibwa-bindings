# Third-party software

`minibwa-sys` vendors the C source of **minibwa** under
`minibwa-sys/vendor/minibwa/` (the pinned upstream commit is recorded in
`minibwa-sys/vendor/COMMIT`). That tree, and the suffix-array library it bundles,
carry their own licenses, reproduced below. Upstream's complete third-party
notices are in `minibwa-sys/vendor/minibwa/LICENSE.txt`.

## minibwa

- Project: <https://github.com/lh3/minibwa> (Heng Li)
- License: MIT — Copyright (c) Dana-Farber Cancer Institute
- Vendored in full at `minibwa-sys/vendor/minibwa/`.
- The pinned commit tracks lh3 `master` plus the AVX2/AVX-512 `ksw_extd2`
  runtime-dispatch patch (lh3 PR #20), carried on the `nh13/minibwa` fork until
  it merges upstream.

## libsais (compiled)

- License: Apache-2.0 — Copyright (c) Ilya Grebnov
- The suffix-array constructor used for index building (`libsais.c` /
  `libsais64.c`). This is the default, permissive index-build path.

## GPL components (vendored, NOT compiled by default)

- `bwtgen.c` — GPL-2 (BWT-SW, Wong Chi Kwong)
- `QSufSort.c` — HPND-style (N. Jesper Larsson)

These files are present in the vendored tree but are **not** compiled in the
default build, which uses the Apache-2.0 `libsais` path. They are compiled only
when the `minibwa-sys` `gpl-bwtgen` Cargo feature is enabled. **Enabling
`gpl-bwtgen` changes the effective license of the resulting binary to GPL-2.**
Do not enable it for artifacts distributed under this repository's MIT license.

This repository itself is licensed under the MIT License (see `LICENSE`).
