# Contributing

## Prerequisites
- Rust stable (see `rust-toolchain.toml`); a C99 compiler; zlib headers.
- For Python: `pixi` (or a venv with `maturin` + `pytest`).

## Workflow
- `cargo ci-fmt && cargo ci-lint && cargo ci-test` before pushing.
- Install the pre-commit hook: `scripts/install-hooks.sh`.
- Conventional Commits; feature branches only; never commit to `main`/`dev`.

## Updating vendored minibwa
- `scripts/refresh-minibwa.sh <clean-commit> [local-source-path]`, then run the
  full test suite — the bindgen layout tests catch struct drift.

## Features
- `minibwa-sys/gpl-bwtgen`: vendors GPL `bwtgen`/`QSufSort`. **Changes the
  effective license to GPL.** Off by default.
- `minibwa-sys/openmp`: threaded `libsais` indexing (`-fopenmp` + `-DLIBSAIS_OPENMP`).

## Release process

The release workflow uses [release-plz](https://release-plz.ctl.dev/) to manage the release PR and [git-cliff](https://git-cliff.org/) to generate the changelog. Python wheels are published via the `pypi.yml` workflow.

1. **Flip `publish = false`**: Before releasing, remove `publish = false` from `minibwa-sys/Cargo.toml` and `minibwa/Cargo.toml`. Do NOT remove it from `minibwa-py/Cargo.toml` (it is a standalone crate outside the workspace, published separately as a wheel).

2. **Version bump**: release-plz opens a release PR that bumps the version in `minibwa-sys/Cargo.toml` and `minibwa/Cargo.toml` together (they are always versioned in lockstep). The `publish.yml` job verifies this lockstep before publishing.

3. **minibwa-py version**: release-plz does NOT touch `minibwa-py` (it is excluded from the workspace). Bump the version manually in both `minibwa-py/Cargo.toml` and `minibwa-py/pyproject.toml` to match the workspace version before merging the release PR.

4. **Publish order**: `minibwa-sys` is published first (no workspace deps), then `minibwa` (depends on `minibwa-sys`). The `publish.yml` workflow enforces this order and skips any crate already at the target version on crates.io.

5. **Python wheels**: Once the GitHub release is created by `publish.yml`, the `release: published` event triggers `pypi.yml`, which builds manylinux_2_28 wheels for x86_64 and aarch64 Linux plus an arm64 macOS wheel, then uploads via OIDC Trusted Publishing to PyPI.

6. **GPL feature must not be enabled in published artifacts**: The `gpl-bwtgen` feature changes the effective license to GPL. Never enable it in published crates or wheels. The feature flag exists solely for local opt-in use.

7. **PyPI Trusted Publishing setup** (one-time): The PyPI project `minibwa` must have a Trusted Publisher configured pointing at the `fg-labs/minibwa-bindings` repo, the `pypi.yml` workflow file, and the `pypi` GitHub environment before the first wheel upload.
