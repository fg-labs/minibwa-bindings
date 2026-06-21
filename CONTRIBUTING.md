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

The release workflow uses [release-plz](https://release-plz.ctl.dev/) to manage the release PR and [git-cliff](https://git-cliff.org/) to generate the changelog. Rust crates publish to crates.io via OIDC Trusted Publishing (no stored token); Python wheels publish to PyPI via the `pypi.yml` workflow, also via Trusted Publishing.

The workspace crates (`minibwa-sys`, `minibwa`) are publishable. `minibwa-py` keeps `publish = false`: it is a standalone crate outside the workspace and ships only as a PyPI wheel, never to crates.io.

### One-time setup (before the first release)

crates.io Trusted Publishing has no "pending publisher" equivalent to PyPI's — a trusted publisher can only be configured for a crate that already exists. So each new crate must be bootstrapped with one manual publish before CI can take over:

1. **Bootstrap crates.io** (manual, once per crate). From a clean checkout of the version to release, with a crates.io API token in the environment, publish in dependency order:

   ```bash
   cargo publish -p minibwa-sys   # no workspace deps; publish first
   cargo publish -p minibwa       # depends on minibwa-sys
   ```

   This creates the crates and records you as owner.

2. **Add crates.io Trusted Publishing** for **both** `minibwa-sys` and `minibwa` (crate Settings → Trusted Publishing on crates.io): repository `fg-labs/minibwa-bindings`, workflow `publish.yml`, environment left blank (the `publish` job sets none). After this, CI publishes subsequent versions with no token.

3. **PyPI Trusted Publishing**: configure a publisher (a "pending publisher" works before the project exists) for project `minibwa` pointing at repo `fg-labs/minibwa-bindings`, workflow `pypi.yml`, environment `pypi`.

4. **GitHub `pypi` environment**: create an environment named `pypi` in the repo settings; the `pypi.yml` upload job runs in it.

5. **First wheels**: because the bootstrap in step 1 already published the version to crates.io, the next push to `main` makes `publish.yml` skip the (already-published) crates and therefore *not* create the GitHub release — so `pypi.yml` will not fire for that version. Cut the release once by hand to publish the first wheels:

   ```bash
   gh release create v<VERSION> --generate-notes --latest
   ```

   From the next version onward, the automated flow below handles everything.

### Per-release (automated)

1. **Version bump**: release-plz opens a release PR that bumps the version in `minibwa-sys/Cargo.toml` and `minibwa/Cargo.toml` together (they are always versioned in lockstep). The `publish.yml` job verifies this lockstep before publishing.

2. **minibwa-py version**: release-plz does NOT touch `minibwa-py` (it is excluded from the workspace). Bump the version manually in both `minibwa-py/Cargo.toml` and `minibwa-py/pyproject.toml` to match the workspace version before merging the release PR.

3. **Publish order**: merging the release PR triggers `publish.yml`, which publishes `minibwa-sys` first (no workspace deps), then `minibwa` (depends on `minibwa-sys`), skipping any crate already at the target version on crates.io, then creates the GitHub release.

4. **Python wheels**: the GitHub release's `release: published` event triggers `pypi.yml`, which builds manylinux_2_28 wheels for x86_64 and aarch64 Linux plus an arm64 macOS wheel, then uploads via Trusted Publishing to PyPI.

5. **GPL feature must not be enabled in published artifacts**: the `gpl-bwtgen` feature changes the effective license to GPL. Never enable it in published crates or wheels. The feature flag exists solely for local opt-in use.
