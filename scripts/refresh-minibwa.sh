#!/usr/bin/env bash
# Refresh the vendored minibwa C snapshot.
# Usage: scripts/refresh-minibwa.sh <commit> [local-source-path]
set -euo pipefail

COMMIT="${1:?usage: refresh-minibwa.sh <commit> [local-source-path]}"
SRC="${2:-}"
DEST="$(cd "$(dirname "$0")/.." && pwd)/minibwa-sys/vendor/minibwa"

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

if [[ -n "$SRC" ]]; then
  git -C "$SRC" worktree add --detach "$tmp/src" "$COMMIT"
  trap 'git -C "$SRC" worktree remove --force "$tmp/src" 2>/dev/null || true; rm -rf "$tmp"' EXIT
  CHECKOUT="$tmp/src"
else
  # lh3 is the canonical upstream; the nh13 fork carries patches not yet merged
  # upstream (e.g. the AVX2/AVX-512 ksw_extd2 runtime dispatch, lh3 PR #20), so
  # fetch both remotes to resolve <commit> regardless of which one it lives on.
  git clone -q https://github.com/lh3/minibwa "$tmp/src"
  git -C "$tmp/src" remote add fork https://github.com/nh13/minibwa
  git -C "$tmp/src" fetch -q fork
  git -C "$tmp/src" checkout --detach "$COMMIT"
  CHECKOUT="$tmp/src"
fi

# Reject a dirty tree (non-reproducible vendor).
if [[ -n "$(git -C "$CHECKOUT" status --porcelain)" ]]; then
  echo "ERROR: source tree at $COMMIT is dirty; refusing to vendor." >&2
  exit 1
fi

rm -rf "$DEST"
mkdir -p "$DEST"
# Copy all C sources/headers at top level; drop subtrees and artifacts we do not compile.
rsync -a \
  --exclude 'mimalloc/' --exclude 'api-test/' --exclude 'test/' --exclude 'tex/' \
  --exclude '.git' --exclude '.github/' --exclude '.gitignore' \
  --exclude 'dev.md' --exclude 'minibwa.1' \
  --exclude '*.o' --exclude '*.a' --exclude '/minibwa' \
  "$CHECKOUT"/ "$DEST"/

echo "$COMMIT" > "$(dirname "$DEST")/COMMIT"
echo "Vendored minibwa @ $COMMIT into $DEST"
