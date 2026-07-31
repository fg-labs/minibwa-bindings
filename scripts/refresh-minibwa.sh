#!/usr/bin/env bash
# Refresh the vendored minibwa C snapshot.
#
# Usage: scripts/refresh-minibwa.sh <commit> [--patch <sha>]... [--src <path>]
#
#   <commit>        upstream base commit (normally a release tag's commit)
#   --patch <sha>   fork commit to cherry-pick onto the base, repeatable, in
#                   apply order; recorded in minibwa-sys/vendor/PATCHES
#   --src <path>    local minibwa clone to vendor from instead of cloning
#
# The pin is recorded as base (COMMIT) plus patches (PATCHES) rather than as a
# single already-rebased SHA, so it stays reproducible: a rebased fork commit
# exists only wherever it was created, while base + patch list replays anywhere.
set -euo pipefail

usage() {
  echo "usage: refresh-minibwa.sh <commit> [--patch <sha>]... [--src <path>]" >&2
  exit 2
}

[[ $# -ge 1 ]] || usage
COMMIT="$1"
shift
[[ "$COMMIT" == -* ]] && usage

SRC=""
PATCHES=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --patch)
      [[ $# -ge 2 ]] || usage
      PATCHES+=("$2")
      shift 2
      ;;
    --src)
      [[ $# -ge 2 ]] || usage
      SRC="$2"
      shift 2
      ;;
    *)
      echo "ERROR: unexpected argument: $1" >&2
      usage
      ;;
  esac
done

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DEST="$ROOT/minibwa-sys/vendor/minibwa"
META="$(dirname "$DEST")"

tmp="$(mktemp -d)"
CHECKOUT="$tmp/src"
trap 'rm -rf "$tmp"' EXIT

if [[ -n "$SRC" ]]; then
  # Clone rather than `git worktree add`: cherry-picking below writes commits, and
  # a clone keeps that entirely out of the caller's repository.
  git clone -q --shared --no-checkout "$SRC" "$CHECKOUT"
else
  # lh3 is the canonical upstream; the nh13 fork carries patches not yet merged
  # upstream (e.g. the AVX2/AVX-512 ksw_extd2 runtime dispatch, lh3 PR #20), so
  # fetch both remotes to resolve <commit> and any --patch regardless of which
  # one they live on.
  git clone -q --no-checkout https://github.com/lh3/minibwa "$CHECKOUT"
  git -C "$CHECKOUT" remote add fork https://github.com/nh13/minibwa
  git -C "$CHECKOUT" fetch -q --tags fork
fi

git -C "$CHECKOUT" -c advice.detachedHead=false checkout -q --detach "$COMMIT"
BASE="$(git -C "$CHECKOUT" rev-parse HEAD)"

# Replay the carried fork patches one at a time so a failure can name the
# offending commit. A committer identity is needed only because cherry-pick
# writes commits; nothing here is ever pushed.
for sha in ${PATCHES[@]+"${PATCHES[@]}"}; do
  if git -C "$CHECKOUT" \
      -c user.name=refresh-minibwa -c user.email=refresh-minibwa@invalid \
      cherry-pick "$sha"; then
    continue
  fi
  # cherry-pick also stops non-zero when the patch is already contained in the
  # base — the expected outcome once upstream merges a carried patch. That
  # leaves nothing staged, which distinguishes it from a real conflict.
  if git -C "$CHECKOUT" diff --cached --quiet && git -C "$CHECKOUT" diff --quiet; then
    git -C "$CHECKOUT" cherry-pick --abort 2>/dev/null || true
    echo "ERROR: patch $sha is already contained in $COMMIT; drop it from --patch." >&2
  else
    git -C "$CHECKOUT" cherry-pick --abort 2>/dev/null || true
    echo "ERROR: cherry-picking $sha onto $COMMIT conflicts; rebase it first." >&2
  fi
  exit 1
done

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

echo "$BASE" > "$META/COMMIT"

# Resolve patches to full SHAs so PATCHES is unambiguous years later.
resolved=()
for sha in ${PATCHES[@]+"${PATCHES[@]}"}; do
  resolved+=("$(git -C "$CHECKOUT" rev-parse "$sha")")
done

{
  echo "# Fork commits cherry-picked onto the upstream base recorded in COMMIT,"
  echo "# in apply order. Reproduce the vendored tree with:"
  echo "#"
  printf '#   scripts/refresh-minibwa.sh %s' "$BASE"
  for sha in ${resolved[@]+"${resolved[@]}"}; do
    printf ' --patch %s' "$sha"
  done
  printf '\n#\n'
  for sha in ${resolved[@]+"${resolved[@]}"}; do
    printf '%s  %s\n' "$sha" "$(git -C "$CHECKOUT" log -1 --format=%s "$sha")"
  done
} > "$META/PATCHES"

echo "Vendored minibwa @ $BASE (+${#PATCHES[@]} patch(es)) into $DEST"
