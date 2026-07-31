# minibwa-upstream-watch

A local composite action that keeps the vendored minibwa C snapshot honest about
upstream releases. Driven nightly by `.github/workflows/upstream-watch.yml`.

## What it does

1. Reads the pin: `minibwa-sys/vendor/COMMIT` (upstream base) plus
   `minibwa-sys/vendor/PATCHES` (fork commits carried on top, optional).
2. Clones `lh3/minibwa` with `nh13/minibwa` as a second remote, since the pin may
   reference a fork commit upstream cannot resolve on its own.
3. Determines the vendored release line as **the newest release tag that is an
   ancestor of the pin**. This is pure git ancestry rather than a version string,
   because `MB_VERSION` describes the commit, not the release the snapshot came
   from — a pin taken between releases reports a version that was never released.
4. Compares that against upstream's newest published release. Drafts and
   prereleases are ignored.
5. If behind, replays the carried fork patches onto the new release tag:
   - **clean** — re-vendors and opens or refreshes a draft PR on
     `chore/update-vendored-minibwa`;
   - **conflict** — opens no PR, and the tracking issue says so prominently.
   Patches upstream has since merged are recognised by patch id (`git cherry`)
   and dropped rather than re-applied.
6. Files or refreshes a single tracking issue, deduplicated on the
   `upstream-minibwa` label plus a marker comment in the body. When upstream
   moves on, the issue's title and body are rewritten in place and a comment
   records the bump. When the pin catches up — whether or not it was this PR
   that landed it — both the issue and any open refresh PR are closed.

## Why base + patches

Recording a single already-rebased SHA would be unreproducible: that commit
exists only wherever the rebase happened. Base plus an ordered patch list can be
replayed by anyone, which is what `scripts/refresh-minibwa.sh` does.

## Safety

The refresh branch is force-pushed only when every commit it adds on top of
`main` is the bot's. A human follow-up commit — likely, since a C API bump often
needs a shim change on the same branch — makes the action leave the branch alone
and say so. The ownership check fails **closed**: if the branch's remote-tracking
ref cannot be read, the push is refused rather than assumed safe. When a refresh
is skipped for either reason, the PR body is left untouched, so it never claims a
re-vendor that did not happen.

Patch shas recorded in `PATCHES` are always the original fork commits. Replaying
a pin mints new shas local to the runner, and those must never be written to the
pin — they would resolve nowhere.

## Running it locally

```sh
python3 .github/actions/minibwa-upstream-watch/watch.py \
  --repo fg-labs/minibwa-bindings --repo-root . --dry-run
```

This renders the issue and PR it would file and writes nothing, using whatever
`gh` is authenticated as. The same is available in CI via the workflow's
`dry-run` input.

## Tests

```sh
cd .github/actions/minibwa-upstream-watch && python3 -m unittest discover
```

Stdlib only. Ancestry, cherry-pick, and patch-id behaviour are exercised against
throwaway git repositories built in a temp directory, and `refresh-minibwa.sh` is
tested end to end against a fake upstream — no fixtures are committed.
