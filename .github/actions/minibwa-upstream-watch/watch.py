#!/usr/bin/env python3
"""Watch upstream minibwa for new releases and prepare the vendor refresh.

The vendored C snapshot is pinned by ``minibwa-sys/vendor/COMMIT`` (an upstream
base commit) plus an optional ``minibwa-sys/vendor/PATCHES`` listing fork commits
cherry-picked on top. This script decides whether that pin is behind upstream's
newest release and, when it is:

* cherry-picks the carried fork patches onto the new release tag and, if that is
  clean, opens or refreshes a draft PR with the re-vendored tree;
* opens or refreshes a single tracking issue describing the gap.

A cherry-pick conflict is not an error: it suppresses the PR and makes the issue
say loudly that the refresh needs a human.

The module is split into a pure core (parsing, selection, rendering, planning)
and thin ``Git``/``GitHub`` subprocess seams so the core can be unit-tested
without a network.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Callable, Iterable, Sequence

# Identifies the tracking issue this script owns. Kept in the issue body so the
# issue is still recognisable if someone edits the title or strips the label.
ISSUE_MARKER = "<!-- minibwa-upstream-watch -->"

# A single long-lived branch, so the refresh PR updates in place rather than
# accumulating one stale PR per upstream release.
BRANCH = "chore/update-vendored-minibwa"

COMMIT_FILE = Path("minibwa-sys/vendor/COMMIT")
PATCHES_FILE = Path("minibwa-sys/vendor/PATCHES")
REFRESH_SCRIPT = Path("scripts/refresh-minibwa.sh")

# Upstream names releases "Minibwa-0.6 (r416)"; we prefer that to the bare tag.
RELEASE_NAME_RE = re.compile(r"minibwa[-\s]*([0-9][0-9.]*)\s*\(\s*(r\d+)\s*\)", re.IGNORECASE)

# Release notes can be long; issues stay readable if each is capped.
MAX_NOTE_CHARS = 3000


class WatchError(RuntimeError):
    """A failure that should fail the workflow, as opposed to a normal outcome."""


# --------------------------------------------------------------------------- #
# Data
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class Release:
    """A published upstream GitHub release."""

    tag: str
    name: str
    body: str
    published_at: str
    url: str
    commit: str = ""

    @property
    def label(self) -> str:
        """Human label, e.g. ``0.6 (r416)``, falling back to the tag."""
        match = RELEASE_NAME_RE.search(self.name or "")
        if match:
            return f"{match.group(1)} ({match.group(2)})"
        return self.tag

    @property
    def version(self) -> str:
        """Version string for commit subjects, e.g. ``0.6-r416``."""
        match = RELEASE_NAME_RE.search(self.name or "")
        if match:
            return f"{match.group(1)}-{match.group(2)}"
        return self.tag.removeprefix("v")


@dataclass(frozen=True)
class Patch:
    """A fork commit carried on top of the vendored upstream base."""

    sha: str
    subject: str

    def line(self) -> str:
        """Render as a ``PATCHES`` file line."""
        return f"{self.sha}  {self.subject}"


@dataclass(frozen=True)
class Pin:
    """The current vendor pin: an upstream base plus carried fork patches."""

    base: str
    patches: tuple[Patch, ...] = ()


@dataclass(frozen=True)
class Plan:
    """What the nightly run intends to do."""

    latest: Release
    current: Release | None
    intervening: tuple[Release, ...] = ()
    carried: tuple[Patch, ...] = ()
    upstreamed: tuple[Patch, ...] = ()
    conflict: Patch | None = None
    candidate: str | None = None

    @property
    def behind(self) -> bool:
        return self.current is None or self.current.tag != self.latest.tag

    @property
    def blocked(self) -> bool:
        return self.conflict is not None


# --------------------------------------------------------------------------- #
# Pure core
# --------------------------------------------------------------------------- #


def parse_patches(text: str) -> tuple[Patch, ...]:
    """Parse a ``PATCHES`` file.

    Each meaningful line is ``<sha>  <subject>``. Blank lines and ``#`` comments
    are ignored, so the file can carry explanatory notes. A missing or empty file
    (i.e. empty ``text``) means no carried patches.
    """
    patches: list[Patch] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        sha, _, subject = line.partition(" ")
        if not re.fullmatch(r"[0-9a-f]{7,40}", sha):
            raise WatchError(f"malformed line in {PATCHES_FILE}: {raw!r}")
        patches.append(Patch(sha=sha, subject=subject.strip()))
    return tuple(patches)


def published_releases(payload: Iterable[dict]) -> list[Release]:
    """Select real releases from the GitHub API payload, newest first.

    Drafts and prereleases are excluded so an upstream release candidate cannot
    trigger a vendor bump.
    """
    releases = [
        Release(
            tag=item["tagName"],
            name=item.get("name") or "",
            body=item.get("body") or "",
            published_at=item.get("publishedAt") or "",
            url=item.get("url") or "",
        )
        for item in payload
        if not item.get("isDraft") and not item.get("isPrerelease")
    ]
    releases.sort(key=lambda r: r.published_at, reverse=True)
    return releases


def newest_ancestor_release(
    releases: Sequence[Release],
    head: str,
    is_ancestor: Callable[[str, str], bool],
) -> Release | None:
    """Return the newest release whose tag is an ancestor of ``head``.

    This is how the vendored "release line" is determined: purely by git
    ancestry, so it stays correct when the pin is a fork commit that no version
    string describes accurately.
    """
    for release in releases:
        if release.commit and is_ancestor(release.commit, head):
            return release
    return None


def releases_between(releases: Sequence[Release], current: Release | None) -> tuple[Release, ...]:
    """Releases newer than ``current``, newest first.

    With no known current release every release is "intervening", which is the
    right answer for a pin that predates all tags.
    """
    if current is None:
        return tuple(releases)
    out: list[Release] = []
    for release in releases:
        if release.tag == current.tag:
            break
        out.append(release)
    return tuple(out)


def truncate(text: str, limit: int = MAX_NOTE_CHARS) -> str:
    """Trim release notes to keep the rendered issue readable.

    Closes a fence left open by the cut: an unterminated ``` would otherwise
    swallow every later section of the rendered issue body.
    """
    text = (text or "").strip()
    if len(text) <= limit:
        return text
    cut = text[:limit].rstrip()
    if cut.count("```") % 2:
        cut += "\n```"
    return cut + "\n\n_(truncated)_"


def issue_title(plan: Plan) -> str:
    return f"Update vendored minibwa to {plan.latest.label}"


def pr_title(plan: Plan) -> str:
    return f"chore(deps): update vendored minibwa to {plan.latest.version}"


def commit_subject(plan: Plan) -> str:
    return pr_title(plan)


def reproduce_command(plan: Plan) -> str:
    """The exact local invocation that reproduces the bot's vendored tree."""
    patches = "".join(f" \\\n    --patch {p.sha}" for p in plan.carried)
    return f"{REFRESH_SCRIPT} {plan.latest.commit}{patches}"


def _release_notes_section(plan: Plan) -> list[str]:
    lines = ["## Upstream release notes", ""]
    for release in plan.intervening:
        lines.append(f"### [{release.label}]({release.url})")
        lines.append("")
        notes = truncate(release.body)
        lines.append(notes if notes else "_No release notes._")
        lines.append("")
    return lines


def _patches_section(plan: Plan) -> list[str]:
    lines = ["## Carried fork patches", ""]
    if not plan.carried and not plan.upstreamed:
        lines += ["None — the pin is plain upstream.", ""]
        return lines
    if plan.carried:
        lines.append("Re-applied on top of the new release tag:")
        lines.append("")
        for patch in plan.carried:
            lines.append(f"- `{patch.sha[:8]}` {patch.subject}")
        lines.append("")
    if plan.upstreamed:
        lines.append("Now merged upstream (identical patch id) and therefore dropped:")
        lines.append("")
        for patch in plan.upstreamed:
            lines.append(f"- `{patch.sha[:8]}` {patch.subject}")
        lines.append("")
    return lines


def render_issue_body(plan: Plan, pr_number: int | None, stale_pr: int | None = None) -> str:
    """Render the tracking issue body."""
    current = plan.current.label if plan.current else "unknown (pin predates all releases)"
    lines = [
        ISSUE_MARKER,
        "",
        f"Vendored minibwa is **{current}**; upstream's newest release is "
        f"**[{plan.latest.label}]({plan.latest.url})**.",
        "",
        "| | Release | Commit |",
        "| --- | --- | --- |",
        f"| Vendored | {current} | `{(plan.current.commit[:8] if plan.current else '—')}` |",
        f"| Latest | {plan.latest.label} | `{plan.latest.commit[:8]}` |",
        "",
    ]

    if plan.blocked:
        assert plan.conflict is not None
        lines += [
            "## :warning: Automated refresh is blocked",
            "",
            "Re-applying the carried fork patches onto "
            f"`{plan.latest.tag}` conflicts, so **no pull request was opened**. "
            "The conflicting patch is:",
            "",
            f"- `{plan.conflict.sha[:8]}` {plan.conflict.subject}",
            "",
            "Resolve it by hand — rebase the patch onto the new tag, then run the "
            "refresh with the rebased commit:",
            "",
            "```sh",
            f"{REFRESH_SCRIPT} {plan.latest.commit} --patch <rebased-sha>",
            "```",
            "",
        ]
        if stale_pr is not None:
            lines += [
                f"An earlier automated refresh PR (#{stale_pr}) is still open and now "
                "targets an older release; it will not be updated while this conflict "
                "stands.",
                "",
            ]
    elif pr_number is not None:
        lines += [
            "## Prepared refresh",
            "",
            f"The re-vendored tree is ready for review in #{pr_number}.",
            "",
        ]

    lines += _patches_section(plan)
    lines += _release_notes_section(plan)

    if not plan.blocked:
        lines += [
            "## Reproduce locally",
            "",
            "```sh",
            reproduce_command(plan),
            "```",
            "",
        ]

    lines += ["_Filed automatically by the nightly upstream-minibwa watch._"]
    return "\n".join(lines)


def render_pr_body(plan: Plan, issue_number: int | None) -> str:
    """Render the draft PR body."""
    current = plan.current.label if plan.current else "unknown"
    lines = [
        f"Re-vendors the minibwa C snapshot from **{current}** to "
        f"**[{plan.latest.label}]({plan.latest.url})**.",
        "",
    ]
    if issue_number is not None:
        lines += [f"Tracking issue: #{issue_number}", ""]
    lines += _patches_section(plan)
    lines += [
        "## Reproduce locally",
        "",
        "```sh",
        reproduce_command(plan),
        "```",
        "",
        "## Review checklist",
        "",
        "- [ ] `cargo ci-test` passes and the FFI shim still matches upstream's API",
        "- [ ] `minibwa-sys/vendor/PATCHES` lists exactly the patches you expect",
        "- [ ] Upstream release notes reviewed for behaviour changes affecting the bindings",
        "- [ ] `THIRD-PARTY.md` still describes the vendored tree accurately",
        "",
    ]
    lines += _release_notes_section(plan)
    lines += ["_Opened automatically by the nightly upstream-minibwa watch._"]
    return "\n".join(lines)


# --------------------------------------------------------------------------- #
# I/O seams
# --------------------------------------------------------------------------- #


def _run(cmd: Sequence[str], cwd: Path | None = None, env: dict | None = None,
         check: bool = True) -> subprocess.CompletedProcess:
    """Run a command, capturing output, raising ``WatchError`` on failure."""
    result = subprocess.run(
        list(cmd),
        cwd=str(cwd) if cwd else None,
        env=env,
        capture_output=True,
        text=True,
    )
    if check and result.returncode != 0:
        raise WatchError(
            f"command failed ({result.returncode}): {' '.join(cmd)}\n"
            f"stdout: {result.stdout}\nstderr: {result.stderr}"
        )
    return result


class Git:
    """Thin wrapper around git in one working tree."""

    def __init__(self, root: Path) -> None:
        self.root = root

    def run(self, *args: str, check: bool = True) -> subprocess.CompletedProcess:
        return _run(["git", *args], cwd=self.root, check=check)

    def out(self, *args: str) -> str:
        return self.run(*args).stdout.strip()

    def is_ancestor(self, maybe_ancestor: str, head: str) -> bool:
        return self.run(
            "merge-base", "--is-ancestor", maybe_ancestor, head, check=False
        ).returncode == 0

    def resolve_tag_commit(self, tag: str) -> str:
        """Commit SHA for a tag, dereferencing annotated tags."""
        result = self.run("rev-list", "-n", "1", tag, check=False)
        return result.stdout.strip() if result.returncode == 0 else ""

    def subject(self, sha: str) -> str:
        return self.out("log", "-1", "--format=%s", sha)

    def cherry(self, upstream: str, head: str) -> list[Patch]:
        """Commits in ``head`` not upstream, by patch id (``git cherry`` ``+``)."""
        patches: list[Patch] = []
        for line in self.out("cherry", upstream, head).splitlines():
            mark, _, sha = line.partition(" ")
            sha = sha.strip()
            if mark == "+" and sha:
                patches.append(Patch(sha=sha, subject=self.subject(sha)))
        return patches


class GitHub:
    """Thin wrapper around the ``gh`` CLI for one repository."""

    def __init__(self, repo: str, token: str, dry_run: bool = False) -> None:
        self.repo = repo
        self.dry_run = dry_run
        # Without an explicit token, fall back to whatever `gh` is authenticated
        # as — which is what makes a local `--dry-run` work.
        self.env = {**os.environ, "GH_TOKEN": token} if token else dict(os.environ)

    def run(self, *args: str, check: bool = True) -> subprocess.CompletedProcess:
        return _run(["gh", *args], env=self.env, check=check)

    def releases(self, upstream_repo: str, limit: int = 30) -> list[dict]:
        out = self.run(
            "release", "list", "--repo", upstream_repo, "--limit", str(limit),
            "--json", "tagName,name,isDraft,isPrerelease,publishedAt",
        ).stdout
        items = json.loads(out)
        # `gh release list` omits bodies; fetch them per release.
        for item in items:
            if item.get("isDraft") or item.get("isPrerelease"):
                continue
            detail = json.loads(
                self.run(
                    "release", "view", item["tagName"], "--repo", upstream_repo,
                    "--json", "body,url",
                ).stdout
            )
            item["body"] = detail.get("body") or ""
            item["url"] = detail.get("url") or ""
        return items

    def ensure_label(self, label: str, description: str, color: str) -> None:
        if self.dry_run:
            return
        self.run(
            "label", "create", label, "--repo", self.repo,
            "--description", description, "--color", color, "--force",
        )

    def find_watch_issue(self, label: str) -> dict | None:
        out = self.run(
            "issue", "list", "--repo", self.repo, "--state", "open",
            "--label", label, "--limit", "50", "--json", "number,title,body",
        ).stdout
        for issue in json.loads(out):
            if ISSUE_MARKER in (issue.get("body") or ""):
                return issue
        return None

    def find_refresh_pr(self, branch: str) -> dict | None:
        out = self.run(
            "pr", "list", "--repo", self.repo, "--state", "open", "--head", branch,
            "--limit", "10", "--json", "number,headRefName,headRefOid",
        ).stdout
        prs = json.loads(out)
        return prs[0] if prs else None


# --------------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------------- #


def read_pin(repo_root: Path) -> Pin:
    """Read the vendor pin, tolerating a missing ``PATCHES`` file."""
    commit_path = repo_root / COMMIT_FILE
    if not commit_path.is_file():
        raise WatchError(f"missing {COMMIT_FILE}")
    base = commit_path.read_text().strip()
    if not re.fullmatch(r"[0-9a-f]{7,40}", base):
        raise WatchError(f"{COMMIT_FILE} does not contain a commit sha: {base!r}")
    patches_path = repo_root / PATCHES_FILE
    text = patches_path.read_text() if patches_path.is_file() else ""
    return Pin(base=base, patches=parse_patches(text))


def clone_sources(upstream_repo: str, fork_repo: str, dest: Path) -> Git:
    """Clone upstream and add the fork as a second remote.

    Mirrors ``scripts/refresh-minibwa.sh``: the pin may reference a fork commit
    that upstream alone cannot resolve.
    """
    _run(["git", "clone", "--quiet", f"https://github.com/{upstream_repo}", str(dest)])
    git = Git(dest)
    git.run("remote", "add", "fork", f"https://github.com/{fork_repo}")
    git.run("fetch", "--quiet", "--tags", "fork")
    git.run("config", "user.name", "minibwa-upstream-watch")
    git.run("config", "user.email", "minibwa-upstream-watch@invalid")
    return git


def pin_tip(git: Git, pin: Pin) -> tuple[str, dict[str, Patch]]:
    """Reconstruct the commit the vendored tree corresponds to.

    With no carried patches the base *is* the tip. Otherwise the patches are
    replayed onto the base — the same picks that produced the vendored tree, onto
    the same base, so this is deterministic and cannot conflict.

    Returns the tip and a map from each *replayed* commit back to the original
    fork commit it reproduces. Replaying mints new shas that exist only in this
    clone, so callers must translate through this map before recording a sha
    anywhere durable.
    """
    if not pin.patches:
        return pin.base, {}
    git.run("checkout", "--quiet", "--detach", pin.base)
    replayed: dict[str, Patch] = {}
    for patch in pin.patches:
        result = git.run("cherry-pick", patch.sha, check=False)
        if result.returncode != 0:
            git.run("cherry-pick", "--abort", check=False)
            raise WatchError(
                f"replaying recorded patch {patch.sha} onto {pin.base} failed; "
                f"{PATCHES_FILE} and {COMMIT_FILE} disagree with the vendored tree"
            )
        replayed[git.out("rev-parse", "HEAD")] = patch
    return git.out("rev-parse", "HEAD"), replayed


def build_plan(git: Git, pin: Pin, releases: Sequence[Release]) -> Plan:
    """Decide where we stand relative to upstream and what can be prepared."""
    if not releases:
        raise WatchError("upstream has no published releases")

    resolved = [replace(r, commit=git.resolve_tag_commit(r.tag)) for r in releases]
    # A release whose tag the clone cannot resolve must not be silently skipped:
    # dropping the newest one would quietly compare against the wrong release.
    if not resolved[0].commit:
        raise WatchError(
            f"cannot resolve tag {releases[0].tag} in the clone; refusing to "
            "compare against an older release"
        )
    resolved = [r for r in resolved if r.commit]
    latest = resolved[0]

    tip, replayed = pin_tip(git, pin)
    current = newest_ancestor_release(resolved, tip, git.is_ancestor)
    plan = Plan(
        latest=latest,
        current=current,
        intervening=releases_between(resolved, current),
    )
    if not plan.behind:
        return plan

    # Everything in the pin that upstream still lacks, by patch id — so a fork
    # patch upstream merged in squashed form is recognised and not re-applied.
    # Translate each survivor back to the original fork commit: `tip` may be a
    # replay whose shas exist only in this clone, and recording those in PATCHES
    # would leave a pin nobody can resolve.
    carried = tuple(replayed.get(p.sha, p) for p in git.cherry(latest.commit, tip))
    # Compare full shas: PATCHES may record an abbreviation, which would never
    # match `git cherry`'s 40-character output.
    carried_full = {git.out("rev-parse", p.sha) for p in carried}
    upstreamed = tuple(
        p for p in pin.patches if git.out("rev-parse", p.sha) not in carried_full
    )

    git.run("checkout", "--quiet", "--detach", latest.commit)
    for patch in carried:
        result = git.run("cherry-pick", patch.sha, check=False)
        if result.returncode != 0:
            git.run("cherry-pick", "--abort", check=False)
            return replace(plan, carried=tuple(carried), upstreamed=upstreamed,
                           conflict=patch)
    return replace(
        plan,
        carried=tuple(carried),
        upstreamed=upstreamed,
        candidate=git.out("rev-parse", "HEAD"),
    )


def vendor_refresh(repo_root: Path, source: Path, plan: Plan) -> None:
    """Run the refresh script against the new release plus the carried patches.

    The script is given the release tag as the base and re-derives the tree by
    cherry-picking the patches itself, so ``COMMIT``/``PATCHES`` record a pin any
    human can reproduce — rather than the candidate SHA, which exists only in this
    runner's temporary clone.
    """
    cmd = [str(repo_root / REFRESH_SCRIPT), plan.latest.commit, "--src", str(source)]
    for patch in plan.carried:
        cmd += ["--patch", patch.sha]
    _run(cmd, cwd=repo_root)


def branch_is_bot_owned(repo: Git, branch: str, base: str, bot_email: str) -> bool:
    """True if every commit the branch adds on top of ``base`` is the bot's.

    Guards against force-pushing over a human's follow-up fix — a likely event,
    since a C API bump often needs a shim change on the same branch.
    """
    result = repo.run("rev-list", f"{base}..origin/{branch}", "--format=%ae",
                      "--no-commit-header", check=False)
    if result.returncode != 0:
        # A missing or unreadable remote-tracking ref tells us nothing about who
        # owns the branch, so it must not authorize a force-push. Fail closed.
        return False
    authors = {line.strip() for line in result.stdout.splitlines() if line.strip()}
    return authors.issubset({bot_email})


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True, help="owner/name of this repository")
    parser.add_argument("--upstream-repo", default="lh3/minibwa")
    parser.add_argument("--fork-repo", default="nh13/minibwa")
    parser.add_argument("--label", default="upstream-minibwa")
    parser.add_argument("--pr-label", default="dependencies")
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--branch", default=BRANCH)
    parser.add_argument("--bot-email", default="")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--output", default="", help="path to write GitHub Action outputs")
    return parser.parse_args(argv)


def emit_outputs(path: str, **values: str) -> None:
    if not path:
        return
    with open(path, "a", encoding="utf-8") as handle:
        for key, value in values.items():
            handle.write(f"{key}={value}\n")


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    token = os.environ.get("GH_TOKEN") or os.environ.get("GITHUB_TOKEN") or ""
    if not token and not args.dry_run:
        raise WatchError("GH_TOKEN is required unless --dry-run is set")

    repo_root = Path(args.repo_root).resolve()
    gh = GitHub(args.repo, token, dry_run=args.dry_run)
    pin = read_pin(repo_root)

    with tempfile.TemporaryDirectory() as tmp:
        source = Path(tmp) / "minibwa"
        git = clone_sources(args.upstream_repo, args.fork_repo, source)
        releases = published_releases(gh.releases(args.upstream_repo))
        plan = build_plan(git, pin, releases)

        if not plan.behind:
            print(f"vendored minibwa is current at {plan.latest.label}")
            _close_stale_issue(gh, args, plan)
            _close_stale_pull_request(gh, args, plan)
            emit_outputs(args.output, behind="false",
                         **{"current-release": plan.latest.tag,
                            "latest-release": plan.latest.tag,
                            "cherry-pick-status": "none"})
            return 0

        status = "conflict" if plan.blocked else ("clean" if plan.carried else "none")
        print(f"behind: {plan.current.label if plan.current else 'unknown'} -> "
              f"{plan.latest.label} (cherry-pick: {status})")

        if args.dry_run:
            print("\n===== ISSUE =====")
            print(issue_title(plan))
            print(render_issue_body(plan, pr_number=None))
            if not plan.blocked:
                print("\n===== PULL REQUEST =====")
                print(pr_title(plan))
                print(render_pr_body(plan, issue_number=None))
            emit_outputs(args.output, behind="true",
                         **{"current-release": plan.current.tag if plan.current else "",
                            "latest-release": plan.latest.tag,
                            "cherry-pick-status": status})
            return 0

        pr_number, pr_refreshed = None, False
        if not plan.blocked:
            pr_number, pr_refreshed = _prepare_pull_request(
                gh, git, repo_root, source, plan, args
            )
        stale_pr = None
        if plan.blocked:
            existing = gh.find_refresh_pr(args.branch)
            stale_pr = existing["number"] if existing else None
        issue_number = _file_or_update_issue(gh, args, plan, pr_number, stale_pr)
        # Only restate the PR body when this run actually re-vendored: on a
        # skipped refresh the description would claim a bump that never happened
        # and discard any human edits.
        if pr_refreshed and pr_number is not None and issue_number is not None:
            _sync_pr_body(gh, args, plan, pr_number, issue_number)

        emit_outputs(args.output, behind="true",
                     **{"current-release": plan.current.tag if plan.current else "",
                        "latest-release": plan.latest.tag,
                        "cherry-pick-status": status,
                        "issue-number": str(issue_number or ""),
                        "pr-number": str(pr_number or "")})
    return 0


def _close_stale_issue(gh: GitHub, args: argparse.Namespace, plan: Plan) -> None:
    """Resolve the tracking issue once the refresh has landed."""
    issue = gh.find_watch_issue(args.label)
    if issue is None or gh.dry_run:
        return
    gh.run("issue", "comment", str(issue["number"]), "--repo", gh.repo, "--body",
           f"Vendored minibwa is now current at **{plan.latest.label}** — closing.")
    gh.run("issue", "close", str(issue["number"]), "--repo", gh.repo)


def _close_stale_pull_request(gh: GitHub, args: argparse.Namespace, plan: Plan) -> None:
    """Retire the refresh PR once the pin has caught up.

    The pin can catch up without this PR merging — a maintainer may refresh by
    hand on their own branch — which would otherwise leave a draft PR open
    forever holding a redundant bump against a moved ``main``.
    """
    pull = gh.find_refresh_pr(args.branch)
    if pull is None or gh.dry_run:
        return
    gh.run("pr", "comment", str(pull["number"]), "--repo", gh.repo, "--body",
           f"Vendored minibwa is now current at **{plan.latest.label}** — "
           "this refresh is redundant, closing.")
    gh.run("pr", "close", str(pull["number"]), "--repo", gh.repo)


def _file_or_update_issue(gh: GitHub, args: argparse.Namespace, plan: Plan,
                          pr_number: int | None, stale_pr: int | None) -> int | None:
    """Create the tracking issue, or rewrite it in place and comment."""
    gh.ensure_label(args.label, "Vendored upstream minibwa is behind a release", "0366d6")
    body = render_issue_body(plan, pr_number, stale_pr)
    title = issue_title(plan)
    existing = gh.find_watch_issue(args.label)

    if existing is None:
        out = gh.run("issue", "create", "--repo", gh.repo, "--title", title,
                     "--body", body, "--label", args.label).stdout.strip()
        match = re.search(r"/issues/(\d+)", out)
        return int(match.group(1)) if match else None

    number = existing["number"]
    gh.run("issue", "edit", str(number), "--repo", gh.repo,
           "--title", title, "--body", body)
    if existing.get("title") != title:
        note = (f"Upstream moved on: now tracking **{plan.latest.label}**. "
                "Title and description updated above.")
        if plan.blocked:
            note += " The automated refresh is blocked by a cherry-pick conflict."
        elif pr_number is not None:
            note += f" Refresh PR: #{pr_number}."
        gh.run("issue", "comment", str(number), "--repo", gh.repo, "--body", note)
    return number


def _prepare_pull_request(gh: GitHub, git: Git, repo_root: Path, source: Path,
                          plan: Plan, args: argparse.Namespace) -> tuple[int | None, bool]:
    """Vendor the candidate tree, push the branch, and open or update the PR.

    Returns the PR number and whether this run actually refreshed it. A skipped
    refresh still yields the existing number so the issue can link it, but the
    caller must not then rewrite the PR body to claim a re-vendor that never
    happened.
    """
    repo = Git(repo_root)
    existing = gh.find_refresh_pr(args.branch)

    # `git fetch origin <branch>` alone leaves FETCH_HEAD only: actions/checkout
    # configures remote.origin.fetch for the checked-out branch, so the refspec
    # must be explicit or refs/remotes/origin/<branch> never exists — and both
    # the ownership guard and --force-with-lease depend on it.
    repo.run("fetch", "--quiet", "origin",
             f"+refs/heads/{args.branch}:refs/remotes/origin/{args.branch}", check=False)
    probe = repo.run("rev-parse", "--verify", "--quiet",
                     f"refs/remotes/origin/{args.branch}", check=False)
    remote_sha = probe.stdout.strip() if probe.returncode == 0 else ""

    if remote_sha and args.bot_email:
        if not branch_is_bot_owned(repo, args.branch, "origin/main", args.bot_email):
            print(f"branch {args.branch} carries non-bot commits; refusing to force-push")
            return (existing["number"] if existing else None), False

    repo.run("checkout", "--quiet", "-B", args.branch, "origin/main")
    vendor_refresh(repo_root, source, plan)
    repo.run("add", "--", str(COMMIT_FILE), str(PATCHES_FILE),
             "minibwa-sys/vendor/minibwa")
    if not repo.run("diff", "--cached", "--quiet", check=False).returncode:
        print("vendored tree already matches the candidate; nothing to push")
        return (existing["number"] if existing else None), False

    repo.run("commit", "--quiet", "-m", commit_subject(plan))
    if remote_sha:
        # An explicit lease: the bare form has nothing to compare against on a
        # branch whose remote-tracking ref we just created.
        repo.run("push", f"--force-with-lease={args.branch}:{remote_sha}",
                 "origin", args.branch)
    else:
        repo.run("push", "origin", args.branch)

    if existing is not None:
        gh.run("pr", "edit", str(existing["number"]), "--repo", gh.repo,
               "--title", pr_title(plan))
        return existing["number"], True

    # `gh pr create --label` fails outright on an unknown label, and this runs
    # before any issue-label setup, so ensure it here.
    gh.ensure_label(args.pr_label, "Pull requests that update a dependency file",
                    "0366d6")
    out = gh.run("pr", "create", "--repo", gh.repo, "--draft", "--base", "main",
                 "--head", args.branch, "--title", pr_title(plan),
                 "--body", render_pr_body(plan, issue_number=None),
                 "--label", args.pr_label).stdout.strip()
    match = re.search(r"/pull/(\d+)", out)
    return (int(match.group(1)) if match else None), True


def _sync_pr_body(gh: GitHub, args: argparse.Namespace, plan: Plan,
                  pr_number: int, issue_number: int) -> None:
    """Rewrite the PR body once the issue number is known, so they cross-link."""
    gh.run("pr", "edit", str(pr_number), "--repo", gh.repo,
           "--body", render_pr_body(plan, issue_number))


if __name__ == "__main__":
    try:
        sys.exit(main())
    except WatchError as exc:
        print(f"error: {exc}", file=sys.stderr)
        sys.exit(1)
