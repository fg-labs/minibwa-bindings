#!/usr/bin/env python3
"""Tests for the upstream-minibwa watch action.

Fixtures are built programmatically: the pure helpers are exercised with
fabricated values, and the git-dependent logic against throwaway repositories
created in a temp directory, so nothing is committed as test data.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

import watch
from watch import (
    Git,
    Patch,
    Pin,
    Plan,
    Release,
    WatchError,
    branch_is_bot_owned,
    build_plan,
    newest_ancestor_release,
    parse_patches,
    published_releases,
    releases_between,
    render_issue_body,
    render_pr_body,
    pin_tip,
    reproduce_command,
)

REPO_ROOT = Path(__file__).resolve().parents[3]


def _git(root: Path, *args: str) -> str:
    env = {
        **os.environ,
        "GIT_AUTHOR_NAME": "Test",
        "GIT_AUTHOR_EMAIL": "test@invalid",
        "GIT_COMMITTER_NAME": "Test",
        "GIT_COMMITTER_EMAIL": "test@invalid",
        "GIT_CONFIG_GLOBAL": "/dev/null",
        "GIT_CONFIG_SYSTEM": "/dev/null",
    }
    result = subprocess.run(
        ["git", *args], cwd=root, env=env, capture_output=True, text=True
    )
    if result.returncode != 0:
        raise AssertionError(f"git {' '.join(args)} failed: {result.stderr}")
    return result.stdout.strip()


class GitFixture:
    """A throwaway git repository with helpers for building history."""

    def __init__(self, root: Path) -> None:
        self.root = root
        root.mkdir(parents=True, exist_ok=True)
        _git(root, "init", "--quiet", "--initial-branch", "master")

    def write(self, name: str, content: str) -> None:
        (self.root / name).write_text(content)

    def commit(self, message: str, **files: str) -> str:
        for name, content in files.items():
            self.write(name, content)
        _git(self.root, "add", "-A")
        _git(self.root, "commit", "--quiet", "-m", message)
        return _git(self.root, "rev-parse", "HEAD")

    def tag(self, name: str) -> str:
        _git(self.root, "tag", name)
        return _git(self.root, "rev-parse", name)

    def checkout(self, ref: str) -> None:
        _git(self.root, "-c", "advice.detachedHead=false", "checkout", "--quiet", ref)

    def git(self) -> Git:
        gitwrap = Git(self.root)
        gitwrap.run("config", "user.name", "Test")
        gitwrap.run("config", "user.email", "test@invalid")
        return gitwrap


def release(tag: str, name: str = "", published: str = "") -> Release:
    return Release(
        tag=tag,
        name=name or f"Minibwa-{tag.lstrip('v')} (r{100 + int(float(tag.lstrip('v')) * 10)})",
        body=f"notes for {tag}",
        published_at=published or f"2026-0{tag[-1]}-01T00:00:00Z",
        url=f"https://example.invalid/{tag}",
    )


# --------------------------------------------------------------------------- #
# Pure helpers
# --------------------------------------------------------------------------- #


class ParsePatchesTest(unittest.TestCase):
    def test_parses_sha_and_subject(self) -> None:
        text = "abc1234  add AVX2 dispatch\ndef5678  another patch\n"
        self.assertEqual(
            parse_patches(text),
            (Patch("abc1234", "add AVX2 dispatch"), Patch("def5678", "another patch")),
        )

    def test_ignores_comments_and_blanks(self) -> None:
        text = "# a comment\n\n   \nabc1234  subject\n"
        self.assertEqual(parse_patches(text), (Patch("abc1234", "subject"),))

    def test_empty_file_means_no_patches(self) -> None:
        self.assertEqual(parse_patches(""), ())

    def test_rejects_malformed_line(self) -> None:
        with self.assertRaises(WatchError):
            parse_patches("not-a-sha  subject\n")

    def test_parses_the_header_refresh_minibwa_writes(self) -> None:
        """The live writer is the heredoc in scripts/refresh-minibwa.sh."""
        text = (
            "# Fork commits cherry-picked onto the upstream base recorded in COMMIT,\n"
            "# in apply order. Reproduce the vendored tree with:\n"
            "#\n"
            f"#   scripts/refresh-minibwa.sh {'a' * 40} --patch {'b' * 40}\n"
            "#\n"
            f"{'b' * 40}  add AVX2 dispatch\n"
        )
        self.assertEqual(parse_patches(text), (Patch("b" * 40, "add AVX2 dispatch"),))


class ReleaseSelectionTest(unittest.TestCase):
    def test_filters_drafts_and_prereleases(self) -> None:
        payload = [
            {"tagName": "v0.7", "name": "Minibwa-0.7 (r420)", "isDraft": True,
             "isPrerelease": False, "publishedAt": "2026-08-01T00:00:00Z"},
            {"tagName": "v0.6rc", "name": "rc", "isDraft": False,
             "isPrerelease": True, "publishedAt": "2026-07-29T00:00:00Z"},
            {"tagName": "v0.6", "name": "Minibwa-0.6 (r416)", "isDraft": False,
             "isPrerelease": False, "publishedAt": "2026-07-30T00:00:00Z"},
        ]
        selected = published_releases(payload)
        self.assertEqual([r.tag for r in selected], ["v0.6"])

    def test_sorts_newest_first(self) -> None:
        payload = [
            {"tagName": "v0.4", "isDraft": False, "isPrerelease": False,
             "publishedAt": "2026-07-12T00:00:00Z"},
            {"tagName": "v0.6", "isDraft": False, "isPrerelease": False,
             "publishedAt": "2026-07-30T00:00:00Z"},
            {"tagName": "v0.5", "isDraft": False, "isPrerelease": False,
             "publishedAt": "2026-07-26T00:00:00Z"},
        ]
        self.assertEqual([r.tag for r in published_releases(payload)],
                         ["v0.6", "v0.5", "v0.4"])


class NewestAncestorTest(unittest.TestCase):
    def test_picks_newest_ancestor(self) -> None:
        releases = [
            Release("v0.6", "", "", "", "", commit="c6"),
            Release("v0.5", "", "", "", "", commit="c5"),
            Release("v0.4", "", "", "", "", commit="c4"),
        ]
        ancestors = {"c4"}
        found = newest_ancestor_release(
            releases, "pin", lambda a, _h: a in ancestors
        )
        self.assertIsNotNone(found)
        self.assertEqual(found.tag, "v0.4")

    def test_returns_none_when_no_release_is_an_ancestor(self) -> None:
        releases = [Release("v0.6", "", "", "", "", commit="c6")]
        self.assertIsNone(newest_ancestor_release(releases, "pin", lambda *_: False))


class ReleasesBetweenTest(unittest.TestCase):
    def test_returns_releases_newer_than_current(self) -> None:
        releases = [release("v0.6"), release("v0.5"), release("v0.4")]
        between = releases_between(releases, releases[2])
        self.assertEqual([r.tag for r in between], ["v0.6", "v0.5"])

    def test_unknown_current_yields_all(self) -> None:
        releases = [release("v0.6"), release("v0.5")]
        self.assertEqual(len(releases_between(releases, None)), 2)


class ReleaseLabelTest(unittest.TestCase):
    def test_parses_upstream_naming(self) -> None:
        rel = Release("v0.6", "Minibwa-0.6 (r416)", "", "", "")
        self.assertEqual(rel.label, "0.6 (r416)")
        self.assertEqual(rel.version, "0.6-r416")

    def test_falls_back_to_tag(self) -> None:
        rel = Release("v0.6", "some other name", "", "", "")
        self.assertEqual(rel.label, "v0.6")
        self.assertEqual(rel.version, "0.6")


# --------------------------------------------------------------------------- #
# Rendering
# --------------------------------------------------------------------------- #


def sample_plan(**overrides) -> Plan:
    latest = Release("v0.6", "Minibwa-0.6 (r416)", "fixed a thing",
                     "2026-07-30T00:00:00Z", "https://example.invalid/v0.6",
                     commit="f" * 40)
    current = Release("v0.4", "Minibwa-0.4 (r400)", "old notes",
                      "2026-07-12T00:00:00Z", "https://example.invalid/v0.4",
                      commit="a" * 40)
    defaults = dict(
        latest=latest,
        current=current,
        intervening=(latest,),
        carried=(Patch("b" * 40, "add AVX2/AVX-512 to ksw_extd2 with runtime dispatch"),),
        candidate="c" * 40,
    )
    defaults.update(overrides)
    return Plan(**defaults)


class RenderIssueTest(unittest.TestCase):
    def test_includes_dedup_marker(self) -> None:
        self.assertIn(watch.ISSUE_MARKER, render_issue_body(sample_plan(), None))

    def test_links_prepared_pr(self) -> None:
        body = render_issue_body(sample_plan(), pr_number=42)
        self.assertIn("#42", body)
        self.assertNotIn("Automated refresh is blocked", body)

    def test_conflict_is_loud_and_suppresses_pr_section(self) -> None:
        conflict = Patch("b" * 40, "add AVX2/AVX-512 to ksw_extd2 with runtime dispatch")
        body = render_issue_body(sample_plan(conflict=conflict, candidate=None), None)
        self.assertIn("Automated refresh is blocked", body)
        self.assertIn(conflict.sha[:8], body)
        self.assertNotIn("Prepared refresh", body)

    def test_reports_patches_that_landed_upstream(self) -> None:
        merged = Patch("d" * 40, "now merged upstream")
        body = render_issue_body(sample_plan(upstreamed=(merged,)), None)
        self.assertIn("merged upstream", body)
        self.assertIn(merged.sha[:8], body)

    def test_handles_unknown_current_release(self) -> None:
        body = render_issue_body(sample_plan(current=None, intervening=()), None)
        self.assertIn("unknown", body)


class RenderPullRequestTest(unittest.TestCase):
    def test_cross_links_issue_and_shows_reproduce_command(self) -> None:
        plan = sample_plan()
        body = render_pr_body(plan, issue_number=7)
        self.assertIn("#7", body)
        self.assertIn(reproduce_command(plan), body)

    def test_reproduce_command_lists_every_patch(self) -> None:
        plan = sample_plan(carried=(Patch("b" * 40, "one"), Patch("c" * 40, "two")))
        command = reproduce_command(plan)
        self.assertIn(plan.latest.commit, command)
        self.assertIn("b" * 40, command)
        self.assertIn("c" * 40, command)


# --------------------------------------------------------------------------- #
# git-backed behaviour
# --------------------------------------------------------------------------- #


class BuildPlanTest(unittest.TestCase):
    """Exercises ancestry, patch carrying, and cherry-pick outcomes on real repos."""

    def setUp(self) -> None:
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.fixture = GitFixture(Path(self.tmp.name) / "minibwa")
        self.fixture.commit("initial", **{"core.c": "line1\nline2\nline3\n"})
        self.v04 = self.fixture.tag("v0.4")

    def _add_upstream_release(self, tag: str, content: str) -> str:
        self.fixture.checkout("master")
        self.fixture.commit(f"release {tag}", **{"core.c": content})
        return self.fixture.tag(tag)

    def _add_fork_patch(self, content: str) -> str:
        """A patch branched off v0.4, mimicking the carried AVX2 commit."""
        self.fixture.checkout(self.v04)
        sha = self.fixture.commit("add AVX2 dispatch", **{"simd.c": content})
        _git(self.fixture.root, "branch", "-f", "forkwork", sha)
        return sha

    def test_up_to_date_pin_is_not_behind(self) -> None:
        releases = [release("v0.4")]
        plan = build_plan(self.fixture.git(), Pin(base=self.v04), releases)
        self.assertFalse(plan.behind)
        self.assertEqual(plan.current.tag, "v0.4")

    def test_behind_with_no_patches_needs_no_cherry_pick(self) -> None:
        self._add_upstream_release("v0.6", "line1\nline2 changed\nline3\n")
        releases = [release("v0.6"), release("v0.4")]
        plan = build_plan(self.fixture.git(), Pin(base=self.v04), releases)
        self.assertTrue(plan.behind)
        self.assertEqual(plan.current.tag, "v0.4")
        self.assertEqual(plan.carried, ())
        self.assertIsNone(plan.conflict)
        self.assertIsNotNone(plan.candidate)

    def test_clean_cherry_pick_carries_the_fork_patch(self) -> None:
        patch = self._add_fork_patch("avx2\n")
        self._add_upstream_release("v0.6", "line1\nline2 changed\nline3\n")
        releases = [release("v0.6"), release("v0.4")]
        plan = build_plan(self.fixture.git(), Pin(base=patch), releases)
        self.assertTrue(plan.behind)
        self.assertEqual([p.sha for p in plan.carried], [patch])
        self.assertIsNone(plan.conflict)
        self.assertIsNotNone(plan.candidate)

    def test_conflicting_cherry_pick_reports_the_offending_patch(self) -> None:
        self.fixture.checkout(self.v04)
        patch = self.fixture.commit("fork edit", **{"core.c": "line1\nFORK\nline3\n"})
        self._add_upstream_release("v0.6", "line1\nUPSTREAM\nline3\n")
        releases = [release("v0.6"), release("v0.4")]
        plan = build_plan(self.fixture.git(), Pin(base=patch), releases)
        self.assertTrue(plan.blocked)
        self.assertIsNotNone(plan.conflict)
        self.assertEqual(plan.conflict.sha, patch)
        self.assertIsNone(plan.candidate)

    def test_patch_merged_upstream_is_dropped_not_reapplied(self) -> None:
        patch = self._add_fork_patch("avx2\n")
        # Upstream lands the identical change as its own commit: same patch id,
        # different sha. `git cherry` must recognise it and not re-apply.
        self.fixture.checkout("master")
        self.fixture.commit("upstream adds avx2", **{"simd.c": "avx2\n"})
        self._add_upstream_release("v0.6", "line1\nline2 changed\nline3\n")
        releases = [release("v0.6"), release("v0.4")]
        plan = build_plan(
            self.fixture.git(),
            Pin(base=self.v04, patches=(Patch(patch, "add AVX2 dispatch"),)),
            releases,
        )
        self.assertTrue(plan.behind)
        self.assertEqual(plan.carried, ())
        self.assertEqual([p.sha for p in plan.upstreamed], [patch])
        self.assertIsNone(plan.conflict)

    def test_declared_patch_keeps_its_original_sha(self) -> None:
        """Replaying a pin mints new shas; they must never reach the plan.

        Recording a replayed sha in PATCHES would pin a commit that exists only
        in the runner's temp clone, and would also make the patch look merged
        upstream because it no longer matches anything in `pin.patches`.
        """
        patch = self._add_fork_patch("avx2\n")
        self._add_upstream_release("v0.6", "line1\nline2 changed\nline3\n")
        pin = Pin(base=self.v04, patches=(Patch(patch, "add AVX2 dispatch"),))
        plan = build_plan(self.fixture.git(), pin, [release("v0.6"), release("v0.4")])
        self.assertEqual([p.sha for p in plan.carried], [patch])
        self.assertEqual(plan.upstreamed, ())

    def test_abbreviated_declared_patch_is_not_reported_upstreamed(self) -> None:
        patch = self._add_fork_patch("avx2\n")
        self._add_upstream_release("v0.6", "line1\nline2 changed\nline3\n")
        pin = Pin(base=self.v04, patches=(Patch(patch[:8], "add AVX2 dispatch"),))
        plan = build_plan(self.fixture.git(), pin, [release("v0.6"), release("v0.4")])
        self.assertEqual(plan.upstreamed, ())
        self.assertEqual(len(plan.carried), 1)

    def test_pin_tip_replays_declared_patches(self) -> None:
        patch = self._add_fork_patch("avx2\n")
        git = self.fixture.git()
        tip, replayed = pin_tip(git, Pin(base=self.v04, patches=(Patch(patch, "add AVX2"),)))
        self.assertNotEqual(tip, self.v04)
        self.assertEqual(replayed, {tip: Patch(patch, "add AVX2")})
        self.assertTrue((self.fixture.root / "simd.c").exists())

    def test_pin_tip_without_patches_is_the_base(self) -> None:
        self.assertEqual(pin_tip(self.fixture.git(), Pin(base=self.v04)), (self.v04, {}))


class BranchOwnershipTest(unittest.TestCase):
    """The force-push guard that protects a human's follow-up commits."""

    def setUp(self) -> None:
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.fixture = GitFixture(Path(self.tmp.name) / "repo")
        self.fixture.commit("initial", **{"a.txt": "a\n"})
        _git(self.fixture.root, "branch", "main")

    def _branch_with_author(self, branch: str, email: str) -> None:
        _git(self.fixture.root, "checkout", "--quiet", "-B", branch, "main")
        (self.fixture.root / "b.txt").write_text("b\n")
        _git(self.fixture.root, "add", "-A")
        _git(self.fixture.root, "-c", f"user.email={email}", "-c", "user.name=X",
             "commit", "--quiet", "--author", f"X <{email}>", "-m", "change")
        # Stand in for the remote-tracking ref the real code inspects.
        _git(self.fixture.root, "update-ref", f"refs/remotes/origin/{branch}",
             _git(self.fixture.root, "rev-parse", "HEAD"))
        _git(self.fixture.root, "update-ref", "refs/remotes/origin/main",
             _git(self.fixture.root, "rev-parse", "main"))

    def test_bot_only_branch_is_owned(self) -> None:
        self._branch_with_author("bot-branch", "bot@invalid")
        self.assertTrue(
            branch_is_bot_owned(Git(self.fixture.root), "bot-branch",
                                "origin/main", "bot@invalid")
        )

    def test_missing_remote_ref_fails_closed(self) -> None:
        """No remote-tracking ref means no evidence — never authorize the push."""
        self.assertFalse(
            branch_is_bot_owned(Git(self.fixture.root), "never-fetched",
                                "origin/main", "bot@invalid")
        )

    def test_human_commit_blocks_force_push(self) -> None:
        self._branch_with_author("bot-branch", "human@invalid")
        self.assertFalse(
            branch_is_bot_owned(Git(self.fixture.root), "bot-branch",
                                "origin/main", "bot@invalid")
        )


class RefreshScriptTest(unittest.TestCase):
    """End-to-end check of scripts/refresh-minibwa.sh against a fake upstream."""

    def setUp(self) -> None:
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        root = Path(self.tmp.name)

        self.source = GitFixture(root / "minibwa")
        self.base = self.source.commit(
            "initial", **{"minibwa.h": "#define MB_VERSION \"0.4\"\n", "main.c": "int main(){}\n"}
        )
        self.patch = self.source.commit("add simd", **{"simd.c": "avx2\n"})

        # A throwaway copy of the project layout the script expects.
        self.project = root / "project"
        (self.project / "scripts").mkdir(parents=True)
        (self.project / "minibwa-sys" / "vendor").mkdir(parents=True)
        shutil.copy(REPO_ROOT / "scripts" / "refresh-minibwa.sh",
                    self.project / "scripts" / "refresh-minibwa.sh")

    def _run_refresh(self, *args: str) -> subprocess.CompletedProcess:
        return subprocess.run(
            ["bash", str(self.project / "scripts" / "refresh-minibwa.sh"), *args],
            capture_output=True, text=True,
        )

    def test_records_base_and_no_patches(self) -> None:
        result = self._run_refresh(self.base, "--src", str(self.source.root))
        self.assertEqual(result.returncode, 0, result.stderr)
        vendor = self.project / "minibwa-sys" / "vendor"
        self.assertEqual((vendor / "COMMIT").read_text().strip(), self.base)
        self.assertEqual(parse_patches((vendor / "PATCHES").read_text()), ())
        self.assertTrue((vendor / "minibwa" / "main.c").exists())

    def test_records_carried_patches(self) -> None:
        result = self._run_refresh(
            self.base, "--patch", self.patch, "--src", str(self.source.root)
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        vendor = self.project / "minibwa-sys" / "vendor"
        self.assertEqual((vendor / "COMMIT").read_text().strip(), self.base)
        patches = parse_patches((vendor / "PATCHES").read_text())
        self.assertEqual([p.sha for p in patches], [self.patch])
        self.assertEqual(patches[0].subject, "add simd")
        # The patch's content must actually be in the vendored tree.
        self.assertTrue((vendor / "minibwa" / "simd.c").exists())

    def test_already_applied_patch_is_not_called_a_conflict(self) -> None:
        """The common path once upstream merges a carried patch."""
        result = self._run_refresh(
            self.patch, "--patch", self.patch, "--src", str(self.source.root)
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("already contained", result.stderr)
        self.assertNotIn("conflicts", result.stderr)

    def test_genuine_conflict_is_reported_as_a_conflict(self) -> None:
        base = self.source.commit("upstream edit", **{"main.c": "int main(){return 1;}\n"})
        self.source.checkout(self.base)
        rival = self.source.commit("fork edit", **{"main.c": "int main(){return 2;}\n"})
        result = self._run_refresh(base, "--patch", rival, "--src", str(self.source.root))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("conflicts", result.stderr)
        self.assertIn(rival[:8], result.stderr)

    def test_rejects_unknown_flag(self) -> None:
        result = self._run_refresh(self.base, "--bogus")
        self.assertNotEqual(result.returncode, 0)

    def test_requires_a_commit(self) -> None:
        self.assertNotEqual(self._run_refresh().returncode, 0)


if __name__ == "__main__":
    unittest.main()
