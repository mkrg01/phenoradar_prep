"""Release planning and failure recovery using local Git and substituted services."""
import importlib.util
import json
from pathlib import Path
import subprocess
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("release_automation", ROOT / ".github/scripts/release.py")
release = importlib.util.module_from_spec(spec)
spec.loader.exec_module(release)
REPOSITORY = "mkrg01/phenoradar_prep"
IMAGE = release.IMAGE_REPOSITORY
DIGEST = "sha256:" + "a" * 64
REBUILT_DIGEST = "sha256:" + "b" * 64


@pytest.fixture
def release_repo(tmp_path):
    remote = tmp_path / "remote.git"
    root = tmp_path / "checkout"
    root.mkdir()
    release.command(["git", "init", "--bare", str(remote)], root)
    release.git("init", "-b", "main", root=root)
    for key, value in [("user.name", "Release test"), ("user.email", "test@example.invalid"),
                       ("commit.gpgsign", "false"), ("tag.gpgsign", "false")]:
        release.git("config", key, value, root=root)
    (root / "VERSION").write_text("0.1.0\n")
    release.git("add", "VERSION", root=root)
    release.git("commit", "-m", "Initial release", root=root)
    release.git("tag", "v0.1.0", root=root)
    (root / "VERSION").write_text("0.2.0\n")
    release.git("commit", "-am", "Release 0.2.0", root=root)
    release.git("remote", "add", "origin", str(remote), root=root)
    release.git("push", "origin", "main", "--tags", root=root)
    return root, remote


def test_plan_selects_new_version_and_skips_ordinary_commits(release_repo):
    root, _ = release_repo
    assert release.plan(REPOSITORY, root)["tag"] == "v0.2.0"
    assert release.plan(REPOSITORY, root)["release"] == "true"
    (root / "VERSION").write_text("0.1.0\n")
    release.git("commit", "-am", "Continue work on the current version", root=root)
    assert release.plan(REPOSITORY, root)["release"] == "false"


def test_plan_rejects_version_regression_and_wrong_registry(release_repo):
    root, _ = release_repo
    (root / "VERSION").write_text("0.0.9\n")
    release.git("commit", "-am", "Invalid version regression", root=root)
    with pytest.raises(ValueError, match="newer"):
        release.plan(REPOSITORY, root)
    with pytest.raises(ValueError, match="For a fork"):
        release.plan("someone/fork", root)


def test_release_requires_committed_sources(release_repo):
    root, _ = release_repo
    (root / "VERSION").write_text("0.3.0\n")
    with pytest.raises(ValueError, match="Commit all release changes"):
        release.plan(REPOSITORY, root)


@pytest.fixture
def services(release_repo, monkeypatch):
    root, _ = release_repo
    real_command = release.command
    revision = release.git("rev-parse", "HEAD", root=root)
    state = SimpleNamespace(images={}, local={"phenoradar_prep:ci": DIGEST}, calls=[],
                            labels={}, public=True, draft=None, manifest=None,
                            fail_upload=False, fail_tag_push=False)
    labels = {"org.opencontainers.image.version": "0.2.0",
              "org.opencontainers.image.revision": revision}
    state.labels[DIGEST] = labels.copy()
    state.labels[REBUILT_DIGEST] = labels.copy()

    def fake(args, root=root, *, check=True):
        state.calls.append(list(args))
        stdout, stderr, code = "", "", 0
        if args[0] == "git":
            if args[1:3] == ["push", "origin"] and state.fail_tag_push:
                state.fail_tag_push = False
                raise RuntimeError("Simulated tag push failure")
            return real_command(args, root, check=check)
        if args[0] == "docker":
            if args[1] == "--config":
                if not state.public:
                    code, stderr = 1, "denied"
            elif args[1:3] == ["manifest", "inspect"]:
                digest = state.images.get(args[-1])
                if digest:
                    stdout = json.dumps({"Descriptor": {"digest": digest}})
                else:
                    code, stderr = 1, "no such manifest"
            elif args[1:3] == ["image", "inspect"]:
                stdout = json.dumps(state.labels[state.local[args[-1]]])
            elif args[1] == "tag":
                state.local[args[3]] = state.local[args[2]]
            elif args[1] == "push":
                state.images[args[2]] = state.local[args[2]]
            elif args[1] == "pull":
                state.local[args[2]] = args[2].split("@", 1)[1]
            elif args[1] != "run":
                raise AssertionError(args)
        elif args[:2] == ["gh", "api"]:
            if state.draft is None:
                code, stderr = 1, "gh: Not Found (HTTP 404)"
            else:
                stdout = json.dumps({"draft": state.draft})
        elif args[:3] == ["gh", "release", "create"]:
            state.draft = True
        elif args[:3] == ["gh", "release", "upload"]:
            if state.fail_upload:
                state.fail_upload = False
                raise RuntimeError("Simulated asset upload failure")
            state.manifest = json.loads(Path(args[4]).read_text())
        elif args[:3] == ["gh", "release", "edit"]:
            state.draft = False
        else:
            raise AssertionError(args)
        result = subprocess.CompletedProcess(args, code, stdout, stderr)
        if check and code:
            raise RuntimeError(stderr)
        return result

    monkeypatch.setattr(release, "command", fake)
    return state


def test_publish_only_pushes_tag_and_records_exact_image(release_repo, services):
    root, remote = release_repo
    before = release.git("rev-parse", "main", root=remote)
    release.publish(REPOSITORY, "phenoradar_prep:ci", root)
    assert release.git("rev-parse", "main", root=remote) == before
    assert release.git("rev-list", "-n", "1", "v0.2.0", root=remote) == before
    assert services.manifest == {"version": "0.2.0", "revision": before,
                                 "image": f"docker://{IMAGE}@{DIGEST}"}
    assert services.images[f"{IMAGE}:latest"] == DIGEST
    assert services.draft is False
    assert [args for args in services.calls if args[:2] == ["git", "push"]] == [
        ["git", "push", "origin", "refs/tags/v0.2.0"]]
    assert (root / "VERSION").read_text() == "0.2.0\n"
    assert release.git("status", "--porcelain", root=root) == ""
    assert release.plan(REPOSITORY, root)["release"] == "false"


@pytest.mark.parametrize("failure", ["fail_upload", "fail_tag_push"])
def test_retry_reuses_published_image_after_rebuild(release_repo, services, failure):
    root, _ = release_repo
    setattr(services, failure, True)
    with pytest.raises(RuntimeError, match="Simulated"):
        release.publish(REPOSITORY, "phenoradar_prep:ci", root)
    services.local["phenoradar_prep:ci"] = REBUILT_DIGEST
    release.publish(REPOSITORY, "phenoradar_prep:ci", root)
    assert services.images[f"{IMAGE}:v0.2.0"] == DIGEST
    assert services.images[f"{IMAGE}:latest"] == DIGEST
    assert services.manifest["image"] == f"docker://{IMAGE}@{DIGEST}"
    assert sum(args == ["docker", "push", f"{IMAGE}:v0.2.0"] for args in services.calls) == 1


def test_private_image_stops_before_tag_and_can_resume(release_repo, services):
    root, _ = release_repo
    services.public = False
    with pytest.raises(RuntimeError, match="visibility to public"):
        release.publish(REPOSITORY, "phenoradar_prep:ci", root)
    assert "v0.2.0" not in release.stable_tags(root)
    assert services.draft is None
    services.public = True
    release.publish(REPOSITORY, "phenoradar_prep:ci", root)
    assert services.draft is False


def test_resuming_older_draft_preserves_latest(release_repo, services):
    root, _ = release_repo
    services.fail_upload = True
    with pytest.raises(RuntimeError, match="Simulated"):
        release.publish(REPOSITORY, "phenoradar_prep:ci", root)
    (root / "VERSION").write_text("0.3.0\n")
    release.git("commit", "-am", "Release a newer version", root=root)
    release.git("tag", "v0.3.0", root=root)
    release.git("push", "origin", "main", "--tags", root=root)
    services.images[f"{IMAGE}:latest"] = REBUILT_DIGEST
    release.git("checkout", "v0.2.0", root=root)

    release.publish(REPOSITORY, "phenoradar_prep:ci", root)

    assert services.images[f"{IMAGE}:latest"] == REBUILT_DIGEST
    assert services.draft is False
    assert ["gh", "release", "edit", "v0.2.0", "--repo", REPOSITORY,
            "--draft=false", "--latest=false"] in services.calls


def test_conflicting_image_is_not_overwritten(release_repo, services):
    root, _ = release_repo
    services.images[f"{IMAGE}:v0.2.0"] = REBUILT_DIGEST
    services.labels[REBUILT_DIGEST]["org.opencontainers.image.revision"] = "another-commit"
    with pytest.raises(ValueError, match="refusing to overwrite"):
        release.publish(REPOSITORY, "phenoradar_prep:ci", root)
    assert "v0.2.0" not in release.stable_tags(root)
    assert ["docker", "push", f"{IMAGE}:v0.2.0"] not in services.calls


def test_registry_errors_are_not_treated_as_missing_images(monkeypatch):
    monkeypatch.setattr(release, "command", lambda *args, **kwargs:
                        subprocess.CompletedProcess(args, 1, "", "unauthorized: access denied"))
    with pytest.raises(RuntimeError, match="Cannot inspect"):
        release.image_digest(f"{IMAGE}:v0.2.0")


def test_release_api_errors_are_not_treated_as_missing_releases(monkeypatch):
    monkeypatch.setattr(release, "command", lambda *args, **kwargs:
                        subprocess.CompletedProcess(args, 1, "", "gh: forbidden (HTTP 403)"))
    with pytest.raises(RuntimeError, match="Cannot check"):
        release.release_info(REPOSITORY, "v0.2.0")
