"""Pinned dependency retrieval, including offline cache reuse and failed publication."""
import copy
import io
import json
import shutil
import subprocess
import tarfile
from pathlib import Path

import pytest

import dataset_software as software
from common import sha256

REVISION = "a" * 40
DIGEST = "b" * 64
VERSION = "1.2.3"


@pytest.fixture
def pins(tmp_path):
    return {"cache_dir": str(tmp_path / "software"), "version": VERSION,
            "revision": REVISION, "image_uri": "docker://ghcr.io/kfuku52/genegalleon@sha256:" + DIGEST,
            "repository": None, "image": None, "image_sha256": None, "settings": {}}


def archive(path, version=VERSION, unsafe=False):
    prefix = "genegalleon-" + REVISION + "/"
    with tarfile.open(path, "w:gz") as bundle:
        for name, content in (("VERSION", version + "\n"),
                              (software.ENTRYPOINT, "#!/bin/bash\nexit 0\n"),
                              ("workflow/support/defaults.tsv", "setting\tvalue\n")):
            data = content.encode()
            member = tarfile.TarInfo(prefix + name)
            member.size = len(data)
            member.mode = 0o755 if name.endswith(".sh") else 0o644
            bundle.addfile(member, io.BytesIO(data))
        if unsafe:
            member = tarfile.TarInfo(prefix + "../../escaped")
            member.size = 3
            bundle.addfile(member, io.BytesIO(b"bad"))


@pytest.fixture
def source_download(tmp_path, monkeypatch):
    source = tmp_path / "source.tar.gz"
    archive(source)
    calls = []
    def download(url, destination):
        calls.append(url)
        shutil.copyfile(source, destination)
    monkeypatch.setattr(software, "download", download)
    return source, calls


@pytest.fixture
def runtime(monkeypatch):
    calls = []
    state = {"fail": False, "version": VERSION, "revision": REVISION}
    monkeypatch.setattr(software.shutil, "which", lambda name: "/fake/apptainer" if name == "apptainer" else None)
    monkeypatch.setattr(software, "architecture", lambda: "amd64")
    def run(command, **kwargs):
        calls.append(command)
        assert command[:4] == ["/fake/apptainer", "pull", "--arch", "amd64"]
        image = Path(command[-2])
        image.write_bytes(b"synthetic SIF content")
        assert Path(kwargs["env"]["APPTAINER_TMPDIR"]).is_dir()
        if state["fail"]: raise subprocess.CalledProcessError(1, command)
        return subprocess.CompletedProcess(command, 0)
    def output(command, **kwargs):
        if command[1] == "inspect":
            return json.dumps({"data": {"attributes": {"labels": {
                "org.opencontainers.image.version": state["version"],
                "org.opencontainers.image.revision": state["revision"]}}}})
        assert command[1:] == ["--version"]
        return "apptainer test-runtime\n"
    monkeypatch.setattr(software.subprocess, "run", run)
    monkeypatch.setattr(software.subprocess, "check_output", output)
    return state, calls


@pytest.mark.parametrize("key,value,message", [
    ("version", "latest", "explicit"),
    ("revision", "main", "full 40-character"),
    ("image_uri", "docker://ghcr.io/kfuku52/genegalleon:latest", "OCI SHA256"),
    ("image_uri", "https://example.org/genegalleon.sif", "require"),
    ("image_sha256", "invalid", "checksum"),
])
def test_moving_or_unverifiable_sources_are_rejected(pins, key, value, message):
    pins[key] = value
    with pytest.raises(ValueError, match=message): software.validate(pins)


def test_download_source_once_then_reuse_offline_and_detect_changes(pins, source_download, monkeypatch):
    _, calls = source_download
    repository, receipt = software.fetch_source(pins)
    assert (repository / "VERSION").read_text().strip() == VERSION
    assert all(Path(r["path"]).is_file() for r in receipt["files"])
    assert len(calls) == 1
    monkeypatch.setattr(software, "download", lambda *args: pytest.fail("network access on cache hit"))
    assert software.fetch_source(pins) == (repository, receipt)
    # Include data/config files as well as executable code in the frozen receipts.
    (repository / "workflow/support/defaults.tsv").write_text("changed")
    with pytest.raises(ValueError, match="registered file changed"):
        software.fetch_source(pins)


def test_incorrect_source_version_does_not_publish_partial_cache(pins, source_download):
    source, calls = source_download
    archive(source, version="9.9.9")
    with pytest.raises(ValueError, match="version mismatch"): software.fetch_source(pins)
    assert not (Path(pins["cache_dir"]) / "source" / REVISION).exists()
    assert not list((Path(pins["cache_dir"]) / "source").glob(".fetch-*"))
    archive(source)
    software.fetch_source(pins)
    assert len(calls) == 2


def test_source_archive_cannot_escape_download_directory(pins, source_download):
    source, _ = source_download
    archive(source, unsafe=True)
    with pytest.raises(ValueError, match="unsafe GeneGalleon source archive"):
        software.fetch_source(pins)
    assert not list(Path(pins["cache_dir"]).rglob("escaped"))
    assert not (Path(pins["cache_dir"]) / "source" / REVISION).exists()


def test_oci_is_pulled_by_digest_and_cached_without_runtime_or_network(pins, runtime, monkeypatch):
    _, calls = runtime
    image, receipt = software.fetch_image(pins)
    assert calls[0][-1] == pins["image_uri"]
    assert receipt["files"][0]["sha256"] == sha256(image)
    assert receipt["identity"]["architecture"] == "amd64"
    monkeypatch.setattr(software.shutil, "which", lambda name: None)
    assert software.fetch_image(pins) == (image, receipt)
    assert len(calls) == 1
    image.write_bytes(b"corrupt")
    with pytest.raises(ValueError, match="registered file changed"): software.fetch_image(pins)


def test_failed_oci_pull_is_not_published_and_can_be_retried(pins, runtime):
    state, calls = runtime
    state["fail"] = True
    with pytest.raises(subprocess.CalledProcessError): software.fetch_image(pins)
    assert not list((Path(pins["cache_dir"]) / "images").glob("*/receipt.json"))
    assert not list((Path(pins["cache_dir"]) / "images").glob(".fetch-*"))
    state["fail"] = False
    image, _ = software.fetch_image(pins)
    assert image.exists() and len(calls) == 2


def test_image_source_mismatch_is_rejected_before_and_after_caching(pins, runtime):
    state, _ = runtime
    state["revision"] = "c" * 40
    with pytest.raises(ValueError, match="image/source mismatch"): software.fetch_image(pins)
    state["revision"] = REVISION
    software.fetch_image(pins)
    other = dict(pins, revision="d" * 40)
    with pytest.raises(ValueError, match="image/source mismatch"): software.fetch_image(other)


def test_direct_sif_requires_matching_file_checksum(pins, runtime, monkeypatch):
    pins.update(image_uri="https://example.org/genegalleon.sif", image_sha256="0" * 64)
    monkeypatch.setattr(software, "download", lambda url, path: Path(path).write_bytes(b"direct SIF"))
    with pytest.raises(ValueError, match="SIF SHA256 mismatch"): software.fetch_image(pins)
    import hashlib
    pins["image_sha256"] = hashlib.sha256(b"direct SIF").hexdigest()
    image, receipt = software.fetch_image(pins)
    assert sha256(image) == pins["image_sha256"]
    assert receipt["identity"]["expected_sha256"] == pins["image_sha256"]


def test_resolved_paths_and_receipts_bind_source_and_image(pins, source_download, runtime):
    original = copy.deepcopy(pins)
    resolved, records, lock = software.resolve(pins)
    assert pins == original
    assert Path(resolved["repository"], software.ENTRYPOINT).is_file()
    assert Path(resolved["image"]).is_file()
    assert lock["source"]["identity"]["revision"] == REVISION
    assert lock["container"]["identity"]["uri"] == pins["image_uri"]
    assert len(records) == 4  # Entry point, defaults, VERSION, and SIF.
    assert software.resolve(pins) == (resolved, records, lock)
