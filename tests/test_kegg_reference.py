import io
import gzip
import json
import shutil
import subprocess
import sys
import tarfile
import urllib.error
from pathlib import Path

import pytest

from common import read_tsv, sha256
from prepare_kegg_reference import download_links, normalize_links, prepare
from verify_kegg_reference import verify
import bootstrap_kegg_reference as bootstrap_module


@pytest.fixture
def kegg_inputs(tmp_path):
    source = tmp_path / "source"
    profiles = source / "profiles"
    profiles.mkdir(parents=True)
    for ko in ["K00001", "K00002"]:
        (profiles / f"{ko}.hmm").write_text(f"HMMER3/f\nNAME  {ko}\nALPH  amino\nHMM A C D\n//\n")
    ko_list = source / "ko_list"
    ko_list.write_text(
        "knum\tthreshold\tscore_type\tprofile_type\tF-measure\tnseq\tnseq_used\talen\tmlen\teff_nseq\tre/pos\tdefinition\n"
        "K00001\t10\tfull\tall\t0.9\t3\t3\t2\t2\t3\t0.5\tfunction one\n"
        "K00002\t-\tdomain\tall\t0.9\t3\t3\t2\t2\t3\t0.5\tfunction two\n"
    )
    module = source / "module_links.tsv"
    module.write_text("ko:K00001\tmd:M00001\nmd:M00002\tko:K00001\nko:K00001\tmd:M00001\n")
    pathway = source / "pathway_links.tsv"
    pathway.write_text("ko:K00001\tpath:map00010\npath:ko00010\tko:K00001\nko:K00002\tpath:ko00020\n")
    return {"profiles_dir": profiles, "ko_list": ko_list,
            "module_links": module, "pathway_links": pathway,
            "reference_dir": tmp_path / "snapshot", "release": "offline-fixture"}


def test_snapshot_is_portable_independent_and_deduplicated(kegg_inputs, tmp_path):
    reference = prepare(**kegg_inputs)
    metadata = verify(reference)
    assert metadata["release"] == "offline-fixture"
    assert metadata["reference_id"] == sha256(reference)
    assert metadata["counts"]["profiles"] == 2
    assert read_tsv(metadata["ko_modules"]) == [
        {"ko": "K00001", "module": "M00001"}, {"ko": "K00001", "module": "M00002"}]
    assert read_tsv(metadata["ko_pathways"]) == [
        {"ko": "K00001", "pathway": "map00010"}, {"ko": "K00002", "pathway": "map00020"}]
    shutil.rmtree(kegg_inputs["profiles_dir"].parent)
    destination = tmp_path / "moved"
    reference.parent.rename(destination)
    moved = verify(destination / "reference.json", full=False)
    assert Path(moved["profiles_dir"]) == destination / "profiles"
    assert moved["reference_id"] == metadata["reference_id"]
    assert moved["verification_mode"] == "sizes"
    verify(destination / "reference.json")


@pytest.mark.parametrize("existing", ["empty", "nonempty", "symlink"])
def test_refuses_any_existing_destination(kegg_inputs, tmp_path, existing):
    root = kegg_inputs["reference_dir"]
    if existing == "symlink":
        root.symlink_to(tmp_path / "missing-target")
    else:
        root.mkdir()
        if existing == "nonempty":
            (root / "keep.txt").write_text("untouched")
    with pytest.raises(FileExistsError, match="already exists"):
        prepare(**kegg_inputs)
    if existing == "nonempty":
        assert (root / "keep.txt").read_text() == "untouched"
    if existing == "symlink":
        assert root.is_symlink()


def test_completed_snapshot_is_never_overwritten(kegg_inputs):
    reference = prepare(**kegg_inputs)
    before = reference.read_bytes()
    with pytest.raises(FileExistsError):
        prepare(**kegg_inputs)
    assert reference.read_bytes() == before
    verify(reference)


@pytest.mark.parametrize("mutation,error", [
    ("empty_profiles", "no extracted"), ("truncated_hmm", "incomplete"),
    ("wrong_name", "mismatched"), ("non_hmm", "HMMER3"),
    ("empty_list", "header"), ("short_list", "header"),
    ("missing_ko", "absent from ko_list"), ("duplicate_ko", "duplicate KO"),
    ("invalid_threshold", "threshold"), ("nan_threshold", "threshold"),
    ("invalid_type", "score_type"), ("duplicate_profile", "duplicate profile"),
])
def test_rejects_incomplete_or_inconsistent_kofam(kegg_inputs, mutation, error):
    profiles, ko_list = kegg_inputs["profiles_dir"], kegg_inputs["ko_list"]
    profile = profiles / "K00001.hmm"
    if mutation == "empty_profiles":
        for path in profiles.glob("*.hmm"):
            path.unlink()
    elif mutation == "truncated_hmm":
        profile.write_text(profile.read_text().replace("//\n", ""))
    elif mutation == "wrong_name":
        profile.write_text(profile.read_text().replace("NAME  K00001", "NAME  K00003"))
    elif mutation == "non_hmm":
        profile.write_text("an HTML error page")
    elif mutation == "empty_list":
        ko_list.write_text("")
    elif mutation == "short_list":
        ko_list.write_text("knum\tthreshold\tscore_type\nK00001\t10\tfull\n")
    elif mutation == "missing_ko":
        ko_list.write_text("\n".join(ko_list.read_text().splitlines()[:-1]) + "\n")
    elif mutation == "duplicate_ko":
        ko_list.write_text(ko_list.read_text() + ko_list.read_text().splitlines()[1] + "\n")
    elif mutation in {"invalid_threshold", "nan_threshold", "invalid_type"}:
        replacement = {"invalid_threshold": "bad", "nan_threshold": "nan", "invalid_type": "10"}[mutation]
        content = ko_list.read_text().replace("K00001\t10", f"K00001\t{replacement}")
        if mutation == "invalid_type":
            content = content.replace("\tfull\t", "\tunknown\t")
        ko_list.write_text(content)
    elif mutation == "duplicate_profile":
        (profiles / "nested").mkdir()
        shutil.copyfile(profile, profiles / "nested" / profile.name)
    with pytest.raises(ValueError, match=error):
        prepare(**kegg_inputs)
    assert not kegg_inputs["reference_dir"].exists()


@pytest.mark.parametrize("kind,content", [
    ("module", ""), ("module", "<html>service unavailable</html>\n"),
    ("module", "ko:KABC01\tmd:M00001\n"),
    ("module", "ko:K00001\tmd:hsa_M00001\n"),
    ("pathway", "ko:K00001\tpath:hsa00010\n"),
    ("pathway", "ko:K00001\tpath:map00010\textra\n"),
])
def test_rejects_bad_mappings_without_publishing(kegg_inputs, kind, content):
    kegg_inputs[f"{kind}_links"].write_text(content)
    with pytest.raises(ValueError, match="invalid|empty"):
        prepare(**kegg_inputs)
    assert not kegg_inputs["reference_dir"].exists()
    assert not list(kegg_inputs["reference_dir"].parent.glob(".snapshot.preparing-*"))


def test_same_size_corruption_requires_full_verification(kegg_inputs):
    reference = prepare(**kegg_inputs)
    profile = reference.parent / "profiles" / "K00001.hmm"
    profile.write_text(profile.read_text().replace("HMM A C D", "HMM A G D"))
    verify(reference, full=False)
    with pytest.raises(ValueError, match="checksum mismatch"):
        verify(reference)


@pytest.mark.parametrize("mutation", ["size", "missing", "extra", "inventory", "symlink"])
def test_quick_verification_rejects_invalid_inventory_or_files(kegg_inputs, mutation, tmp_path):
    reference = prepare(**kegg_inputs)
    profile = reference.parent / "profiles" / "K00001.hmm"
    if mutation == "size":
        profile.write_text(profile.read_text() + "\n")
    elif mutation == "missing":
        profile.unlink()
    elif mutation == "extra":
        (profile.parent / "K00003.hmm").write_text(profile.read_text())
    elif mutation == "inventory":
        inventory = reference.parent / "files.json"
        inventory.write_text(inventory.read_text() + " ")
    else:
        external = tmp_path / "external.hmm"
        profile.rename(external)
        profile.symlink_to(external)
    with pytest.raises(ValueError, match="mismatch|unrecorded|symlink"):
        verify(reference, full=False)


def test_cli_verification_report(kegg_inputs, tmp_path):
    reference = prepare(**kegg_inputs)
    script = Path(__file__).resolve().parents[1] / "workflow/scripts/verify_kegg_reference.py"
    output = tmp_path / "reference_qc.json"
    subprocess.run([sys.executable, str(script), "--reference", str(reference), "--output", str(output)], check=True)
    report = json.loads(output.read_text())
    assert report["reference_id"] == sha256(reference)
    assert report["verification_mode"] == "sha256"
    assert report["verified_files"] == 7


def test_download_retries_and_rate_limits(monkeypatch, tmp_path):
    import prepare_kegg_reference as module
    waits, calls = [], []
    monkeypatch.setattr(module.time, "sleep", waits.append)

    class Response(io.BytesIO):
        status = 200
        headers = {"ETag": "fixture"}

    def urlopen(request, timeout):
        calls.append(request.full_url)
        if len(calls) == 1:
            raise urllib.error.HTTPError(request.full_url, 429, "rate limited", {"Retry-After": "3"}, None)
        return Response(b"ko:K00001\tmd:M00001\n")

    monkeypatch.setattr(module.urllib.request, "urlopen", urlopen)
    destination = tmp_path / "module.tsv"
    source = download_links("https://rest.kegg.jp/link/module/ko", destination)
    assert source["etag"] == "fixture"
    assert waits == [0.5, 3, 0.5]
    assert len(calls) == 2
    assert normalize_links(destination, "module") == [{"ko": "K00001", "module": "M00001"}]


def test_download_does_not_retry_permanent_errors(monkeypatch, tmp_path):
    import prepare_kegg_reference as module
    monkeypatch.setattr(module.time, "sleep", lambda _: None)
    calls = []

    def urlopen(request, timeout):
        calls.append(request.full_url)
        raise urllib.error.HTTPError(request.full_url, 404, "not found", {}, None)

    monkeypatch.setattr(module.urllib.request, "urlopen", urlopen)
    with pytest.raises(urllib.error.HTTPError):
        download_links("https://rest.kegg.jp/link/module/ko", tmp_path / "module.tsv")
    assert len(calls) == 1


@pytest.fixture
def reference_downloads(kegg_inputs, monkeypatch):
    """Serve real compressed fixture files through the bootstrap HTTP boundary."""
    archive = io.BytesIO()
    with tarfile.open(fileobj=archive, mode="w:gz") as handle:
        handle.add(kegg_inputs["profiles_dir"], arcname="profiles")
    payloads = {"profiles": archive.getvalue(),
                "ko_list": gzip.compress(kegg_inputs["ko_list"].read_bytes()),
                "module_links": kegg_inputs["module_links"].read_bytes(),
                "pathway_links": kegg_inputs["pathway_links"].read_bytes()}
    urls = {key: url for key, (url, _) in bootstrap_module.DOWNLOADS.items()}
    calls = []

    class Response(io.BytesIO):
        status = 200

        def __init__(self, payload):
            super().__init__(payload)
            self.headers = {"Content-Length": str(len(payload)), "ETag": "fixture"}

    def fetch(request, timeout):
        calls.append(request.full_url)
        key = next(key for key, url in urls.items() if url == request.full_url)
        return Response(payloads[key])

    monkeypatch.setattr(bootstrap_module, "urlopen", fetch)
    monkeypatch.setattr(bootstrap_module.time, "sleep", lambda _: None)
    return payloads, urls, calls, fetch


def test_bootstrap_downloads_once_and_reuses_portable_snapshot(tmp_path, reference_downloads):
    _, urls, calls, _ = reference_downloads
    root = tmp_path / "resources/kegg/snapshot_v1"
    reference = bootstrap_module.bootstrap(root)
    assert calls == list(urls.values())
    assert verify(reference)["counts"]["profiles"] == 2
    metadata = json.loads(reference.read_text())
    for key, (_, filename) in bootstrap_module.DOWNLOADS.items():
        source = metadata["sources"]["downloads"][key]
        assert source["sha256"] == sha256(root.parent / "downloads" / filename)
        assert source["url"] == urls[key]
    before = {p: (sha256(p), p.stat().st_mtime_ns) for p in root.rglob("*") if p.is_file()}
    shutil.rmtree(root.parent / "downloads")
    shutil.rmtree(tmp_path / "source")
    assert bootstrap_module.bootstrap(root) == reference
    assert calls == list(urls.values())
    assert all((sha256(p), p.stat().st_mtime_ns) == record for p, record in before.items())
    moved = tmp_path / "relocated"
    root.rename(moved)
    verify(moved / "reference.json")


def test_bootstrap_resumes_completed_downloads_after_interruption(tmp_path, reference_downloads, monkeypatch):
    _, urls, calls, fetch = reference_downloads
    root = tmp_path / "kegg/snapshot_v1"

    def interrupted_fetch(request, timeout):
        response = fetch(request, timeout)
        if request.full_url == urls["ko_list"]:
            response.headers["Content-Length"] = str(int(response.headers["Content-Length"]) + 1)
        return response

    monkeypatch.setattr(bootstrap_module, "urlopen", interrupted_fetch)
    with pytest.raises(ValueError, match="incomplete download"):
        bootstrap_module.bootstrap(root)
    assert not root.exists()
    cache = root.parent / "downloads"
    assert (cache / "profiles.tar.gz").is_file()
    assert (cache / "profiles.tar.gz.json").is_file()
    assert not (cache / "ko_list.gz").exists()
    before = (cache / "profiles.tar.gz").stat().st_mtime_ns
    monkeypatch.setattr(bootstrap_module, "urlopen", fetch)
    verify(bootstrap_module.bootstrap(root))
    assert calls.count(urls["profiles"]) == 1
    assert calls.count(urls["ko_list"]) == 2
    assert (cache / "profiles.tar.gz").stat().st_mtime_ns == before


def test_bootstrap_rejects_corrupt_snapshot_without_replacing_it(tmp_path, reference_downloads):
    _, _, calls, _ = reference_downloads
    root = tmp_path / "kegg/snapshot_v1"
    bootstrap_module.bootstrap(root)
    count = len(calls)
    profile = root / "profiles/K00001.hmm"
    profile.write_text("damaged")
    with pytest.raises(ValueError, match="mismatch"):
        bootstrap_module.bootstrap(root)
    assert len(calls) == count
    assert profile.read_text() == "damaged"


@pytest.mark.parametrize("corruption", ["gzip", "traversal", "symlink"])
def test_bootstrap_rejects_invalid_archive_without_publishing(tmp_path, reference_downloads, corruption):
    payloads, _, _, _ = reference_downloads
    if corruption == "gzip":
        payloads["profiles"] = payloads["profiles"][:-8]  # missing gzip footer
    else:
        archive = io.BytesIO()
        with tarfile.open(fileobj=archive, mode="w:gz") as handle:
            member = tarfile.TarInfo("../escaped.hmm" if corruption == "traversal" else "profiles/link")
            if corruption == "symlink":
                member.type, member.linkname = tarfile.SYMTYPE, "/tmp/outside"
            handle.addfile(member)
        payloads["profiles"] = archive.getvalue()
    root = tmp_path / "kegg/snapshot_v1"
    with pytest.raises((EOFError, tarfile.TarError, ValueError)):
        bootstrap_module.bootstrap(root)
    assert not root.exists()
    assert not list(root.parent.glob(".snapshot_v1.extracting-*"))
    assert not list(tmp_path.rglob("escaped.hmm"))
