"""Taxonomy bootstrap tests use a small real taxdump and never contact NCBI."""
import io
import json
import os
import shutil
import sqlite3
import subprocess
import tarfile
from pathlib import Path

import pytest
import yaml

import prepare_taxonomy
from common import file_record

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def tiny_taxdump(tiny_inputs):
    with sqlite3.connect(tiny_inputs["taxonomy_db"]) as db:
        taxa = db.execute("SELECT taxid, parent, spname, rank FROM species").fetchall()
    files = {
        "nodes.dmp": "".join(f"{taxid}\t|\t{parent}\t|\t{rank}\t|\n"
                             for taxid, parent, name, rank in taxa),
        "names.dmp": "".join(f"{taxid}\t|\t{name}\t|\t\t|\tscientific name\t|\n"
                             for taxid, parent, name, rank in taxa),
        "merged.dmp": "99\t|\t42\t|\n",
    }
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w:gz") as archive:
        for name, text in files.items():
            content = text.encode()
            info = tarfile.TarInfo(name)
            info.size = len(content)
            archive.addfile(info, io.BytesIO(content))
    return buffer.getvalue()


def test_download_builds_queryable_frozen_snapshot(tiny_taxdump, tmp_path, monkeypatch):
    calls = []

    def fake_urlopen(url, timeout):
        calls.append((url, timeout))
        return io.BytesIO(tiny_taxdump)

    monkeypatch.setattr(prepare_taxonomy, "urlopen", fake_urlopen)
    destination = tmp_path / "taxonomy" / "taxa.sqlite"
    prepare_taxonomy.prepare(destination)
    assert calls == [(prepare_taxonomy.TAXDUMP_URL, 60)]
    from ete4 import NCBITaxa
    ncbi = NCBITaxa(dbfile=str(destination), update=False)
    try:
        assert ncbi.get_lineage(42) == [1, 2, 42]
        assert ncbi.get_taxid_translator([42]) == {42: "Alpha plant"}
        assert ncbi.get_lineage(99) == [1, 2, 42]
    finally:
        ncbi.db.close()
    record_path = Path(str(destination) + ".json")
    record = json.loads(record_path.read_text())
    assert record["source"] == prepare_taxonomy.TAXDUMP_URL
    assert record["method"] == "ncbi_download"
    assert record["snapshot"] == file_record(destination)
    assert record["taxdump"]["bytes"] == len(tiny_taxdump)
    before = (destination.read_bytes(), destination.stat().st_mtime_ns, record_path.read_bytes())
    prepare_taxonomy.prepare(destination)
    assert len(calls) == 1
    assert before == (destination.read_bytes(), destination.stat().st_mtime_ns, record_path.read_bytes())
    assert not list(destination.parent.glob(".taxa.sqlite.preparing-*"))


def test_local_source_is_offline_and_existing_database_is_preserved(tiny_inputs, tmp_path, monkeypatch):
    def no_download(*args, **kwargs):
        pytest.fail("a local or existing database must not trigger a download")

    monkeypatch.setattr(prepare_taxonomy, "urlopen", no_download)
    source = Path(tiny_inputs["taxonomy_db"])
    original = file_record(source)
    destination = tmp_path / "taxonomy" / "taxa.sqlite"
    prepare_taxonomy.prepare(destination, source)
    assert file_record(source) == original
    with sqlite3.connect(source) as src, sqlite3.connect(destination) as dst:
        assert list(src.iterdump()) == list(dst.iterdump())
    record_path = Path(str(destination) + ".json")
    record = json.loads(record_path.read_text())
    assert record["source"] == str(source)
    assert record["method"] == "sqlite_backup"
    assert record["snapshot"] == file_record(destination)
    before = (destination.read_bytes(), destination.stat().st_mtime_ns)
    record_path.unlink()  # Externally supplied databases need not have a sidecar.
    source.unlink()
    prepare_taxonomy.prepare(destination, source)
    assert before == (destination.read_bytes(), destination.stat().st_mtime_ns)
    assert not record_path.exists()


@pytest.mark.parametrize("failure", ["download", "invalid_archive"])
def test_failed_bootstrap_leaves_no_database_and_can_be_retried(tiny_inputs, tmp_path, monkeypatch, failure):
    def broken_download(url, timeout):
        if failure == "download":
            raise OSError("download interrupted")
        return io.BytesIO(b"invalid taxdump")

    monkeypatch.setattr(prepare_taxonomy, "urlopen", broken_download)
    destination = tmp_path / "taxonomy" / "taxa.sqlite"
    with pytest.raises((OSError, subprocess.CalledProcessError)):
        prepare_taxonomy.prepare(destination)
    assert not destination.exists()
    assert not Path(str(destination) + ".json").exists()
    assert not list(destination.parent.glob(".taxa.sqlite.preparing-*"))
    prepare_taxonomy.prepare(destination, tiny_inputs["taxonomy_db"])
    assert destination.is_file()


def test_missing_database_is_scheduled_without_downloading(tmp_path):
    snakemake = os.environ.get("SNAKEMAKE_BIN") or shutil.which("snakemake")
    if not snakemake:
        pytest.skip("Snakemake is not available")
    destination = tmp_path / "taxonomy" / "taxa.sqlite"
    configfile = tmp_path / "config.yaml"
    configfile.write_text(yaml.safe_dump({"taxonomy": {"database": str(destination)}}))
    result = subprocess.run([
        snakemake, "--snakefile", str(ROOT / "workflow/Snakefile"),
        "--configfile", str(configfile), "--cores", "1", "--dry-run", "--", str(destination),
    ], cwd=ROOT, env={**os.environ, "XDG_CACHE_HOME": str(tmp_path / "cache")},
        capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "rule prepare_taxonomy:" in result.stdout
    assert not destination.exists()
