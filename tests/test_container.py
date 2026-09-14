"""Container generation, environment coverage, and verified ASTRAL export."""
import json
import os
from pathlib import Path
import re
import shutil

import pytest

from common import file_record, write_json
from generate_container import generate
from prepare_phylogeny_tools import SOURCES, prepare

ROOT = Path(__file__).resolve().parents[1]


def test_generation_covers_optional_environments_and_hashes_deploy_scripts(tmp_path):
    snakemake = os.environ.get("SNAKEMAKE_BIN") or shutil.which("snakemake")
    if not snakemake:
        pytest.skip("Snakemake 9.8.0 required for container generation")
    shutil.copytree(ROOT / "workflow", tmp_path / "workflow",
                    ignore=shutil.ignore_patterns("__pycache__"))
    first = generate(tmp_path, snakemake)
    names = set(re.findall(r"COPY workflow/envs/([\w.-]+)\.yaml", first))
    assert names == {p.stem for p in (ROOT / "workflow/envs").glob("*.yaml")}
    assert not (tmp_path / ".snakemake").exists()
    assert "@sha256:" in first and "miniforge3:latest" not in first
    for name in ("dating", "monophy"):
        assert f"COPY workflow/envs/{name}.post-deploy.sh" in first
    before = re.search(r"COPY workflow/envs/dating.yaml (\S+)", first).group(1)
    script = tmp_path / "workflow/envs/dating.post-deploy.sh"
    script.write_text(script.read_text() + "\n# Test a changed deployment recipe.\n")
    second = generate(tmp_path, snakemake)
    after = re.search(r"COPY workflow/envs/dating.yaml (\S+)", second).group(1)
    assert before != after
    assert generate(tmp_path, snakemake) == second


def test_container_generation_refuses_unsupported_pin_files(tmp_path, monkeypatch):
    shutil.copytree(ROOT / "workflow/envs", tmp_path / "workflow/envs")
    (tmp_path / "workflow/envs/dating.linux-64.pin.txt").write_text("@EXPLICIT\n")
    monkeypatch.setattr("generate_container.subprocess.check_output", lambda *a, **k: "9.8.0\n")
    with pytest.raises(ValueError, match="does not apply .pin.txt"):
        generate(tmp_path)


@pytest.mark.parametrize("corrupt", [False, True])
def test_astral_bundle_validates_source_and_binary_without_a_compiler(tmp_path, monkeypatch, corrupt):
    bundle = tmp_path / "bundle"
    (bundle / "bin").mkdir(parents=True)
    binary = bundle / "bin/astral4_int128"
    binary.write_bytes(b"verified test binary")
    write_json(bundle / "aster.json", {"source": SOURCES["aster"],
               "command": ["g++", "-D", "LARGE_DATA"], "executable": file_record(binary)})
    if corrupt:
        binary.write_bytes(b"changed binary")
    monkeypatch.setenv("PHENORADAR_PHYLOGENY_TOOLS", str(bundle))
    monkeypatch.setenv("PATH", "")
    destination = tmp_path / "export"
    if corrupt:
        with pytest.raises(ValueError, match="differs from pinned build"):
            prepare(destination)
        assert not (destination / "bin/astral4_int128").exists()
    else:
        prepare(destination)
        exported = destination / "bin/astral4_int128"
        assert exported.read_bytes() == binary.read_bytes()
        record = json.loads((destination / "aster.json").read_text())
        assert record["executable"] == file_record(exported)
        assert record["bundle_executable"] == file_record(binary)
        before = (destination / "aster.json").read_bytes()
        prepare(destination)
        assert (destination / "aster.json").read_bytes() == before
