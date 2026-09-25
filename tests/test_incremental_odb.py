"""Mixed ODB imports map only missing species and survive dataset removal/readdition."""
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from common import file_record, read_tsv, write_json, write_tsv
from incremental_odb import plan
from merge_odb import merge
from translate_cds import translate

ROOT = Path(__file__).resolve().parents[1]


def snapshot(root, proteins, species):
    root.mkdir(parents=True)
    annotation = root / "annotations.tsv"
    annotation.write_text("#query\tODB_OG\n" + "".join(f"{s}_g1\tOG1\n{s}_g2\tOG2\n" for s in species))
    write_json(root / "snapshot.json", {"schema_version": 1, "version": "v12", "node": 3193,
        "proteins": [{"species": s, "odb_species": s.replace("-", "_"),
                      **file_record(proteins / (s.replace("-", "_") + "_protein.fa"))} for s in species],
        "annotations": dict(file_record(annotation), path="annotations.tsv")})


def test_plan_and_merge_do_not_leak_other_species_from_import(tmp_path):
    proteins = tmp_path / "proteins"
    proteins.mkdir()
    names = ["Alpha_plant", "Beta_plant", "Removed_plant"]
    for s in names: (proteins / f"{s}_protein.fa").write_text(f">{s}_g1\nMK\n>{s}_g2\nMP\n")
    old = tmp_path / "old"
    snapshot(old, proteins, names)
    samples = tmp_path / "samples.tsv"
    rows = [{"species": s, "odb_species": s} for s in names[:2]]
    write_tsv(samples, list(rows[0]), rows)
    out = tmp_path / "plan"
    result = plan(samples, proteins, out, tmp_path / "cache", existing=old)
    assert result["mapped_species"] == []
    merge(samples, out / "chunks.json", tmp_path / "chunks", proteins, tmp_path / "db.sqlite",
          tmp_path / "mapping.tsv", tmp_path / "qc.json", source_plan=out / "plan.json")
    assert not any(r["#query"].startswith("Removed_") for r in read_tsv(tmp_path / "mapping.tsv"))
    assert json.loads((tmp_path / "qc.json").read_text())["excluded_annotation_rows"] == 2


def test_mixed_mapping_and_cross_run_reuse(tiny_inputs, fake_odb, frozen_reference,
                                          command_environment, workflow_project, seed_taxonomy):
    snakemake = os.environ.get("SNAKEMAKE_BIN") or shutil.which("snakemake")
    seqkit = os.environ.get("SEQKIT_BIN") or shutil.which("seqkit")
    if not snakemake or not seqkit: pytest.skip("Snakemake and seqkit required")
    root = workflow_project
    seed_taxonomy(tiny_inputs["taxonomy_db"])
    original_proteins = root / "original_proteins"
    original_proteins.mkdir()
    translate(str(Path(tiny_inputs["cds_dir"]) / "Alpha_plant_longestCDS.fa.gz"),
              original_proteins / "Alpha_plant_protein.fa", original_proteins / "Alpha_plant.json", seqkit=seqkit)
    old = root / "existing"
    snapshot(old, original_proteins, ["Alpha_plant"])
    reference = root / "resources/orthodb/v12_3193"
    reference.parent.mkdir(parents=True)
    reference.symlink_to(frozen_reference, target_is_directory=True)
    events = root / "events.txt"
    env = command_environment({"python": sys.executable, "seqkit": seqkit, "ODB-mapper": fake_odb})
    env["FAKE_ODB_LOG"] = str(events)
    config = {"run_name": "first", "odb": {"incremental": True, "existing_results": str(old), "chunk_size": 1}}
    override = root / "override.yaml"
    def execute():
        override.write_text(yaml.safe_dump(config))
        cmd = [snakemake, "--snakefile", str(ROOT / "workflow/Snakefile"), "--configfile", str(override),
               "--cores", "2", "--resources", "mem_mb=16000", "--set-threads", "odb_map=1",
               "--set-resources", "odb_map:mem_mb=3000"]
        result = subprocess.run(cmd, cwd=root, env=env, capture_output=True, text=True, timeout=120)
        logs = "\n".join(str(p) + ": " + p.read_text()[-3000:] for p in (root / "logs").rglob("*.log"))
        assert result.returncode == 0, result.stdout + result.stderr + logs
        return root / "results" / config["run_name"]
    first = execute()
    assert len(events.read_text().splitlines()) == 1
    assert json.loads((first / "orthogroups/mapping/merge_qc.json").read_text())["mode"] == "mixed"
    assert len(list((root / "resources/odb_cache/v12_3193").glob("*/snapshot.json"))) == 1
    config["run_name"] = "same_species_new_run"
    second = execute()
    assert len(events.read_text().splitlines()) == 1
    assert json.loads((second / "orthogroups/mapping/merge_qc.json").read_text())["mode"] == "existing"
    metadata = root / "input/metadata.tsv"
    all_rows = read_tsv(metadata)
    write_tsv(metadata, list(all_rows[0]), [r for r in all_rows if r["scientific_name"] == "Beta sp-X"])
    config["run_name"] = "removed_alpha"
    removed = execute()
    assert {r["species"] for r in read_tsv(removed / "orthogroups/expression/tpm.tsv")} == {"Beta_sp-X"}
    assert not (removed / "proteins/Alpha_plant_protein.fa").exists()
    assert len(events.read_text().splitlines()) == 1
    write_tsv(metadata, list(all_rows[0]), all_rows)
    config["run_name"] = "restored"
    execute()
    assert len(events.read_text().splitlines()) == 1
