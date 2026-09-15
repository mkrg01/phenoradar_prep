"""Real rule planning uses 100-species chunks and CPU-dependent ODB batches."""
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest
import yaml

from common import write_tsv

ROOT = Path(__file__).resolve().parents[1]


def test_hundred_species_chunks_and_cpu_dependent_batches(
    tmp_path, workflow_project, tiny_inputs, seed_taxonomy, fake_odb,
    frozen_reference, command_environment,
):
    snakemake = os.environ.get("SNAKEMAKE_BIN") or shutil.which("snakemake")
    if not snakemake:
        pytest.skip("Snakemake required")
    seed_taxonomy(tiny_inputs["taxonomy_db"])
    source = workflow_project / "input"
    metadata, busco = [], []
    species = [f"Species_{i:03d}" for i in range(101)]
    for i, name in enumerate(species):
        run = f"R{i:03d}"
        with gzip.open(source / "cds" / f"{name}_longestCDS.fa.gz", "wt") as handle:
            handle.write(f">{name}_g1\nATGAAATAA\n")
        write_tsv(source / "quant" / name / run / f"{run}_abundance.tsv", ["target_id", "tpm"],
                  [{"target_id": f"{name}_g1", "tpm": 100}])
        metadata.append({"scientific_name": name.replace("_", " "), "run": run, "taxid": 42})
        busco.append({"Species": name.replace("_", " "), "busco_cds_single": 1,
                      "busco_cds_duplicated": 0, "busco_cds_fragmented": 0,
                      "busco_cds_missing": 0, "busco_cds_total": 1})
    write_tsv(source / "metadata.tsv", list(metadata[0]), list(reversed(metadata)))
    write_tsv(source / "busco/summary.tsv", list(busco[0]), busco)
    override = tmp_path / "override.yaml"
    override.write_text(yaml.safe_dump({"run_name": "test"}))
    reference = workflow_project / "resources/orthodb/v12_3193"
    reference.parent.mkdir(parents=True)
    reference.symlink_to(frozen_reference, target_is_directory=True)
    events = tmp_path / "map_events"
    env = {**command_environment({"python": sys.executable, "ODB-mapper": fake_odb}),
           "FAKE_ODB_LOG": str(events), "FAKE_ODB_BARRIER_COUNT": "2"}
    base = [snakemake, "--snakefile", str(ROOT / "workflow/Snakefile"), "--configfile", str(override),
            "--cores", "2", "--resources", "mem_mb=16000", "--printshellcmds",
            "--set-resources", "odb_map:mem_mb=3000"]

    def execute(targets, options=(), threads=1):
        result = subprocess.run(base + ["--set-threads", f"odb_map={threads}", *options, "--", *targets],
                                cwd=workflow_project, env=env, text=True, capture_output=True, timeout=120)
        assert result.returncode == 0, result.stdout + result.stderr + "\n" + "\n".join(
            p.read_text()[-3000:] for p in (tmp_path / "logs").rglob("*.log"))
        return result.stdout + result.stderr

    execute(["results/test/orthogroups/mapping/manifests"])
    mapping = tmp_path / "results/test/orthogroups/mapping"
    chunks = json.loads((mapping / "manifests/chunks.json").read_text())
    assert [len(row["species"]) for row in chunks] == [100, 1]
    assert [name for row in chunks for name in row["species"]] == species

    # Supply completed translations so the test exercises scheduling and mapping only.
    for row in chunks:
        manifest = mapping / "manifests" / f'{row["chunk"]}.fs'
        assert manifest.read_text().splitlines() == row["proteins"]
        for name, filename in zip(row["species"], row["proteins"]):
            path = Path(filename)
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(f">{name}_g1\nMK\n")
            path.with_suffix(".json").write_text("{}\n")
    targets = [str((mapping / "chunks" / row["chunk"] / f'{row["chunk"]}.og.annotations')
                   .relative_to(workflow_project)) for row in chunks]
    # An 8-thread override is capped by the two available CPUs; batch size follows the cap.
    dry = execute(targets, options=("--dry-run",), threads=8)
    assert dry.count("rule odb_map:") == 2
    assert dry.count("--jobs 2 --batch-size 8") == 2
    execute(targets)
    assert len(events.read_text().splitlines()) == 2  # Both MAP jobs overlapped at the barrier.
    for row in chunks:
        out = mapping / "chunks" / row["chunk"]
        identity = json.loads((out / "provenance.json").read_text())["identity"]
        assert identity["jobs"] == 1 and identity["batch_size"] == 4
        configuration = (out / "orthologer_conf.sh").read_text()
        assert "OP_NJOBMAX_BATCH=4\n" in configuration
        assert "OP_NJOBMAX_LOCAL=1\n" in configuration
    before = {path: path.stat().st_mtime_ns for path in mapping.rglob("provenance.json")}
    assert "Nothing to be done" in execute(targets)
    assert len(events.read_text().splitlines()) == 2
    assert all(path.stat().st_mtime_ns == stamp for path, stamp in before.items())
