"""Exercise incremental datasets with synthetic RNA-seq artifacts and scheduler doubles."""
import copy
import gzip
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from common import read_tsv, write_json, write_tsv
from dataset import (gg_environment, load, materialize, plan, prepare, status, submit, worker)
from dataset_assets import (COUNTS, identities, import_existing, register_busco, register_quant,
                            register_reference, resolve)

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def dataset_project(tmp_path, tiny_inputs):
    root = tmp_path / "project"
    root.mkdir()
    shutil.copytree(ROOT / "config", root / "config")
    shutil.copytree(ROOT / "profiles", root / "profiles")
    shutil.copytree(ROOT / "workflow", root / "workflow", ignore=shutil.ignore_patterns("__pycache__"))
    shutil.copy2(ROOT / "run_pipeline.sh", root / "run_pipeline.sh")
    shutil.copytree(Path(tiny_inputs["metadata"]).parent, root / "input")
    metadata = root / "input/metadata.tsv"
    rows = [r for r in read_tsv(metadata) if r["run"] != "A2"]
    write_tsv(metadata, list(rows[0]), rows)
    full = root / "input/busco/full"
    full.mkdir(parents=True, exist_ok=True)
    for row in read_tsv(root / "input/busco/summary.tsv"):
        species = row["Species"].replace(" ", "_")
        statuses = []
        for key, state in zip(COUNTS[:-1], ("Complete", "Duplicated", "Fragmented", "Missing")):
            statuses.extend([state] * int(row[key]))
        (full / f"{species}.busco.full.tsv").write_text("# The lineage dataset is: embryophyta_odb12\n" + "".join(
            f"M{i}\t{state}\t{species}_g1\t100\t3\n" for i,state in enumerate(statuses)))
    return root


def imported(root):
    store = root / "resources/dataset_assets"
    import_existing(store, root / "input", root / "input/metadata.tsv")
    return store


@pytest.mark.parametrize("override", [None, "pilot_20260925"])
def test_prepare_cli_uses_config_name_with_optional_override(dataset_project, override):
    root = dataset_project
    imported(root)
    cfg = yaml.safe_load((root / "config/build.yaml").read_text())
    cfg["name"] = "angiosperm_leaf_20260925"
    config = root / "config/named.yaml"
    config.write_text(yaml.safe_dump(cfg))
    command = [sys.executable, str(root / "workflow/scripts/dataset.py"),
               "prepare", "--config", "config/named.yaml"]
    if override is not None:
        command.extend(["--name", override])
    result = subprocess.run(command, cwd=root, text=True, capture_output=True)
    assert result.returncode == 0, result.stderr
    name = override or cfg["name"]
    path = root / "builds" / name
    assert result.stdout.strip() == str(path)
    frozen = load(path, check_code=True)
    assert frozen["name"] == frozen["config"]["name"] == name
    assert frozen["analysis"]["run_name"] == "build_" + name
    assert {row["assembly"] for row in status(path)} == {"reuse"}
    # Renaming in the source config cannot rename an existing preparation.
    cfg["name"] = "angiosperm_leaf_20260926"
    config.write_text(yaml.safe_dump(cfg))
    assert load(path, check_code=True)["name"] == name
    with pytest.raises(ValueError, match="name already exists"):
        prepare(root, name, config)


def test_prepare_legacy_config_requires_explicit_name(dataset_project):
    root = dataset_project
    imported(root)
    config = root / "config/build.yaml"
    cfg = yaml.safe_load(config.read_text())
    cfg.pop("name", None)
    config.write_text(yaml.safe_dump(cfg))
    with pytest.raises(ValueError, match="set name in build config or pass --name"):
        prepare(root, None, config)
    path = prepare(root, "explicit_20260925", config)
    assert load(path)["name"] == "explicit_20260925"


@pytest.mark.parametrize("name", ["", "../escape", "leaf/20260925", "two words", 20260925, True])
def test_prepare_rejects_invalid_config_name(dataset_project, name):
    root = dataset_project
    config = root / "config/build.yaml"
    cfg = yaml.safe_load(config.read_text())
    cfg["name"] = name
    config.write_text(yaml.safe_dump(cfg))
    with pytest.raises(ValueError, match="build name must be a simple directory name"):
        prepare(root, None, config)
    assert not (root / "builds").exists()


def test_manual_metadata_unique_species_and_normalized_collisions(tmp_path):
    path = tmp_path / "metadata.tsv"
    fields = ["scientific_name", "run", "taxid"]
    for names in (("Alpha plant", "Alpha plant"), ("Beta sp-X", "Beta sp_X")):
        write_tsv(path, fields, [dict(zip(fields, [name, f"R{i}", "42"])) for i, name in enumerate(names)])
        with pytest.raises(ValueError, match="one row/run per species"):
            identities(path)


def test_legacy_import_and_changed_metadata_only_reuses_products(dataset_project):
    root = dataset_project
    imported(root)
    cfg = root / "config/build.yaml"
    report = plan(root, cfg)[-1]
    assert [(r["assembly"], r["busco"], r["quant"]) for r in report] == [
        ("reuse", "reuse", "reuse"), ("reuse", "reuse", "reuse"), ("reuse", "reuse", "reuse")]
    path = prepare(root, "base", cfg)
    assert submit(path, until="quant", dry_run=True) == []
    assert not (path / "genegalleon").exists()
    assert load(path)["analysis"]["odb"]["incremental"] is True
    assert load(path)["software_lock"] is None  # Reuse-only preparation does not acquire software.


def test_removal_readdition_and_frozen_membership(dataset_project):
    root = dataset_project
    store = imported(root)
    cfg = root / "config/build.yaml"
    first = prepare(root, "base", cfg)
    first_input = materialize(first)
    assert {r["scientific_name"] for r in read_tsv(first_input / "metadata.tsv")} == {"Alpha plant", "Beta sp-X", "Gamma plant"}
    metadata = root / "input/metadata.tsv"
    original = read_tsv(metadata)
    write_tsv(metadata, list(original[0]), [original[0]])
    # Existing batches are unaffected by later edits to the source metadata.
    assert len(load(first)["items"]) == 3
    second = prepare(root, "removed", cfg)
    second_input = materialize(second)
    assert [r["scientific_name"] for r in read_tsv(second_input / "metadata.tsv")] == ["Alpha plant"]
    assert not (second_input / "cds/Beta_sp-X_longestCDS.fa.gz").exists()
    assert [p.name for p in (second_input / "quant").iterdir()] == ["Alpha_plant"]
    assert (store / "Beta_sp-X").exists()
    write_tsv(metadata, list(original[0]), original)
    third = prepare(root, "restored", cfg)
    assert submit(third, until="quant", dry_run=True) == []
    assert len(read_tsv(materialize(third) / "metadata.tsv")) == 3
    assert len(read_tsv(first_input / "metadata.tsv")) == 3


def test_changed_run_only_requires_quant_and_missing_is_not_silently_dropped(dataset_project):
    root = dataset_project
    imported(root)
    metadata = root / "input/metadata.tsv"
    rows = read_tsv(metadata)
    rows[0]["run"] = "Anew"
    write_tsv(metadata, list(rows[0]), rows)
    report = plan(root, root / "config/build.yaml")[-1]
    assert [report[0][s] for s in ("assembly", "busco", "quant")] == ["reuse", "reuse", "pending"]
    fake_genegalleon(root)
    path = prepare(root, "newrun", root / "config/build.yaml")
    with pytest.raises(ValueError, match="dataset incomplete: Alpha_plant: quant"):
        materialize(path)
    assert not (path / "input").exists()


def test_modified_registered_cds_is_a_conflict(dataset_project):
    root = dataset_project
    imported(root)
    cds = root / "input/cds/Alpha_plant_longestCDS.fa.gz"
    with gzip.open(cds, "wt") as handle: handle.write(">Alpha_plant_g1\nATGCCC\n")
    report = plan(root, root / "config/build.yaml")[-1]
    assert report[0]["assembly"] == "conflict"
    assert "registered file changed" in report[0]["reason"]


def test_frozen_metadata_and_configuration_are_verified(dataset_project):
    root = dataset_project
    imported(root)
    path = prepare(root, "immutable", root / "config/build.yaml")
    (path / "metadata.tsv").write_text("modified")
    with pytest.raises(ValueError, match="registered file changed"):
        status(path)


def fake_genegalleon(root):
    repo = root / "fake_gg"
    (repo / "workflow").mkdir(parents=True)
    (repo / "genegalleon.sif").write_text("test container identifier")
    implementation = r'''
import csv, gzip, json, os, sys
from pathlib import Path
work = Path(os.environ["gg_workspace_dir"])
files = sorted((work / "input/amalgkit_metadata").glob("*.tsv"))
metadata = files[int(os.environ["GG_ARRAY_TASK_ID"]) - 1]
row = next(csv.DictReader(metadata.open(), delimiter="\t"))
species = row["scientific_name"].replace(" ", "_")
run = row["run"]
prefix = "GG_TRANSCRIPTOME_"
out = work / "output/transcriptome_assembly"
out.mkdir(parents=True, exist_ok=True)
with (work / "events.jsonl").open("a") as handle:
    handle.write(json.dumps({"species":species, "env":{k:v for k,v in os.environ.items() if k.startswith(prefix)}})+"\n")
assert os.environ["LC_ALL"] == os.environ["SINGULARITYENV_LC_ALL"] == os.environ["APPTAINERENV_LC_ALL"] == "C"
assert os.environ[prefix+"RUN_MULTISPECIES_SUMMARY"] == "0"
assert os.environ[prefix+"REMOVE_AMALGKIT_FASTQ_AFTER_COMPLETION"] == "0"
if os.environ[prefix+"RUN_AMALGKIT_GETFASTQ"] == "1" and os.environ.get("FAKE_GG_FAIL_DOWNLOAD"):
    sys.exit(29)
if os.environ[prefix+"RUN_ASSEMBLY"] == "1":
    cds = out / "longest_cds" / (species+"_longestCDS.fa.gz")
    cds.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(cds,"wt") as h:
        h.write(">"+species+"_g1\nATGAAATAA\n>"+species+"_g2\nATGCCCTAA\n")
    if os.environ.get("FAKE_GG_FAIL_ASSEMBLY"): sys.exit(23)
if os.environ[prefix+"RUN_BUSCO_LONGEST_CDS"] == "1":
    for directory, suffix, content in [("busco_full_longest_cds","full.tsv", "# The lineage dataset is: embryophyta_odb12 (test)\nB1\tComplete\t"+species+"_g1:0-9\t100\t3\nB2\tComplete\t"+species+"_g2:0-9\t100\t3\n"), ("busco_short_longest_cds","short.txt","C:100%[S:100%,D:0%],F:0%,M:0%,n:2\n")]:
        path = out / directory / (species+"_busco."+suffix)
        path.parent.mkdir(parents=True, exist_ok=True); path.write_text(content)
if os.environ[prefix+"RUN_AMALGKIT_QUANT"] == "1":
    abundance = out / "amalgkit_quant" / species / run / (run+"_abundance.tsv")
    abundance.parent.mkdir(parents=True, exist_ok=True)
    abundance.write_text("target_id\ttpm\n"+species+"_g1\t500000\n"+species+"_g2\t500000\n")
    if not os.environ.get("FAKE_GG_INCOMPLETE_MERGE"):
        for suffix in ("eff_length","est_counts","tpm","metadata"):
            path = out / "amalgkit_merge" / species / (species+"_"+suffix+".tsv")
            path.parent.mkdir(parents=True, exist_ok=True); path.write_text("synthetic\n")
'''
    script = repo / "workflow/gg_transcriptome_generation_entrypoint.sh"
    script.write_text("#!/usr/bin/env bash\nexec " + sys.executable + " - <<'PY'\n" + implementation + "\nPY\n")
    cfg = yaml.safe_load((root / "config/build.yaml").read_text())
    cfg["genegalleon"]["repository"] = str(repo)
    (root / "config/build.yaml").write_text(yaml.safe_dump(cfg))
    return repo


def new_dataset(root, names=("New plant",), array_size=None):
    fake_genegalleon(root)
    if array_size is not None:
        cfg = yaml.safe_load((root / "config/build.yaml").read_text())
        cfg["slurm"]["array_size"] = array_size
        (root / "config/build.yaml").write_text(yaml.safe_dump(cfg))
    fields = ["scientific_name", "run", "taxid"]
    write_tsv(root / "input/new.tsv", fields, [dict(zip(fields, [name, f"SRR{i+1}", "42"])) for i,name in enumerate(names)])
    return prepare(root, "addition", root / "config/build.yaml", "input/new.tsv")


@pytest.mark.parametrize("key", ["amalgkit_rrna_filter", "amalgkit_contam_filter"])
@pytest.mark.parametrize("value", ["yes", "no"])
def test_read_filter_settings_are_not_configurable(dataset_project, key, value):
    config = dataset_project / "config/build.yaml"
    cfg = yaml.safe_load(config.read_text())
    cfg["genegalleon"]["settings"][key] = value
    config.write_text(yaml.safe_dump(cfg))
    with pytest.raises(ValueError, match=f"managed/invalid GeneGalleon setting: {key}"):
        plan(dataset_project, config)


def test_staged_workers_reuse_and_native_array_filename_order(dataset_project, monkeypatch):
    root = dataset_project
    for key in ("AMALGKIT_RRNA_FILTER", "AMALGKIT_CONTAM_FILTER"):
        monkeypatch.setenv("GG_TRANSCRIPTOME_" + key, "yes")
    # Prefix species sort order differs from the order of *_metadata.tsv filenames.
    path = new_dataset(root, ("New plant", "New plant alba"))
    commands = submit(path, until="busco", dry_run=True)
    assert len(commands) == 2
    assert "--array=1,2%5" in commands[0]
    assert "--dependency=afterok:JOB_ID_assembly" in commands[1]
    for stage in ("assembly", "busco", "quant"):
        for index in (1, 2): worker(path, stage, index)
    events = [json.loads(line) for line in (path / "genegalleon/events.jsonl").read_text().splitlines()]
    assert [e["species"] for e in events] == ["New_plant", "New_plant_alba"] * 3
    for event in events:
        assert event["env"]["GG_TRANSCRIPTOME_AMALGKIT_RRNA_FILTER"] == "no"
        assert event["env"]["GG_TRANSCRIPTOME_AMALGKIT_CONTAM_FILTER"] == "no"
    assert all(r["quant"] == "reuse" for r in status(path))
    worker(path, "assembly", 1)
    assert len((path / "genegalleon/events.jsonl").read_text().splitlines()) == 6
    assert submit(path, until="quant", dry_run=True) == []
    assert len(read_tsv(materialize(path) / "metadata.tsv")) == 2


def test_failed_stage_is_not_registered_and_retry_is_limited(dataset_project, monkeypatch):
    path = new_dataset(dataset_project)
    submit(path, until="quant", dry_run=True)
    monkeypatch.setenv("FAKE_GG_FAIL_ASSEMBLY", "1")
    with pytest.raises(subprocess.CalledProcessError): worker(path, "assembly", 1)
    assert status(path)[0]["assembly"] == "pending"
    monkeypatch.delenv("FAKE_GG_FAIL_ASSEMBLY")
    worker(path, "assembly", 1)
    assert status(path)[0]["assembly"] == "reuse"
    assert list((path / "jobs/incomplete/New_plant/assembly").rglob("*.gz"))
    worker(path, "busco", 1)
    monkeypatch.setenv("FAKE_GG_INCOMPLETE_MERGE", "1")
    with pytest.raises(ValueError, match="quant/merge did not finish"): worker(path, "quant", 1)
    assert status(path)[0]["quant"] == "pending"
    monkeypatch.delenv("FAKE_GG_INCOMPLETE_MERGE")
    worker(path, "quant", 1)
    assert status(path)[0]["quant"] == "reuse"


def test_completed_new_species_reused_by_next_dataset(dataset_project):
    root = dataset_project
    path = new_dataset(root)
    submit(path, until="quant", dry_run=True)
    for stage in ("assembly", "busco", "quant"): worker(path, stage, 1)
    second = prepare(root, "next", root / "config/build.yaml", "input/new.tsv")
    assert submit(second, until="quant", dry_run=True) == []
    assert read_tsv(materialize(second) / "metadata.tsv")[0]["run"] == "SRR1"


def test_partial_pilot_never_exports_incomplete_dataset(dataset_project):
    path = new_dataset(dataset_project, ("New plant", "Other plant"))
    subset = dataset_project / "pilot.txt"
    subset.write_text("Other_plant\n")
    commands = submit(path, until="mapping", species=subset, dry_run=True)
    assert len(commands) == 3
    assert all("--array=2%5" in cmd for cmd in commands)
    assert not any("mapping" in cmd[-1] for cmd in commands)
    for stage in ("assembly", "busco", "quant"): worker(path, stage, 2)
    with pytest.raises(ValueError, match="dataset incomplete: New_plant"):
        materialize(path)
    assert "--array=1%5" in submit(path, until="assembly", dry_run=True)[0]


def test_slurm_submission_dependencies_and_duplicate_submission_guard(dataset_project, monkeypatch):
    path = new_dataset(dataset_project)
    calls = []
    active = False
    def scheduler(command, **kwargs):
        nonlocal active
        calls.append(command)
        if command[0] == "squeue": return "1001\n" if active else ""
        return str(1000 + len([c for c in calls if c[0] == "sbatch"])) + ";cluster\n"
    monkeypatch.setattr(subprocess, "check_output", scheduler)
    submit(path, until="mapping")
    batch = json.loads((path / "jobs/submission_0001.json").read_text())
    assert len(batch["jobs"]) == 4
    assert "--dependency=afterok:1001" in batch["jobs"][1]["command"]
    assert "--dependency=afterok:1003" in batch["jobs"][3]["command"]
    active = True
    with pytest.raises(ValueError, match="queued/running"):
        submit(path, until="mapping")
    assert len([c for c in calls if c[0] == "sbatch"]) == 4


def test_private_relative_reads_are_frozen_and_reuse_detects_changed_bytes(dataset_project):
    root = dataset_project
    fake_genegalleon(root)
    reads = root / "input/local.fastq"
    reads.write_text("@r1\nATGC\n+\nIIII\n")
    fields = ["scientific_name", "run", "taxid", "private_file", "lib_layout", "read1_path"]
    write_tsv(root / "input/private.tsv", fields,
              [dict(zip(fields, ["Private plant", "LOCAL1", "42", "yes", "single", "local.fastq"]))])
    cfg = root / "config/build.yaml"
    path = prepare(root, "private", cfg, "input/private.tsv")
    submit(path, until="quant", dry_run=True)
    staged = path / "genegalleon/input/reads/Private_plant/read1_path.fastq"
    assert staged.read_bytes() == reads.read_bytes()
    assert read_tsv(path / "genegalleon/input/amalgkit_metadata/Private_plant_metadata.tsv")[0]["read1_path"] == "/workspace/input/reads/Private_plant/read1_path.fastq"
    original = reads.read_bytes()
    reads.write_text("changed after submission")
    with pytest.raises(ValueError, match="registered file changed"):
        worker(path, "assembly", 1)
    reads.write_bytes(original)
    for stage in ("assembly", "busco", "quant"): worker(path, stage, 1)
    assert plan(root, cfg, "input/private.tsv")[-1][0]["quant"] == "reuse"
    second = prepare(root, "private_reused", cfg, "input/private.tsv")
    assert submit(second, until="quant", dry_run=True) == []
    reads.write_text("different reads under same run")
    report = plan(root, cfg, "input/private.tsv")[-1][0]
    assert report["quant"] == "conflict"
    assert "registered file changed" in report["reason"]
    reads.unlink()
    assert plan(root, cfg, "input/private.tsv")[-1][0]["quant"] == "reuse"


def test_retry_quarantine_does_not_touch_another_species_with_same_prefix(dataset_project, monkeypatch):
    path = new_dataset(dataset_project, ("New plant", "New plant alba"))
    submit(path, until="assembly", dry_run=True)
    monkeypatch.setenv("FAKE_GG_FAIL_ASSEMBLY", "1")
    with pytest.raises(subprocess.CalledProcessError): worker(path, "assembly", 1)
    monkeypatch.delenv("FAKE_GG_FAIL_ASSEMBLY")
    worker(path, "assembly", 2)
    other = path / "genegalleon/output/transcriptome_assembly/assembled_transcripts_with_isoforms/New_plant_alba_isoform.fa.gz"
    other.parent.mkdir(parents=True, exist_ok=True)
    other.write_bytes(b"completed output of another species")
    worker(path, "assembly", 1)
    assert other.read_bytes() == b"completed output of another species"
    assert all(r["assembly"] == "reuse" for r in status(path))


def test_split_slurm_arrays_preserve_species_identity_and_bound_concurrency(dataset_project, monkeypatch):
    path = new_dataset(dataset_project, ("Alpha new", "Beta new", "Gamma new"), array_size=2)
    commands = submit(path, until="assembly", dry_run=True)
    assert len(commands) == 2
    assert "--array=1,2%5" in commands[0]
    assert "--array=1%5" in commands[1]
    assert "--dependency=afterok:JOB_ID_assembly" in commands[1]
    # Execute the generated high-index batch locally with the GeneGalleon double.
    # Local Slurm index 1 must select logical species 3, not species 1.
    subprocess.run(["bash", commands[1][-1]], cwd=dataset_project,
                   env=dict(os.environ, SLURM_ARRAY_TASK_ID="1"), check=True)
    assert [r["assembly"] for r in status(path)] == ["pending", "pending", "reuse"]
    calls = []
    def scheduler(command, **kwargs):
        calls.append(command)
        assert command[0] == "sbatch"
        return str(1000 + len(calls)) + "\n"
    monkeypatch.setattr(subprocess, "check_output", scheduler)
    submit(path, until="busco")
    batch = json.loads((path / "jobs/submission_0002.json").read_text())
    assert [j["array_offset"] for j in batch["jobs"]] == [0, 0, 2]
    assert [j["indices"] for j in batch["jobs"]] == [[1, 2], [1, 2], [3]]
    assert "--dependency=afterok:1001" in calls[1]
    assert "--dependency=afterok:1002" in calls[2]


def test_extra_native_metadata_file_cannot_shift_species_array_index(dataset_project):
    path = new_dataset(dataset_project)
    submit(path, until="assembly", dry_run=True)
    (path / "genegalleon/input/amalgkit_metadata/Aardvark_backup.tsv").write_text("unexpected metadata")
    with pytest.raises(ValueError, match="metadata file set changed"):
        worker(path, "assembly", 1)
    assert not (path / "genegalleon/events.jsonl").exists()


def test_prepare_resolves_automatic_dependencies_once_and_freezes_the_lock(dataset_project, monkeypatch):
    import dataset_software
    root = dataset_project
    repository = fake_genegalleon(root)
    config_path = root / "config/build.yaml"
    cfg = yaml.safe_load(config_path.read_text())
    cfg["genegalleon"]["repository"] = None
    cfg["genegalleon"]["image"] = None
    config_path.write_text(yaml.safe_dump(cfg))
    fields = ["scientific_name", "run", "taxid"]
    write_tsv(root / "input/automatic.tsv", fields, [dict(zip(fields, ["New plant", "SRR1", "42"]))])
    calls = []
    def source(config):
        calls.append("source")
        return repository, {"kind": "downloaded_source", "identity": {"revision": config["revision"]},
                            "files": dataset_software.source_records(repository)}
    def image(config):
        calls.append("image")
        path = repository / "genegalleon.sif"
        return path, {"kind": "downloaded_image", "identity": {"uri": config["image_uri"]},
                      "files": [dataset_software.record(path)]}
    monkeypatch.setattr(dataset_software, "fetch_source", source)
    monkeypatch.setattr(dataset_software, "fetch_image", image)
    # Planning must not fetch dependencies even when every upstream stage is missing.
    assert plan(root, config_path, "input/automatic.tsv")[-1][0]["assembly"] == "pending"
    assert calls == []
    path = prepare(root, "automatic", config_path, "input/automatic.tsv")
    manifest = load(path, check_code=True)
    assert calls == ["source", "image"]
    assert manifest["config"]["genegalleon"]["repository"] == str(repository)
    assert manifest["software_lock"]["source"]["identity"]["revision"] == cfg["genegalleon"]["revision"]
    cfg["genegalleon"]["revision"] = "f" * 40
    config_path.write_text(yaml.safe_dump(cfg))
    submit(path, until="assembly", dry_run=True)
    worker(path, "assembly", 1)
    assert calls == ["source", "image"]
    assert load(path)["software_lock"] == manifest["software_lock"]
