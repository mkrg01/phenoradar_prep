import gzip
import json
import os
import shutil
import sqlite3
import subprocess
import sys
from pathlib import Path

import pytest

from aggregate_tpm import aggregate
from common import file_record, read_tsv, write_json, write_tsv
from make_manifests import make
from merge_odb import annotation_pairs, merge as merge_odb
from merge_tpm import merge as merge_tpm
from prepare_metadata import prepare
from run_odb_chunk import run as run_odb
from snapshot_taxonomy import snapshot
from verify_odb_reference import verify

ROOT = Path(__file__).resolve().parents[1]


def test_selection_and_explicit_run_paths(tiny_inputs, tmp_path):
    out = tmp_path / "metadata"
    prepare(**tiny_inputs, outdir=out)
    rows = read_tsv(out / "samples.tsv")
    assert [r["run"] for r in rows] == ["A1", "A2", "B1"]
    assert rows[-1]["species"] == "Beta_sp-X"
    assert rows[-1]["odb_species"] == "Beta_sp_X"
    assert all(Path(r["abundance"]).is_file() for r in rows)
    assert {r["kingdom"] for r in read_tsv(out / "metadata_high_busco.tsv")} == {"Plants"}
    assert json.loads((out / "selection.json").read_text())["selected_species"] == 2
    subset = tmp_path / "subset.txt"
    subset.write_text("Beta_sp-X\n")
    prepare(**tiny_inputs, outdir=out, species_list=subset)
    assert [r["run"] for r in read_tsv(out / "samples.tsv")] == ["B1"]


@pytest.mark.parametrize("threshold", [0, 0.5])
def test_missing_busco_species_are_filtered_and_recorded(tiny_inputs, tmp_path, threshold):
    rows = read_tsv(tiny_inputs["busco"])[1:]  # Alpha plant has two runs.
    write_tsv(tiny_inputs["busco"], list(rows[0]), rows)
    # Excluded species must not require CDS, abundance, or taxonomy records.
    (Path(tiny_inputs["cds_dir"]) / "Alpha_plant_longestCDS.fa.gz").unlink()
    shutil.rmtree(Path(tiny_inputs["quant_dir"]) / "Alpha_plant")
    with sqlite3.connect(tiny_inputs["taxonomy_db"]) as db:
        db.execute("DELETE FROM species WHERE taxid = 42")
    out = tmp_path / "metadata"
    prepare(**tiny_inputs, outdir=out, threshold=threshold)
    assert [r["run"] for r in read_tsv(out / "samples.tsv")] == (["B1", "G1"] if threshold == 0 else ["B1"])
    all_rows = read_tsv(out / "metadata_all.tsv")
    assert len(all_rows) == 4
    excluded = [r for r in all_rows if r["species"] == "Alpha_plant"]
    assert len(excluded) == 2
    assert all(r["selected"] == "False" and r["busco_percent"] == "" for r in excluded)
    report = json.loads((out / "selection.json").read_text())
    assert report["input_species"] == 3
    assert report["missing_busco_species"] == ["Alpha plant"]
    assert report["missing_busco_runs"] == 2
    assert report["unknown_taxids"] == []
    # Explicit species lists retain their requirement that every request is eligible.
    subset = tmp_path / "subset.txt"
    subset.write_text("Alpha_plant\nBeta_sp-X\n")
    with pytest.raises(ValueError, match="requested species absent or below BUSCO threshold"):
        prepare(**tiny_inputs, outdir=out, species_list=subset)


def test_no_species_with_busco_still_fails_selection(tiny_inputs, tmp_path):
    rows = read_tsv(tiny_inputs["busco"])
    for row in rows:
        row["Species"] = "Unmatched " + row["Species"]
    write_tsv(tiny_inputs["busco"], list(rows[0]), rows)
    with pytest.raises(ValueError, match="no species passed selection"):
        prepare(**tiny_inputs, outdir=tmp_path / "metadata")


@pytest.mark.parametrize("problem", ["duplicate_run", "zero_total", "missing_quant", "missing_taxonomy"])
def test_bad_inputs_fail_early(tiny_inputs, tmp_path, problem):
    if problem == "duplicate_run":
        rows = read_tsv(tiny_inputs["metadata"])
        rows[1]["run"] = rows[0]["run"]
        write_tsv(tiny_inputs["metadata"], list(rows[0]), rows)
    elif problem == "zero_total":
        rows = read_tsv(tiny_inputs["busco"])
        rows[0]["busco_cds_total"] = 0
        write_tsv(tiny_inputs["busco"], list(rows[0]), rows)
    elif problem == "missing_quant":
        next(Path(tiny_inputs["quant_dir"]).rglob("A1_abundance.tsv")).unlink()
    else:
        Path(tiny_inputs["taxonomy_db"]).unlink()
    with pytest.raises(ValueError):
        prepare(**tiny_inputs, outdir=tmp_path / "out")


def test_snapshot_is_offline_and_does_not_overwrite(tiny_inputs, tmp_path):
    out = tmp_path / "snapshot.sqlite"
    snapshot(tiny_inputs["taxonomy_db"], out)
    with sqlite3.connect(out) as db:
        assert db.execute("SELECT count(*) FROM species").fetchone()[0] == 5
    with pytest.raises(ValueError, match="already exists"):
        snapshot(tiny_inputs["taxonomy_db"], out)


def test_annotation_variants(tmp_path):
    path = tmp_path / "annotations"
    path.write_text("# comment\n#query\tODB_OG\tevalue\ngene1\tOG1\t0\n")
    assert list(annotation_pairs(path)) == [("gene1", "OG1")]
    path.write_text("#cluster_id\tgene_id\tscore\nOG1\tgene1\t99\n")
    assert list(annotation_pairs(path)) == [("gene1", "OG1")]
    path.write_text("# No data found with taxids from example\n")
    assert list(annotation_pairs(path)) == []
    path.write_text("unknown\theader\ngene1\tOG1\n")
    with pytest.raises(ValueError):
        list(annotation_pairs(path))


def test_reference_checksum_detects_same_size_change(frozen_reference):
    marker = frozen_reference / "reference.json"
    verify(marker)
    data = frozen_reference / "odbmapper/v12/data/tiny.db"
    data.write_text("X" * data.stat().st_size)
    with pytest.raises(ValueError, match="checksum mismatch"):
        verify(marker)


@pytest.mark.parametrize("query", ["unrecognized_id", "Beta_sp-X_g1"])
def test_merge_rejects_queries_outside_chunk(tmp_path, query):
    proteins = tmp_path / "proteins"
    proteins.mkdir()
    rows = []
    for species in ["Alpha_plant", "Beta_sp-X"]:
        odb_species = species.replace("-", "_")
        rows.append({"species": species, "odb_species": odb_species})
        (proteins / f"{odb_species}_protein.fa").write_text(f">{species}_g1\nMK*\n")
    samples = tmp_path / "samples.tsv"
    write_tsv(samples, list(rows[0]), rows)
    manifests = tmp_path / "manifests"
    make(samples, proteins, manifests, chunk_size=1)
    chunks = tmp_path / "chunks"
    for i, value in enumerate([query, "Beta_sp-X_g1"]):
        label = f"chunk_{i:03d}"
        root = chunks / label
        root.mkdir(parents=True)
        annotated = root / f"{label}.og.annotations"
        annotated.write_text(f"#query\tODB_OG\n{value}\tOG1\n")
        write_json(root / "provenance.json", {"results": [file_record(annotated)]})
    with pytest.raises(ValueError, match="does not belong"):
        merge_odb(samples, manifests / "chunks.json", chunks, proteins,
                  tmp_path / "db.sqlite", tmp_path / "map.tsv", tmp_path / "qc.json")


def test_aggregation_ambiguity_and_run_preservation(tiny_inputs, tmp_path):
    meta = tmp_path / "metadata"
    prepare(**tiny_inputs, outdir=meta)
    database = tmp_path / "mappings.sqlite"
    with sqlite3.connect(database) as db:
        db.executescript("CREATE TABLE genes(query TEXT PRIMARY KEY, species TEXT); CREATE TABLE mappings(query TEXT, og TEXT);")
        for name in ["Alpha_plant", "Beta_sp-X"]:
            db.executemany("INSERT INTO genes VALUES (?, ?)", [(f"{name}_g{i}", name) for i in [1, 2, 3]])
            db.executemany("INSERT INTO mappings VALUES (?, ?)", [(f"{name}_g1", "OG1"), (f"{name}_g2", "OG2")])
    runs = tmp_path / "runs"
    for run in ["A1", "A2", "B1"]:
        aggregate(meta / "samples.tsv", run, database, runs / f"{run}.tsv", runs / f"{run}.qc.json")
    rows = read_tsv(runs / "A1.tsv")
    assert [float(r["tpm_sum"]) for r in rows] == [20, 30]
    assert [float(r["tpm"]) for r in rows] == [400000, 600000]
    report = json.loads((runs / "A1.qc.json").read_text())
    assert report["mapped_tpm_fraction"] == 0.5
    merge_tpm(meta / "samples.tsv", runs, tmp_path / "final")
    assert [r["run"] for r in read_tsv(tmp_path / "final/tpm_wide.tsv")] == ["A1", "A2", "B1"]
    with sqlite3.connect(database) as db:
        db.execute("INSERT INTO mappings VALUES ('Alpha_plant_g1', 'OG3')")
    with pytest.raises(ValueError, match="multiple OGs"):
        aggregate(meta / "samples.tsv", "A1", database, runs / "x.tsv", runs / "x.json")
    aggregate(meta / "samples.tsv", "A1", database, runs / "split.tsv", runs / "split.json", "split")
    assert sum(float(r["tpm_sum"]) for r in read_tsv(runs / "split.tsv")) == 50
    aggregate(meta / "samples.tsv", "A1", database, runs / "drop.tsv", runs / "drop.json", "drop")
    assert [r["orthogroup"] for r in read_tsv(runs / "drop.tsv")] == ["OG2"]


def test_odb_resume_is_bound_to_input_contents(fake_odb, frozen_reference, tmp_path, monkeypatch):
    protein = tmp_path / "protein.fa"
    protein.write_text(">Alpha_plant_g1\nMK*\n")
    manifest = tmp_path / "chunk_000.fs"
    manifest.write_text(str(protein) + "\n")
    fail = tmp_path / "fail_once"
    fail.touch()
    events = tmp_path / "events.txt"
    monkeypatch.setenv("FAKE_ODB_FAIL_ONCE", str(fail))
    monkeypatch.setenv("FAKE_ODB_LOG", str(events))
    args = dict(manifest=manifest, reference=frozen_reference / "reference.json", output_dir=tmp_path / "out",
                work_dir=tmp_path / "work", label="chunk_000", command=str(fake_odb), jobs=1,
                batch_size=1, min_free_gb=0, allow_nonlocal=True, keep_work=True)
    with pytest.raises(subprocess.CalledProcessError):
        run_odb(**args)
    run_odb(**args)
    old = json.loads((tmp_path / "out/provenance.json").read_text())["fingerprint"]
    assert len(set(events.read_text().splitlines())) == 1  # failed work resumed
    protein.write_text(">Alpha_plant_g1\nMKQ*\n")
    run_odb(**args)
    new = json.loads((tmp_path / "out/provenance.json").read_text())["fingerprint"]
    assert old != new
    assert len(set(events.read_text().splitlines())) == 2


@pytest.fixture
def existing_odb(tmp_path):
    root = tmp_path / "existing"
    root.mkdir()
    proteins = tmp_path / "original_proteins"
    proteins.mkdir()
    records, annotations = [], ["#query\tODB_OG\n"]
    for species in ["Alpha_plant", "Beta_sp-X"]:
        odb_species = species.replace("-", "_")
        protein = proteins / f"{odb_species}_protein.fa"
        protein.write_text("".join(f">{species}_g{i}\nMK*\n" for i in [1, 2, 3]))
        records.append({"species": species, "odb_species": odb_species, **file_record(protein)})
        for i in [1, 2]:
            annotations.extend([f"{species}_g{i}\tOG{i}\n"] * 2)
    annotation = root / "annotations.tsv"
    annotation.write_text("".join(annotations))
    write_json(root / "snapshot.json", {"schema_version": 1, "version": "v12", "node": 3193,
               "proteins": records, "annotations": {**file_record(annotation), "path": "annotations.tsv"}})
    return root, proteins


@pytest.mark.parametrize("problem,message", [
    ("protein", "protein differs"), ("annotation", "result changed"),
    ("node", "version/node differs"), ("species", "missing selected species"),
    ("query", "does not belong"),
])
def test_existing_odb_rejects_incompatible_inputs(existing_odb, tmp_path, problem, message):
    root, proteins = existing_odb
    snapshot = json.loads((root / "snapshot.json").read_text())
    rows = [{k: r[k] for k in ["species", "odb_species"]} for r in snapshot["proteins"]]
    samples = tmp_path / "samples.tsv"
    write_tsv(samples, list(rows[0]), rows)
    if problem == "protein":
        path = proteins / "Alpha_plant_protein.fa"
        path.write_text(path.read_text().replace("MK*", "MQ*"))
    elif problem in {"annotation", "query"}:
        path = root / "annotations.tsv"
        path.write_text(path.read_text() + "unknown_gene\tOG1\n")
        if problem == "query":
            snapshot["annotations"]["sha256"] = file_record(path)["sha256"]
    elif problem == "node":
        snapshot["node"] = 1
    else:
        snapshot["proteins"].pop()
    write_json(root / "snapshot.json", snapshot)
    with pytest.raises(ValueError, match=message):
        merge_odb(samples, "unused", "unused", proteins, tmp_path / "db.sqlite",
                  tmp_path / "map.tsv", tmp_path / "qc.json", existing=root)
    assert not (tmp_path / "db.sqlite").exists()


@pytest.mark.parametrize("reuse", [False, True])
def test_snakemake_end_to_end_and_incremental_rerun(tiny_inputs, fake_odb, frozen_reference, existing_odb, tmp_path, reuse, command_environment, workflow_project):
    snakemake = os.environ.get("SNAKEMAKE_BIN") or shutil.which("snakemake")
    seqkit = os.environ.get("SEQKIT_BIN") or shutil.which("seqkit")
    if not snakemake or not seqkit:
        pytest.skip("set SNAKEMAKE_BIN and SEQKIT_BIN to run the workflow integration test")
    import yaml
    config = {
        "analysis": "test", "inputs": {k: tiny_inputs[k] for k in ["metadata", "busco", "cds_dir", "quant_dir"]},
        "taxonomy": {"source": tiny_inputs["taxonomy_db"]},
        "odb": {"chunk_size": 1, "threads": 1, "batch_size": 1, "mem_gb": 3,
                "min_free_gb": 0, "allow_nonlocal": True, "keep_work": True},
    }
    configfile = tmp_path / "config.yaml"
    if reuse:
        config["odb"]["existing_results"] = str(existing_odb[0])
        # Leave the fixed reference absent; ODB must never execute in import mode.
    else:
        reference = workflow_project / "resources/orthodb/v12_3193"
        reference.parent.mkdir(parents=True)
        reference.symlink_to(frozen_reference, target_is_directory=True)
    configfile.write_text(yaml.safe_dump(config))
    # Supply standard command names as the rule environments do in production.
    env = command_environment({"python": sys.executable, "seqkit": seqkit,
                               "ODB-mapper": tmp_path / "absent_odb_command" if reuse else fake_odb})
    env["FAKE_ODB_LOG"] = str(tmp_path / "events.txt")
    base = [snakemake, "--snakefile", str(ROOT / "workflow/Snakefile"), "--configfile", str(configfile),
            "--cores", "2", "--resources", "mem_mb=16000"]
    def execute(extra=()):
        result = subprocess.run(base + list(extra), cwd=workflow_project, env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if result.returncode:
            logs = "\n".join(f"{p}:\n{p.read_text()}" for p in (tmp_path / "logs").rglob("*.log"))
            pytest.fail(result.stdout + "\n" + logs)
        return result.stdout
    execute(["--", "prepare"])
    taxonomy = workflow_project / "resources/taxonomy/taxa.sqlite"
    taxonomy_before = file_record(taxonomy), taxonomy.stat().st_mtime_ns
    assert json.loads(Path(str(taxonomy) + ".json").read_text())["method"] == "sqlite_backup"
    # Once prepared, the snapshot no longer depends on the bootstrap source.
    Path(tiny_inputs["taxonomy_db"]).unlink()
    plan = execute(["--dry-run"])
    assert "rule prepare_taxonomy:" not in plan
    assert plan.count("rule odb_map:") == (0 if reuse else 2)
    assert "rule prepare_odb_reference:" not in plan
    assert plan.count("mem_mb=3000") == (0 if reuse else 2)
    assert not any("<TBD>" in line for line in plan.splitlines() if "input:" in line)
    if not reuse:
        # Each MAP waits for both chunks to start: the CPU and memory budgets
        # permit overlap, and no separate ODB limit may serialize them.
        env["FAKE_ODB_BARRIER_COUNT"] = "2"
    execute()
    env.pop("FAKE_ODB_BARRIER_COUNT", None)
    out = tmp_path / "results/test"
    assert len(read_tsv(out / "tpm/tpm_wide.tsv")) == 3
    events = tmp_path / "events.txt"
    assert (len(events.read_text().splitlines()) if events.exists() else 0) == (0 if reuse else 2)
    assert json.loads((out / "odb/merged/merge_qc.json").read_text())["duplicate_pairs_removed"] == 4
    assert "Nothing to be done" in execute(["--dry-run"])
    abundance = Path(tiny_inputs["quant_dir"]) / "Alpha_plant/A1/A1_abundance.tsv"
    rows = read_tsv(abundance)
    rows[0]["tpm"] = 40
    write_tsv(abundance, list(rows[0]), rows)
    dry = execute(["--dry-run"])
    assert "rule aggregate_tpm" in dry
    assert "rule odb_map:" not in dry
    assert "rule merge_odb:" not in dry
    execute()
    assert (len(events.read_text().splitlines()) if events.exists() else 0) == (0 if reuse else 2)
    # Shrinking selection must rebuild the checkpoint DAG and omit stale runs/chunks.
    subset = tmp_path / "subset.txt"
    subset.write_text("Beta_sp-X\n")
    config["selection"] = {"species_list": str(subset)}
    configfile.write_text(yaml.safe_dump(config))
    execute()
    assert [r["run"] for r in read_tsv(out / "tpm/tpm_wide.tsv")] == ["B1"]
    if reuse:
        assert not events.exists()
        qc = json.loads((out / "odb/merged/merge_qc.json").read_text())
        assert qc["mode"] == "existing"
        assert qc["excluded_annotation_rows"] == 4
        assert qc["duplicate_pairs_removed"] == 2
    assert (file_record(taxonomy), taxonomy.stat().st_mtime_ns) == taxonomy_before


@pytest.mark.parametrize("memory,message", [
    ({"mem_gb": 0}, "odb.mem_gb must be a positive integer"),
    ({"mem_gb": True}, "odb.mem_gb must be a positive integer"),
    ({"mem_gb": 1.5}, "odb.mem_gb must be a positive integer"),
    ({"mem_mb": 256000}, "Replace odb.mem_mb with odb.mem_gb"),
])
def test_invalid_memory_config_fails_before_work(tmp_path, memory, message):
    snakemake = os.environ.get("SNAKEMAKE_BIN") or shutil.which("snakemake")
    if not snakemake:
        pytest.skip("Snakemake is not available")
    import yaml
    configfile = tmp_path / "memory.yaml"
    configfile.write_text(yaml.safe_dump({"odb": memory}))
    env = os.environ.copy()
    env["XDG_CACHE_HOME"] = str(tmp_path / "cache")
    result = subprocess.run([snakemake, "--snakefile", str(ROOT / "workflow/Snakefile"),
                             "--configfile", str(configfile), "--cores", "1", "--dry-run"],
                            cwd=ROOT, env=env, capture_output=True, text=True, timeout=30)
    assert result.returncode != 0
    assert message in result.stdout + result.stderr
