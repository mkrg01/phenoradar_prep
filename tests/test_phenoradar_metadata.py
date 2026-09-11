"""Consumer metadata retains the cohort and rejects conflicting annotations."""
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest
import yaml

from common import read_tsv, write_tsv
from phenoradar_metadata import fields, prepare, read_base, with_pairs


@pytest.fixture
def inputs(tmp_path):
    samples = tmp_path / "samples.tsv"
    metadata = tmp_path / "metadata.tsv"
    traits = tmp_path / "traits.tsv"
    rows = [{"species": species, "run": run, "taxid": taxid}
            for species, run, taxid in [("Plant_A", "A1", "1"), ("Plant_A", "A2", "1"),
                                       ("Plant_B-x", "B1", "2"), ("Plant_C", "C1", "3")]]
    write_tsv(samples, list(rows[0]), rows)
    write_tsv(metadata, ["species", "run", "family"],
              [{"species": row["species"], "run": row["run"], "family": "Family A" if row["species"] == "Plant_A" else ""}
               for row in rows])
    write_tsv(traits, ["species", "C4"], [{"species": "Plant A", "C4": "0"},
                                         {"species": "Plant B-x", "C4": "1"},
                                         {"species": "Outside_dataset", "C4": "1"}])
    return {"samples": samples, "metadata": metadata, "traits": traits,
            "output": tmp_path / "species_metadata.tsv"}


def test_prepare_preserves_missing_traits_and_collapses_consistent_runs(inputs):
    rows = prepare(**inputs)
    assert rows == [
        {"species": "Plant_A", "C4": "0", "contrast_pair_id": "", "family": "Family A"},
        {"species": "Plant_B-x", "C4": "1", "contrast_pair_id": "", "family": ""},
        {"species": "Plant_C", "C4": "", "contrast_pair_id": "", "family": ""},
    ]
    assert read_base(inputs["output"]) == rows
    assert inputs["output"].read_text().splitlines()[0] == "species\tC4\tcontrast_pair_id\tfamily"


def test_optional_traits_and_custom_binary_trait(inputs):
    assert all(row["C4"] == "" for row in prepare(**{**inputs, "traits": None}))
    write_tsv(inputs["traits"], ["species", "other"],
              [{"species": "Plant_A", "other": "0"}, {"species": "Plant_C", "other": "NA"}])
    rows = prepare(**inputs, trait="other")
    assert rows[0]["other"] == "0" and rows[2]["other"] == ""
    assert read_base(inputs["output"], "other") == rows
    write_tsv(inputs["traits"], ["species", "other"], [{"species": "Plant_A", "other": "C3"}])
    with pytest.raises(ValueError, match="0, 1 or missing"):
        prepare(**inputs, trait="other")


def test_explicit_missing_phenotype_source_fails(inputs):
    inputs["traits"].unlink()
    with pytest.raises(FileNotFoundError):
        prepare(**inputs)


def test_partial_pairs_preserve_base_and_replace_previous_pairs(inputs, tmp_path):
    rows = prepare(**inputs)
    pairs = tmp_path / "pairs.tsv"
    write_tsv(pairs, ["species", "C4", "contrast_pair_id", "role"],
              [{"species": name, "C4": value, "contrast_pair_id": "pair1", "role": "observed"}
               for name, value in [("Plant_A", "0"), ("Plant_B-x", "1")]])
    joined = with_pairs(read_base(inputs["output"]), pairs)
    assert [row["contrast_pair_id"] for row in joined] == ["pair1", "pair1", ""]
    assert all(row["contrast_pair_id"] == "" for row in rows)
    assert prepare(**inputs, pairs=pairs) == joined
    write_tsv(pairs, ["species", "C4", "contrast_pair_id"], [])
    assert with_pairs(joined, pairs) == rows


@pytest.mark.parametrize("damage, message", [
    ("duplicate_run", "duplicate run"),
    ("identity", "conflicting sample identity"),
    ("family", "conflicting family"),
    ("unknown_species", "unknown species"),
    ("missing_species", "species missing"),
    ("wrong_run_species", "run/species differs"),
    ("missing_run", "runs missing"),
])
def test_bad_selected_inputs_preserve_previous_output(inputs, damage, message):
    prepare(**inputs)
    before = inputs["output"].read_bytes()
    key = "samples" if damage in {"duplicate_run", "identity"} else "metadata"
    rows = read_tsv(inputs[key])
    if damage == "duplicate_run":
        rows[1]["run"] = rows[0]["run"]
    elif damage == "identity":
        rows[1]["taxid"] = "99"
    elif damage == "family":
        rows[1]["family"] = "Other family"
    elif damage == "unknown_species":
        rows[-1]["species"] = "Unknown"
    elif damage == "missing_species":
        rows.pop()
    elif damage == "wrong_run_species":
        rows[-1]["run"] = "B1"
    elif damage == "missing_run":
        rows.pop(1)
    write_tsv(inputs[key], list(rows[0]), rows)
    with pytest.raises(ValueError, match=message):
        prepare(**inputs)
    assert inputs["output"].read_bytes() == before


@pytest.mark.parametrize("damage, message", [
    ("unknown", "unknown species"), ("duplicate", "duplicate species"),
    ("trait", "trait differs"), ("malformed", "malformed TSV row"),
])
def test_bad_pair_merges_are_rejected(inputs, tmp_path, damage, message):
    rows = prepare(**inputs)
    pairs = tmp_path / "pairs.tsv"
    pair = {"species": "Plant_A", "C4": "0", "contrast_pair_id": "pair1"}
    if damage == "unknown":
        pair["species"] = "Outside_dataset"
    elif damage == "trait":
        pair["C4"] = "1"
    write_tsv(pairs, list(pair), [pair, pair] if damage == "duplicate" else [pair])
    if damage == "malformed":
        pairs.write_text("species\tC4\tcontrast_pair_id\nPlant_A\t0\n")
    with pytest.raises(ValueError, match=message):
        with_pairs(rows, pairs)


@pytest.mark.parametrize("trait", ["species", "family", "contrast_pair_id", "", " C4", "C4\n"])
def test_trait_column_cannot_conflict_with_output(trait):
    with pytest.raises(ValueError, match="non-reserved output column"):
        fields(trait)


def test_stored_base_rejects_duplicates_and_bad_headers(inputs):
    rows = prepare(**inputs)
    write_tsv(inputs["output"], fields(), rows + [rows[0]])
    with pytest.raises(ValueError, match="duplicate species"):
        read_base(inputs["output"])
    inputs["output"].write_text("species\tC4\tcontrast_pair_id\tfamily\tfamily\n")
    with pytest.raises(ValueError, match="duplicate TSV columns"):
        read_base(inputs["output"])


def test_explicit_metadata_target_backfills_completed_selection_without_raw_inputs(workflow_project, command_environment):
    snakemake = os.environ.get("SNAKEMAKE_BIN") or shutil.which("snakemake")
    if not snakemake:
        pytest.skip("Snakemake required")
    project = workflow_project
    metadata_dir = project / "results/test/metadata"
    raw = project / "unavailable_inputs"
    samples = [{"species": species, "scientific_name": species.replace("_", " "),
                "odb_species": species.replace("-", "_"), "taxid": str(i + 1), "run": f"R{i}",
                "cds": str(raw / f"{species}.faa"), "abundance": str(raw / f"R{i}.tsv")}
               for i, species in enumerate(["Plant_A", "Plant_B-x"])]
    write_tsv(metadata_dir / "samples.tsv", list(samples[0]), samples)
    write_tsv(metadata_dir / "metadata_high_busco.tsv", ["species", "run", "family", "C4"],
              [{"species": row["species"], "run": row["run"], "family": "Plantaceae", "C4": "1"}
               for row in samples])
    before = {path: (path.read_bytes(), path.stat().st_mtime_ns) for path in metadata_dir.iterdir()}
    cfg = project / "override.yaml"
    cfg.write_text(yaml.safe_dump({"analysis": "test", "inputs": {
        "metadata": str(raw / "metadata.tsv"), "species_trait": str(raw / "traits.tsv"),
        "busco": str(raw / "busco.tsv"), "cds_dir": str(raw / "cds"), "quant_dir": str(raw / "quant"),
    }}))
    root = Path(__file__).resolve().parents[1]
    process = subprocess.run(
        [snakemake, "--snakefile", str(root / "workflow/Snakefile"), "--configfile", str(cfg),
         "--cores", "1", "--", "phenoradar_metadata"],
        cwd=project, env=command_environment({"python": sys.executable}), text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    log = project / "logs/test/phenoradar_metadata.log"
    assert process.returncode == 0, process.stdout + (log.read_text() if log.exists() else "")
    assert "select_metadata" not in process.stdout
    assert all((path.read_bytes(), path.stat().st_mtime_ns) == previous for path, previous in before.items())
    assert read_base(metadata_dir / "species_metadata.tsv") == [
        {"species": row["species"], "C4": "", "contrast_pair_id": "", "family": "Plantaceae"}
        for row in samples
    ]
    assert {path.name for path in metadata_dir.iterdir()} == {
        "samples.tsv", "metadata_high_busco.tsv", "species_metadata.tsv",
    }
    assert not raw.exists()
    assert not (project / "resources").exists()
    assert not (project / "work").exists()
