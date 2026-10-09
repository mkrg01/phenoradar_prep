"""Validate AMALGKIT taxids before freezing inputs or invoking GeneGalleon."""
import subprocess

import pytest
import yaml

from common import read_tsv, write_tsv
from dataset import load, plan, prepare, status, submit, worker
from dataset_assets import identities
from test_datasets import dataset_project, fake_genegalleon, native_events


RANK_FIELDS = ("taxid_domain", "taxid_kingdom", "taxid_phylum", "taxid_class",
               "taxid_order", "taxid_family", "taxid_genus", "taxid_species")


def metadata_row(**values):
    return dict(scientific_name="Aldrovanda vesiculosa", run="SRR1979677",
                taxid="173386", **values)


def test_rank_ids_are_normalized_without_changing_other_metadata(tmp_path):
    path = tmp_path / "metadata.tsv"
    row = metadata_row(**dict(zip(RANK_FIELDS, (
        "2759.0", "33090.00", "35493", "03398.0", " 3524.00 ",
        "", "9007199254740993.0", "173386.0"))),
        host_taxid="999.0", total_bases="123.50", sample_attribute_taxid="custom.0")
    write_tsv(path, list(row), [row])
    original = path.read_bytes()

    fields, items = identities(path)
    normalized = items[0]["row"]
    assert [normalized[field] for field in RANK_FIELDS] == [
        "2759", "33090", "35493", "3398", "3524", "", "9007199254740993", "173386"]
    assert all(normalized[field] == value for field, value in row.items() if field not in RANK_FIELDS)
    assert fields == [*row, "species_id", "analysis_sample_id"]
    assert path.read_bytes() == original


@pytest.mark.parametrize("value", ["", " \t "])
def test_blank_species_rank_remains_missing(tmp_path, value):
    path = tmp_path / "metadata.tsv"
    row = metadata_row(taxid_species=value)
    write_tsv(path, list(row), [row])
    assert identities(path)[1][0]["row"]["taxid_species"] == ""


def test_species_rank_keeps_its_own_taxid(tmp_path):
    path = tmp_path / "metadata.tsv"
    row = metadata_row(taxid_species="4360.00")
    write_tsv(path, list(row), [row])
    normalized = identities(path)[1][0]["row"]
    assert normalized["taxid"] == "173386"
    assert normalized["taxid_species"] == "4360"


@pytest.mark.parametrize("field", RANK_FIELDS)
def test_fractional_rank_ids_are_rejected_with_sample_context(tmp_path, field):
    path = tmp_path / "metadata.tsv"
    row = metadata_row(**{field: "173386.5"})
    write_tsv(path, list(row), [row])
    with pytest.raises(ValueError) as error:
        identities(path)
    assert field in str(error.value)
    assert "SRR1979677" in str(error.value)
    assert "173386.5" in str(error.value)


@pytest.mark.parametrize("value", ["0", "0.0", "-1", "NaN", "inf", "bad", "1e5", "173386.", "１２３"])
def test_invalid_species_rank_ids_are_rejected(tmp_path, value):
    path = tmp_path / "metadata.tsv"
    row = metadata_row(taxid_species=value)
    write_tsv(path, list(row), [row])
    with pytest.raises(ValueError, match="positive integer taxid_species"):
        identities(path)


@pytest.mark.parametrize("value,expected", [
    ("173386", "173386"), ("173386.0", "173386"), ("173386.00", "173386"),
    (" 0173386.00 ", "173386"), ("9007199254740993.0", "9007199254740993")])
def test_required_taxid_accepts_integer_decimal_notation(tmp_path, value, expected):
    path = tmp_path / "metadata.tsv"
    row = metadata_row(taxid_species="173386.0")
    row["taxid"] = value
    write_tsv(path, list(row), [row])
    original = path.read_bytes()
    normalized = identities(path)[1][0]["row"]
    assert normalized["taxid"] == expected
    assert normalized["taxid_species"] == "173386"
    assert path.read_bytes() == original


@pytest.mark.parametrize("value", ["", " \t ", "0", "0.0", "-1", "173386.5", "NaN", "inf", "bad"])
def test_required_taxid_rejects_missing_or_invalid_values(tmp_path, value):
    path = tmp_path / "metadata.tsv"
    row = metadata_row(taxid_species="173386.0")
    row["taxid"] = value
    write_tsv(path, list(row), [row])
    with pytest.raises(ValueError, match="positive taxid required"):
        identities(path)


def test_different_taxid_notations_are_compared_by_normalized_value(tmp_path):
    path = tmp_path / "metadata.tsv"
    rows = [metadata_row(), metadata_row()]
    rows[1].update(run="SRR2", taxid="173386.0")
    write_tsv(path, list(rows[0]), rows)
    assert {item["row"]["taxid"] for item in identities(path)[1]} == {"173386"}
    rows[1]["taxid"] = "173387.0"
    write_tsv(path, list(rows[0]), rows)
    with pytest.raises(ValueError, match="conflicting taxids"):
        identities(path)


def test_plan_rejects_invalid_rank_before_preparing_a_build(dataset_project):
    root = dataset_project
    source = root / "input/rank_taxids.tsv"
    row = metadata_row(taxid_species="173386.5")
    write_tsv(source, list(row), [row])
    config = root / "config/build.yaml"
    with pytest.raises(ValueError, match="positive integer taxid_species"):
        plan(root, config, source)
    with pytest.raises(ValueError, match="positive integer taxid_species"):
        prepare(root, "invalid_rank", config, source)
    assert not (root / "results/invalid_rank").exists()


def test_build_preserves_source_and_passes_normalized_rank_to_worker(dataset_project):
    root = dataset_project
    fake_genegalleon(root)
    source = root / "input/rank_taxids.tsv"
    row = metadata_row(taxid_species="173386.0", taxid_family="4360.00")
    row["taxid"] = "173386.00"
    write_tsv(source, list(row), [row])
    original = source.read_bytes()
    build = prepare(root, "normalized_ranks", root / "config/build.yaml", source)
    manifest = load(build, check_code=True)
    frozen = read_tsv(build / "metadata.tsv")
    assert frozen == [manifest["items"][0]["row"]]
    assert frozen[0]["taxid"] == "173386"
    assert frozen[0]["taxid_species"] == "173386"
    assert frozen[0]["taxid_family"] == "4360"
    assert (build / "source_metadata.tsv").read_bytes() == original

    submit(build, until="quant", dry_run=True)
    staged = build / "work/genegalleon/Aldrovanda_vesiculosa_SRR1979677/input/amalgkit_metadata/Aldrovanda_vesiculosa_metadata.tsv"
    assert read_tsv(staged) == frozen
    worker(build, 1)
    assert status(build)[0]["quant"] == "reuse"
    assert source.read_bytes() == original
    assert (build / "source_metadata.tsv").read_bytes() == original


def test_old_failed_build_checkpoints_survive_rank_normalization(dataset_project, monkeypatch):
    root = dataset_project
    fake_genegalleon(root)
    source = root / "input/rank_taxids.tsv"
    row = metadata_row(taxid_species="173386.0")
    write_tsv(source, list(row), [row])
    config = root / "config/build.yaml"
    # Model a saved build made before optional rank IDs were normalized.
    with monkeypatch.context() as old_behavior:
        old_behavior.setattr("dataset_assets.normalize_taxonomy_ids", lambda row: dict(row))
        previous = prepare(root, "old_ranks", config, source)
    previous_manifest = (previous / "build.json").read_bytes()
    submit(previous, until="quant", dry_run=True)
    with monkeypatch.context() as failed_quant:
        failed_quant.setenv("FAKE_GG_FAIL_QUANT", "1")
        with pytest.raises(subprocess.CalledProcessError):
            worker(previous, 1)
    assert [status(previous)[0][stage] for stage in ("assembly", "busco", "quant")] == [
        "reuse", "reuse", "pending"]

    cfg = yaml.safe_load(config.read_text())
    cfg["reuse_from"] = str(previous)
    config.write_text(yaml.safe_dump(cfg))
    row["taxid"] = "173386.00"
    write_tsv(source, list(row), [row])
    updated = prepare(root, "new_ranks", config, source)
    assert [status(updated)[0][stage] for stage in ("assembly", "busco", "quant")] == [
        "reuse", "reuse", "pending"]
    assert read_tsv(updated / "metadata.tsv")[0]["taxid_species"] == "173386"
    assert read_tsv(updated / "metadata.tsv")[0]["taxid"] == "173386"
    submit(updated, until="quant", dry_run=True)
    worker(updated, 1)
    assert status(updated)[0]["quant"] == "reuse"
    env = native_events(updated)[0]["env"]
    assert env["GG_TRANSCRIPTOME_RUN_ASSEMBLY"] == "0"
    assert env["GG_TRANSCRIPTOME_RUN_BUSCO_LONGEST_CDS"] == "0"
    assert env["GG_TRANSCRIPTOME_RUN_AMALGKIT_QUANT"] == "1"
    assert (previous / "build.json").read_bytes() == previous_manifest
    assert read_tsv(previous / "metadata.tsv")[0]["taxid_species"] == "173386.0"
