"""Resolve BUSCO full tables by exact species ID without suffix configuration."""
import gzip
from pathlib import Path

import pytest

from busco_phylogeny import busco_full_path, busco_table, plan
from common import read_tsv, write_tsv


LINEAGE = "embryophyta_odb12"
FULL_TABLE = (f"# The lineage dataset is: {LINEAGE}\n"
              "# Busco id\tStatus\tSequence\tScore\tLength\n"
              "1at1\tComplete\tg1\t100\t100\n")


@pytest.mark.parametrize("filename", [
    "Plant_alpha.busco.full.tsv", "Plant_alpha.tsv",
    "Plant_alpha/full_table.tsv", f"Plant_alpha/run_{LINEAGE}/full_table.tsv",
])
@pytest.mark.parametrize("compressed", [False, True])
def test_full_tables_are_detected_and_read(tmp_path, filename, compressed):
    path = tmp_path / (filename + (".gz" if compressed else ""))
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if compressed else open
    with opener(path, "wt") as handle:
        handle.write(FULL_TABLE)
    found = busco_full_path(tmp_path, "Plant_alpha", LINEAGE)
    assert found == path
    assert busco_table(found, LINEAGE)[1] == {"1at1"}


@pytest.mark.parametrize("extra", [
    "Plant_alpha.busco.full.tsv.gz", "Plant_alpha.tsv", "Plant_alpha/full_table.tsv",
])
def test_ambiguous_tables_name_all_candidates(tmp_path, extra):
    first = tmp_path / "Plant_alpha.busco.full.tsv"
    first.write_text(FULL_TABLE)
    second = tmp_path / extra
    second.parent.mkdir(parents=True, exist_ok=True)
    second.write_text(FULL_TABLE)
    with pytest.raises(ValueError, match="ambiguous BUSCO full tables for Plant_alpha") as error:
        busco_full_path(tmp_path, "Plant_alpha", LINEAGE)
    assert str(first) in str(error.value)
    assert str(second) in str(error.value)


@pytest.mark.parametrize("unrelated", [
    "Plant_alpha_extra.busco.full.tsv", "Other_species.busco.full.tsv",
    "Plant_alpha/run_embryophyta_odb10/full_table.tsv", "summary.tsv",
])
def test_missing_table_does_not_use_other_species_or_lineages(tmp_path, unrelated):
    path = tmp_path / unrelated
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(FULL_TABLE)
    with pytest.raises(ValueError, match="BUSCO full table missing for Plant_alpha"):
        busco_full_path(tmp_path, "Plant_alpha", LINEAGE)


def test_plan_uses_cds_paths_from_samples_and_tracks_discovered_tables(tmp_path):
    rows = []
    for species in ["Plant_alpha", "Plant_beta", "Plant_gamma", "Plant_delta"]:
        cds = tmp_path / f"{species}.original.fna"
        cds.write_text(">g1\n" + "ATG" * 100 + "\n")
        (tmp_path / f"{species}.tsv").write_text(FULL_TABLE)
        rows.append({"species": species, "cds": str(cds)})
    samples = tmp_path / "samples.tsv"
    write_tsv(samples, ["species", "cds"], rows)
    settings = {"busco_full_dir": str(tmp_path), "outgroup": "Plant_alpha",
                "lineage": LINEAGE, "min_taxa": 4, "max_markers": 1}
    out = tmp_path / "plan"
    plan(samples, out, settings)
    expected = {row["species"]: row["cds"] for row in rows}
    for row in read_tsv(out / "species.tsv"):
        assert row["sequences"] == expected[row["species"]]
        assert Path(row["busco_table"]) == tmp_path / f'{row["species"]}.tsv'
    (tmp_path / "Plant_alpha.tsv").unlink()
    moved = tmp_path / "Plant_alpha.busco.full.tsv.gz"
    with gzip.open(moved, "wt") as handle:
        handle.write(FULL_TABLE)
    plan(samples, out, settings)
    assert Path(read_tsv(out / "species.tsv")[0]["busco_table"]) == moved
