"""Fixed contrast anchors expand to disjoint clades and retain full sample identity."""
import json
from pathlib import Path

import pytest
from ete4 import Tree

from common import file_record, read_tsv, sha256, write_json, write_tsv
from consumer_metadata import validate as validate_consumer
from expanded_pair_clades import COLUMN, annotate, export_analysis, validate_result

NEWICK = "(O:1,(((A:1,B:1):1,X:1):1,((C:1,D:1):1,Y:1):1):1,Z:1);"


def metadata():
    rows = []
    for name in "ABXCDYOZ":
        rows.append(dict(species=name, species_id=name, C4="1" if name in "BDY" else "0",
                         role="observed", group=name, representative=name, is_representative="1",
                         n_species_in_group="2" if name == "B" else "1", n_samples_in_group="2" if name in "AB" else "1",
                         contrast_pair_id="1" if name in "AB" else "2" if name in "CD" else ""))
    rows += [{**rows[0], "species": "A_run2", "is_representative": "0"},
             {**rows[1], "species": "B2", "species_id": "B2", "is_representative": "0"},
             dict(species="U", species_id="U", C4="", role="missing_trait", group="", representative="",
                  is_representative="0", n_species_in_group="", n_samples_in_group="", contrast_pair_id="")]
    return rows


@pytest.mark.parametrize("newick", [NEWICK, "(Z:1,((Y:1,(D:1,C:1):1):1,(X:1,(B:1,A:1):1):1):1,O:1);"])
def test_maximal_disjoint_membership_preserves_pairs_traits_and_biological_counts(newick):
    source = metadata()
    original = [dict(row) for row in source]
    rows, clades, report = annotate(Tree(newick, parser=0), source, "C4")
    assert source == original
    assert {r["species"]: r[COLUMN] for r in rows} == {
        **dict.fromkeys(["A", "B", "X", "A_run2", "B2"], "1"),
        **dict.fromkeys(["C", "D", "Y"], "2"), **dict.fromkeys(["O", "Z", "U"], "")}
    assert all({k: v for k, v in row.items() if k != COLUMN} == before for row, before in zip(rows, original))
    assert report["expanded_pair_clades"] == 2 and report["expanded_pair_species"] == 7
    assert report["expanded_pair_samples"] == 8
    assert clades[0]["n_samples"] == 5 and clades[0]["n_species"] == 4
    assert clades[0]["n_species_a"] == clades[0]["n_species_b"] == 2
    assert clades[0]["n_representatives"] == 3


@pytest.mark.parametrize("pair_count", [0, 1])
def test_zero_pairs_stay_empty_and_one_pair_can_expand_to_the_root(pair_count):
    source = metadata()
    for row in source:
        if not pair_count or row["contrast_pair_id"] == "2":
            row["contrast_pair_id"] = ""
    rows, clades, _ = annotate(Tree(NEWICK, parser=0), source, "C4")
    assert len(clades) == pair_count
    assert {r[COLUMN] for r in rows if r["C4"]} == {"1" if pair_count else ""}
    assert rows[-1][COLUMN] == ""
    empty, clades, report = annotate(None, [source[-1]], "C4")
    assert empty[0][COLUMN] == "" and not clades and report["expanded_pair_species"] == 0


@pytest.mark.parametrize("problem", ["non_clade", "one_tip", "unknown_representative", "missing_trait", "duplicate"])
def test_invalid_anchors_or_membership_are_rejected(problem):
    rows = metadata()
    if problem == "non_clade":
        rows[1]["contrast_pair_id"] = ""
        rows[2]["contrast_pair_id"] = "1"
    elif problem == "one_tip": rows[1]["contrast_pair_id"] = ""
    elif problem == "unknown_representative": rows[-2]["representative"] = "unknown"
    elif problem == "missing_trait": rows[-2]["C4"] = ""
    else: rows.append(dict(rows[0]))
    with pytest.raises(ValueError): annotate(Tree(NEWICK, parser=0), rows, "C4")


@pytest.fixture
def saved_analysis(tmp_path):
    source = tmp_path / "analysis"
    rows = metadata()
    write_json(source / "analysis.json", {"pipeline": {"trait": "C4", "phylogeny": {
        "trees": ["representatives"], "contrast_pairs": {"enabled": True}}}})
    write_json(source / "completed.json", {"historical": "unchanged"})
    write_tsv(source / "metadata/species_metadata.tsv", ["species", "C4", "contrast_pair_id", "family"],
              [{**{k: r[k] for k in ["species", "C4", "contrast_pair_id"]}, "family": "Plant family"} for r in rows])
    write_tsv(source / "metadata/samples.tsv", ["species", "scientific_name", "run"],
              [dict(species=r["species"], scientific_name=f"Plant {r['species']}", run=f"Run_{r['species']}") for r in rows])
    contrast = source / "phylogeny/representatives/contrast"
    write_tsv(contrast / "species_metadata.tsv", list(rows[0]), rows)
    (contrast / "summary_tree.nwk").write_text(NEWICK)
    write_tsv(contrast / "contrast_pairs.tsv", ["contrast_pair_id", "representative_a", "representative_b"],
              [dict(contrast_pair_id="1", representative_a="A", representative_b="B"),
               dict(contrast_pair_id="2", representative_a="C", representative_b="D")])
    write_json(contrast / "summary.json", {"trait": "C4", "summary_species": 8, "contrast_pairs": 2})
    return source


def test_saved_analysis_export_is_atomic_preserves_sources_and_validates_provenance(saved_analysis):
    original = {path: path.read_bytes() for path in saved_analysis.rglob("*") if path.is_file()}
    output = saved_analysis.parent / "export"
    result = export_analysis(saved_analysis, output)
    assert result["branches"]["representatives"]["expanded_pair_species"] == 7
    ready = validate_consumer(output / "metadata/species_metadata.tsv")
    assert sum(bool(r[COLUMN]) for r in ready) == 8
    assert sum(bool(r["contrast_pair_id"]) for r in ready) == 6
    assert all(r["family"] == "Plant family" for r in ready)
    contrast = output / "phylogeny/representatives/contrast"
    assert validate_result(contrast, "C4") == result["branches"]["representatives"]
    for stem in ["expanded_pair_tree", "expanded_pair_tree_id"]:
        assert (contrast / f"{stem}.pdf").read_bytes().startswith(b"%PDF")
        assert "Expanded pair clade" in (contrast / f"{stem}.svg").read_text()
    info = json.loads((output / "provenance.json").read_text())
    assert info["reestimated"] is info["pairs_reselected"] is False
    assert all(file_record(r["path"]) == r for r in info["sources"])
    assert all(not Path(r["path"]).is_absolute() and sha256(output / r["path"]) == r["sha256"] for r in info["outputs"])
    assert {path: path.read_bytes() for path in original} == original
    with pytest.raises(ValueError, match="already exists"): export_analysis(saved_analysis, output)


def test_failed_saved_export_publishes_nothing(saved_analysis):
    path = saved_analysis / "phylogeny/representatives/contrast/contrast_pairs.tsv"
    path.write_text(path.read_text().replace("\tA\tB", "\tA\tX"))
    with pytest.raises(ValueError, match="representative differs"):
        export_analysis(saved_analysis)
    assert not (saved_analysis / "expanded_pair_clades").exists()
    assert not list(saved_analysis.glob(".expanded-pairs-*"))
