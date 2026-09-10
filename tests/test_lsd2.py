"""Native LSD2 integration and validation of time units, bounds and rounding."""
import json
import os
from pathlib import Path
import shutil

import pytest

from common import read_tsv, write_json, write_tsv
from date_phylogeny import (date, lsd_dates, newick_text, parse_lsd_report,
                           read_lsd_dates, validate_dated, validate_settings)
from infer_phylogeny import read_tree


@pytest.mark.parametrize("settings", [
    {"variance": True}, {"variance": 3}, {"variance": 1.0},
    {"variance_parameter": float("nan")}, {"variance_parameter": 0},
    {"variance": 0, "variance_parameter": 1}, {"numsites": True},
    {"numsites": 2**31}, {"smoothing": 10}, [],
])
def test_invalid_lsd2_settings(settings):
    with pytest.raises(ValueError):
        validate_settings(settings)


def source_tree(tmp_path):
    path = tmp_path / "source.nwk"
    path.write_text("((A:1,B:1)AB:1,(C:1,D:1)CD:1)ROOT;\n")
    return read_tree(path)


def native_tree(tmp_path, text):
    path = tmp_path / "native.nexus"
    path.write_text("#NEXUS\nBegin trees;\ntree 1 = " + text + "\nEnd;\n")
    return path


NATIVE = '((D[&date="0"]:25,C[&date="0"]:25)[&date="-25"]:75,(B[&date="0"]:33.3333,A[&date="0"]:33.3333)[&date="-33.3333"]:66.6667)[&date="-100"];'


def test_dates_reverse_bounds_and_fix_every_tip_at_present(tmp_path):
    text = lsd_dates(source_tree(tmp_path), [
        {"node": "ROOT", "taxa": "A,D", "min_age_ma": 90, "max_age_ma": 110},
        {"node": "AB", "taxa": "A,B", "min_age_ma": 30, "max_age_ma": 30}])
    assert text.splitlines() == ["6", "A 0", "B 0", "C 0", "D 0",
                                 "mrca(A,D) b(-110,-90)", "mrca(A,B) -30"]


def test_native_rounding_is_reconciled_with_bounds_and_ultrametric_export(tmp_path):
    source = source_tree(tmp_path)
    bounds = [{"node": "ROOT", "min_age_ma": 100, "max_age_ma": 100},
              {"node": "AB", "min_age_ma": 33.333333, "max_age_ma": 33.333333}]
    tree, ages, raw, adjustments = read_lsd_dates(native_tree(tmp_path, NATIVE), source, bounds)
    assert ages["AB"] == 33.333333
    assert raw["AB"] == 33.3333
    assert len(adjustments) == 1
    assert abs(adjustments[0]["adjustment_ma"]) < adjustments[0]["rounding_tolerance_ma"]
    output = tmp_path / "time.nwk"
    output.write_text(newick_text(tree))
    _, checked, height = validate_dated(output, source, bounds)
    assert height == pytest.approx(100, abs=1e-12)
    assert checked["AB"] == pytest.approx(33.333333, abs=1e-12)
    assert tree.common_ancestor(["A", "B"]).name == "AB"
    assert all(tree.get_distance(tree, tip) == pytest.approx(100, abs=1e-12) for tip in tree.leaves())


@pytest.mark.parametrize("text,reason", [
    (NATIVE.replace('D[&date="0"]', 'A[&date="0"]'), "species set"),
    (NATIVE.replace('D[&date="0"]', 'Z[&date="0"]'), "species set"),
    (NATIVE.replace('D[&date="0"]', 'D'), "missing/duplicate"),
    (NATIVE.replace('D[&date="0"]', 'D[&date="-1"]'), "tip dates"),
    (NATIVE.replace('date="-25"', 'date="nan"'), "nonfinite"),
    (NATIVE.replace('date="-25"', 'date="1"'), "future"),
    (NATIVE.replace('date="-25"', 'date="-125"'), "temporal"),
    (NATIVE.replace(':75', ':0.075'), "time units"),
    (NATIVE.replace(':75', ':-75'), "invalid time branches"),
    (NATIVE.replace('D[&date="0"]', 'Q[&date="0"]').replace('B[&date="0"]', 'D[&date="0"]').replace('Q[&date="0"]', 'B[&date="0"]'), "topology"),
])
def test_invalid_native_results_are_rejected(tmp_path, text, reason):
    with pytest.raises(ValueError, match=reason):
        read_lsd_dates(native_tree(tmp_path, text), source_tree(tmp_path), [])


def test_rounding_cannot_hide_a_real_calibration_violation(tmp_path):
    with pytest.raises(ValueError, match="beyond rounding"):
        read_lsd_dates(native_tree(tmp_path, NATIVE), source_tree(tmp_path),
                       [{"node": "AB", "min_age_ma": 34, "max_age_ma": 35}])


@pytest.mark.parametrize("body", ["", "rate 0.1, tMRCA -100 , objective function nan\n",
    "rate 0, tMRCA -100 , objective function 1\n", "rate 0.1, tMRCA 100 , objective function 1\n",
    "rate 0.1, tMRCA -100 , objective function 1\n" * 2])
def test_exit_zero_is_insufficient_without_a_valid_completed_fit(tmp_path, body):
    path = tmp_path / "report"
    path.write_text("LEAST-SQUARE METHODS - v.2.4.4\n" + body)
    with pytest.raises(ValueError):
        parse_lsd_report(path)


def test_missing_site_count_or_empty_process_cannot_publish(tmp_path):
    tree, provenance, calibration = tmp_path / "tree", tmp_path / "tree.json", tmp_path / "bounds"
    tree.write_text("(A:.1,(B:.08,(C:.04,D:.04):.04):.02);\n")
    write_tsv(calibration, ["taxa", "min_age_ma", "max_age_ma", "source"],
              [{"taxa": "A,B", "min_age_ma": 100, "max_age_ma": 100, "source": "synthetic"}])
    write_json(provenance, {"branch_length_unit": "substitutions_per_site", "outgroup": "A", "mean_gene_length": 250})
    with pytest.raises(ValueError, match="total_gene_sites"):
        date(tree, provenance, calibration, tmp_path / "out", "/bin/true")
    with pytest.raises(FileNotFoundError):
        date(tree, provenance, calibration, tmp_path / "out", "/bin/true", {"numsites": 750})
    assert not (tmp_path / "out/species_tree.dated.nwk").exists()
    assert not (tmp_path / "out/provenance.json").exists()
    assert list((tmp_path / "out/lsd2_runs").glob("*/run.log"))


@pytest.mark.parametrize("variance", [0, 1, 2])
def test_real_lsd2_preserves_topology_zero_branches_and_calibrations(tmp_path, variance):
    lsd2 = os.environ.get("LSD2_BIN") or shutil.which("lsd2")
    if not lsd2:
        pytest.skip("set LSD2_BIN for real dating integration")
    def clade(names):
        if len(names) == 1:
            return names[0] + (":0.07" if int(names[0][1:]) % 2 else ":0.03")
        k = len(names) // 2
        return f"({clade(names[:k])},{clade(names[k:])}):0.04"
    tree, provenance, calibration = tmp_path / "tree.nwk", tmp_path / "tree.json", tmp_path / "bounds.tsv"
    tree.write_text(clade([f"T{i:02d}" for i in range(32)]) + ";\n")
    original = read_tree(tree)
    original["T00"].dist = 0
    original.common_ancestor(["T00", "T01"]).dist = 0
    tree.write_text(newick_text(original)); before = tree.read_bytes()
    write_json(provenance, {"branch_length_unit": "substitutions_per_site", "outgroup": "T00", "total_gene_sites": 32000})
    write_tsv(calibration, ["taxa", "min_age_ma", "max_age_ma", "source"], [
        {"taxa": "T00,T31", "min_age_ma": 100, "max_age_ma": 100, "source": "synthetic root"},
        {"taxa": "T16,T31", "min_age_ma": 50, "max_age_ma": 60, "source": "synthetic internal"}])
    out = tmp_path / "dated tree"
    date(tree, provenance, calibration, out, lsd2, {"variance": variance})
    report = json.loads((out / "provenance.json").read_text())
    assert report["method"] == "LSD2 least-squares dating"
    assert report["confidence_intervals"] is False
    assert report["topology_preserved"] is True
    assert report["native_report"]["objective"] >= 0
    assert tree.read_bytes() == before
    dated = read_tree(out / "species_tree.dated.nwk", original.leaf_names())
    assert all(dated.get_distance(dated, leaf) == pytest.approx(100, rel=1e-10) for leaf in dated.leaves())
    node = dated.common_ancestor(["T16", "T31"])
    assert 50 - 1e-8 <= 100 - dated.get_distance(dated, node) <= 60 + 1e-8
    assert len(read_tsv(out / "node_ages.tsv")) == 31
    assert all(abs(float(r["adjustment_ma"])) <= float(r["rounding_tolerance_ma"])
               for r in read_tsv(out / "rounding_adjustments.tsv"))
