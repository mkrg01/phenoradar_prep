"""Dating safeguards and real random CV with heterogeneous branch rates."""
import json
import os
from pathlib import Path
import shutil

import pytest

from common import read_tsv, write_json, write_tsv
from date_phylogeny import (calibration_rows, date, final_objective, newick_text,
                           parse_cv, parse_prime, validate_dated, validate_settings)
from infer_phylogeny import read_tree


@pytest.mark.parametrize("settings", [
    {"smoothing": float("nan")}, {"cv_multiplier": 1}, {"cv_stop": 2000},
    {"replicates": 1}, {"numsites": True}, {"numsites": 2**31},
    {"optimization_iterations": 0}, {"cv_start": True}, {"unknown": 1},
])
def test_invalid_treepl_settings(settings):
    with pytest.raises(ValueError):
        validate_settings(settings)


def test_incomplete_cv_and_successful_exit_without_final_fit_are_errors(tmp_path):
    path = tmp_path / "output"
    path.write_text("chisq: (10) 3\n")
    with pytest.raises(ValueError, match="complete"):
        parse_cv(path, [10, 1])
    path.write_text("chisq: (10) nan\nchisq: (1) 2\n")
    with pytest.raises(ValueError, match="nonfinite"):
        parse_cv(path, [10, 1])
    path.write_text("opt = 1\noptad = 2\n")
    with pytest.raises(ValueError, match="all optimizers"):
        parse_prime(path)
    path.write_text("treepl_final_converged: 1\nafter opt calc: inf\n")
    with pytest.raises(ValueError, match="incomplete"):
        final_objective(path)
    path.write_text("thorough optimization hit 1000 iterations, breaking\n")
    with pytest.raises(ValueError, match="iteration limit"):
        final_objective(path)


def test_validation_checks_topology_bounds_and_ultrametricity(tmp_path):
    source, result = tmp_path / "source", tmp_path / "result"
    source.write_text("((A:1,B:1)AB:1,(C:1,D:1)CD:1)ROOT;\n")
    original = read_tree(source)
    bounds = [{"node": "ROOT", "min_age_ma": 90, "max_age_ma": 110}]
    result.write_text("((D:25,C:25):75,(B:40,A:40):60);\n")
    tree, ages, height = validate_dated(result, original, bounds)
    assert ages == {"ROOT": 100, "CD": 25, "AB": 40}
    assert height == 100
    assert tree.common_ancestor(["A", "B"]).name == "AB"
    result.write_text("((A:50,C:50):50,(B:50,D:50):50);\n")
    with pytest.raises(ValueError, match="topology"):
        validate_dated(result, original, bounds)
    result.write_text("((A:30,B:40):60,(C:25,D:25):75);\n")
    with pytest.raises(ValueError, match="ultrametric"):
        validate_dated(result, original, bounds)
    result.write_text("((A:1,B:1):1,(C:1,D:1):1);\n")
    with pytest.raises(ValueError, match="violates calibration"):
        validate_dated(result, original, bounds)


def test_total_site_count_cannot_be_replaced_by_mean_and_empty_run_cannot_publish(tmp_path):
    tree, provenance, calibration = tmp_path / "tree", tmp_path / "provenance", tmp_path / "bounds"
    tree.write_text("(A:.1,(B:.08,(C:.04,D:.04):.04):.02);\n")
    write_tsv(calibration, ["taxa", "min_age_ma", "max_age_ma", "source"],
              [{"taxa": "A,B", "min_age_ma": 100, "max_age_ma": 100, "source": "synthetic"}])
    write_json(provenance, {"branch_length_unit": "substitutions_per_site", "outgroup": "A", "mean_gene_length": 250})
    with pytest.raises(ValueError, match="total_gene_sites"):
        date(tree, provenance, calibration, tmp_path / "out", "/bin/true")
    with pytest.raises(ValueError, match="all optimizers"):
        date(tree, provenance, calibration, tmp_path / "out", "/bin/true", {"numsites": 750})
    assert not (tmp_path / "out/species_tree.dated.nwk").exists()
    assert not (tmp_path / "out/provenance.json").exists()


def test_real_random_cv_and_calibrations_with_heterogeneous_rates(tmp_path):
    treepl = os.environ.get("TREEPL_BIN") or shutil.which("treePL")
    if not treepl:
        pytest.skip("set TREEPL_BIN to the project treePL build")
    # 32 tips exercise multiple sampled tips per CV group and the sister check.
    from ete4 import Tree
    def clade(names):
        if len(names) == 1:
            return names[0] + (":0.07" if int(names[0][1:]) % 2 else ":0.03")
        k = len(names) // 2
        return f"({clade(names[:k])},{clade(names[k:])}):0.04"
    original = Tree(clade([f"T{i:02d}" for i in range(32)]) + ";", parser=1)
    original["T00"].dist = 0  # track native small-branch correction
    tree, provenance, calibration = tmp_path / "tree.nwk", tmp_path / "tree.json", tmp_path / "bounds.tsv"
    tree.write_text(newick_text(original)); before = tree.read_bytes()
    write_json(provenance, {"branch_length_unit": "substitutions_per_site", "outgroup": "T00", "total_gene_sites": 32000})
    write_tsv(calibration, ["taxa", "min_age_ma", "max_age_ma", "source"], [
        {"taxa": "T00,T31", "min_age_ma": 100, "max_age_ma": 100, "source": "synthetic root"},
        {"taxa": "T16,T31", "min_age_ma": 50, "max_age_ma": 60, "source": "synthetic internal"}])
    out = tmp_path / "dated tree"
    date(tree, provenance, calibration, out, treepl,
         {"cv_start": 10.0, "cv_stop": 1.0, "cv_replicates": 2, "replicates": 2})
    report = json.loads((out / "provenance.json").read_text())
    assert report["numsites"] == 32000
    assert report["numsites_is_effective_sample_size"] is False
    assert report["short_branches_adjusted"] == 1
    assert len(read_tsv(out / "cross_validation.tsv")) == 4
    assert len(read_tsv(out / "optimization_replicates.tsv")) == 2
    assert report["smoothing"] in (10, 1)
    assert report["confidence_intervals"] is False
    assert tree.read_bytes() == before
    dated = read_tree(out / "species_tree.dated.nwk", original.leaf_names())
    for leaf in dated.leaves():
        assert dated.get_distance(dated, leaf) == pytest.approx(100, rel=1e-7)
    assert report["root_age_ma"] == pytest.approx(100, rel=1e-7)
    assert "sampling 3 tips from 32 total tips 10 times" in next((out / "treepl_runs").glob("*/cv_01/run.log")).read_text()
