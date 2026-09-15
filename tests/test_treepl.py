"""Dating safeguards and native leave-one-out CV with heterogeneous branch rates."""
import json
import os
from pathlib import Path
import shutil

import pytest

from common import read_tsv, write_json, write_tsv
from date_phylogeny import (date, final_objective, final_smoothing, newick_text, treepl_build,
                           parse_cv, parse_prime, validate_dated, validate_settings)
from infer_phylogeny import read_tree


@pytest.mark.parametrize("changed", ["patch", "binary", "commit"])
def test_patched_or_unverified_build_is_rejected(tmp_path, changed):
    from common import file_record
    from date_phylogeny import TREEPL_COMMIT, TREEPL_ARCHIVE_SHA256
    binary = tmp_path / "bin/treePL"
    binary.parent.mkdir()
    binary.write_bytes(b"synthetic executable")
    manifest = tmp_path / "share/treepl/build.json"
    record = {"commit": TREEPL_COMMIT, "archive_sha256": TREEPL_ARCHIVE_SHA256,
              "source_patch_applied": False, "executable_sha256": file_record(binary)["sha256"]}
    write_json(manifest, record)
    assert treepl_build(binary)[0] == record
    if changed == "patch":
        record["source_patch_applied"] = True
    elif changed == "commit":
        record["commit"] = "unknown"
    else:
        binary.write_bytes(b"different executable")
    write_json(manifest, record)
    with pytest.raises(ValueError, match="verified unmodified"):
        treepl_build(binary)


def test_native_rounding_is_validated_without_changing_branch_lengths(tmp_path):
    source, native = tmp_path / "source", tmp_path / "native"
    source.write_text("((A:1,B:1)AB:1,(C:1,D:1)CD:1)ROOT;\n")
    original = read_tree(source)
    bounds = [{"node": "ROOT", "min_age_ma": 100, "max_age_ma": 100},
              {"node": "AB", "min_age_ma": 40, "max_age_ma": 40}]
    native.write_text("((B:40,A:40):60,(D:25.000001,C:25):75);\n")
    before = native.read_bytes()
    tree, ages, height = validate_dated(native, original, bounds)
    assert ages["ROOT"] == pytest.approx(100, abs=1e-6) and ages["AB"] == 40
    assert all(tree.get_distance(tree, n) == pytest.approx(100, abs=1e-6) for n in tree.leaves())
    assert tree["D"].dist == 25.000001 and tree["C"].dist == 25
    assert native.read_bytes() == before
    native.write_text("((B:40,A:40):60,(D:25.001,C:25):75);\n")
    with pytest.raises(ValueError, match="beyond.*rounding"):
        validate_dated(native, original, bounds)


def test_upstream_report_needs_no_project_specific_convergence_marker(tmp_path):
    path = tmp_path / "log"
    path.write_text("after opt calc: 123.45\n")
    assert final_objective(path) == 123.45
    path.write_text("after opt calc: 123.45\nafter opt calc: 124\n")
    with pytest.raises(ValueError, match="incomplete"):
        final_objective(path)
    # Native CV logs a full fit for each grid point, followed by its final fit.
    path.write_text("smoothing:10\nafter opt calc: 123.45\nsmoothing:1\nafter opt calc: 124\n")
    assert final_objective(path, expected_fits=2) == 124
    assert final_smoothing(path, expected_fits=2) == 1
    with pytest.raises(ValueError, match="smoothing"):
        final_smoothing(path, expected_fits=3)


@pytest.mark.parametrize("settings", [
    {"smooth": float("nan")}, {"cvmultstep": 1}, {"cvstop": 2000},
    {"replicates": 3}, {"cv_replicates": 3}, {"smooth": False},
    {"lfiter": 0}, {"pliter": -1}, {"cviter": False}, {"cvstart": True},
    {"optimization_iterations": 2}, {"thorough": None}, {"unknown": 1},
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
    path.write_text("after opt calc: inf\n")
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


def test_total_site_count_cannot_be_replaced_by_mean_and_empty_run_cannot_publish(tmp_path, monkeypatch):
    tree, provenance, calibration = tmp_path / "tree", tmp_path / "provenance", tmp_path / "bounds"
    tree.write_text("(A:.1,(B:.08,(C:.04,D:.04):.04):.02);\n")
    write_tsv(calibration, ["taxa", "min_age_ma", "max_age_ma", "source"],
              [{"taxa": "A,B", "min_age_ma": 100, "max_age_ma": 100, "source": "synthetic"}])
    write_json(provenance, {"branch_length_unit": "substitutions_per_site", "outgroup": "A", "mean_gene_length": 250})
    with pytest.raises(ValueError, match="total_gene_sites"):
        date(tree, provenance, calibration, tmp_path / "out", "/bin/true")
    write_json(provenance, {"branch_length_unit": "substitutions_per_site", "outgroup": "A", "total_gene_sites": 750})
    monkeypatch.setattr("date_phylogeny.treepl_build", lambda command: ({}, {}))
    with pytest.raises(ValueError, match="all optimizers"):
        date(tree, provenance, calibration, tmp_path / "out", "/bin/true")
    assert not (tmp_path / "out/species_tree.dated.nwk").exists()
    assert not (tmp_path / "out/provenance.json").exists()


def test_real_native_cv_and_calibrations_with_heterogeneous_rates(tmp_path):
    treepl = os.environ.get("TREEPL_BIN") or shutil.which("treePL")
    if not treepl:
        pytest.skip("set TREEPL_BIN to the unmodified dating-environment build")
    # Rate heterogeneity, zero branches, and internal/root calibrations.
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
    out.mkdir()
    for name in ("optimization_replicates.tsv", "rounding_adjustments.tsv"):
        (out / name).write_text("superseded report\n")
    date(tree, provenance, calibration, out, treepl,
         {"cvstart": 10.0, "cvstop": 1.0})
    report = json.loads((out / "provenance.json").read_text())
    assert report["numsites"] == 32000
    assert report["numsites_is_effective_sample_size"] is False
    assert report["short_branches_adjusted"] == 1
    assert len(read_tsv(out / "cross_validation.tsv")) == 2
    assert not (out / "optimization_replicates.tsv").exists()
    assert not (out / "rounding_adjustments.tsv").exists()
    assert len(report["commands"]) == 2
    assert [Path(c["cwd"]).name for c in report["commands"]] == ["prime", "fit"]
    assert report["smoothing"] in (10, 1)
    assert report["confidence_intervals"] is False
    assert tree.read_bytes() == before
    dated = read_tree(out / "species_tree.dated.nwk", original.leaf_names())
    for leaf in dated.leaves():
        assert dated.get_distance(dated, leaf) == pytest.approx(100, rel=1e-7)
    assert report["root_age_ma"] == pytest.approx(100, rel=1e-7)
    assert report["cv_method"] == "native leave-one-out"
    assert report["smoothing_selection"] == "treePL native CV"
    assert report["smoothing"] == final_smoothing(out / "treepl.log", expected_fits=3)
    assert report["build"]["source_patch_applied"] is False
    assert "\ncv\n" in (out / "treepl.config.txt").read_text()
    assert "randomcv" not in (out / "treepl.config.txt").read_text()
    native = read_tree(out / "treepl.dated.nwk")
    branches = lambda t: {tuple(sorted(n.leaf_names())): n.dist for n in t.traverse() if not n.is_root}
    assert branches(native) == branches(dated)
    assert report["native_time_branch_lengths_preserved"] is True
    assert set(read_tsv(out / "node_ages.tsv")[0]) == {"node", "age_ma"}
