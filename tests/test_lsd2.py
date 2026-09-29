"""Native LSD2 dating, basal-clade retention, and calibration safeguards."""
import json
import os
from pathlib import Path
import shutil

import pytest
from ete4 import Tree

from common import file_record, read_tsv, write_json, write_tsv
from date_phylogeny import (LSD2_COMMIT, LSD2_ARCHIVE_SHA256, date, dates_text,
                           lsd2_build, native_result, normalize_numerical_zeros, validate_dated, validate_settings)
from infer_phylogeny import astral, read_tree
from phylogeny_outgroup import root_on_outgroup, validate_root
from phylogeny_root import prepare_root


@pytest.fixture
def lsd2():
    binary = os.environ.get("LSD2_BIN") or shutil.which("lsd2")
    if not binary:
        pytest.skip("set LSD2_BIN to the verified dating-environment build")
    return binary


def inputs(tmp_path, text="((O1:.02,O2:.02):.03,(A:.04,(B:.02,C:.02):.02):.03);", bounds=(100, 100)):
    tree, qc, cal = tmp_path / "tree.nwk", tmp_path / "tree.json", tmp_path / "bounds.tsv"
    tree.write_text(text + "\n")
    write_json(qc, {"branch_length_unit": "substitutions_per_site", "outgroup": ["O1", "O2"],
                    "total_gene_sites": 32000})
    write_tsv(cal, ["taxa", "min_age_ma", "max_age_ma", "source"], [
        {"taxa": "O1,A", "min_age_ma": bounds[0], "max_age_ma": bounds[1], "source": "synthetic root"}])
    return tree, qc, cal


@pytest.mark.parametrize("changed", ["patch", "binary", "commit"])
def test_unverified_build_is_rejected(tmp_path, changed):
    binary = tmp_path / "bin/lsd2"
    binary.parent.mkdir()
    binary.write_bytes(b"synthetic executable")
    manifest = tmp_path / "share/lsd2/build.json"
    record = {"commit": LSD2_COMMIT, "archive_sha256": LSD2_ARCHIVE_SHA256,
              "source_patch_applied": False, "executable_sha256": file_record(binary)["sha256"]}
    write_json(manifest, record)
    assert lsd2_build(binary)[0] == record
    if changed == "patch":
        record["source_patch_applied"] = True
    elif changed == "commit":
        record["commit"] = "unknown"
    else:
        binary.write_bytes(b"different executable")
    write_json(manifest, record)
    with pytest.raises(ValueError, match="verified unmodified"):
        lsd2_build(binary)


@pytest.mark.parametrize("settings", [{"variance": True}, {"variance": 3}, {"variance_parameter": 0},
    {"variance_parameter": float("nan")}, {"variance": 0, "variance_parameter": .1}, {"smooth": 1},
    {"remove_outgroup": True}, {"root": "a"}, {"collapse": True}])
def test_invalid_settings(settings):
    with pytest.raises(ValueError):
        validate_settings(settings)


def test_positive_age_bounds_become_negative_dates_in_reverse_order():
    assert dates_text([{"taxa": "A,B", "min_age_ma": 90, "max_age_ma": 110},
                       {"taxa": "C,D", "min_age_ma": 40, "max_age_ma": 40}]) == (
        "2\nmrca(A,B) b(-110,-90)\nmrca(C,D) -40\n")


def test_incomplete_native_fit_cannot_publish(tmp_path, monkeypatch):
    tree, qc, cal = inputs(tmp_path)
    monkeypatch.setattr("date_phylogeny.lsd2_build", lambda _: ({}, {}))
    with pytest.raises((ValueError, FileNotFoundError)):
        date(tree, qc, cal, tmp_path / "dated", "/bin/true")
    assert not (tmp_path / "dated/provenance.json").exists()
    assert not (tmp_path / "dated/species_tree.dated.nwk").exists()


def test_native_result_rejects_nonfinite_or_incomplete_fit(tmp_path):
    result = tmp_path / "result"
    for value in ("", " rate nan, tMRCA -100 , objective function 0\nTOTAL ELAPSED TIME: 0\n",
                  " rate 0.001, tMRCA -100 , objective function nan\nTOTAL ELAPSED TIME: 0\n"):
        result.write_text(value)
        with pytest.raises(ValueError):
            native_result(result, 1)


def test_native_rounding_topology_and_bounds(tmp_path):
    source, result = tmp_path / "source", tmp_path / "native"
    source.write_text("((A:1,B:1)AB:1,(C:1,D:1)CD:1)ROOT;")
    original = read_tree(source)
    bounds = [{"node": "ROOT", "min_age_ma": 100, "max_age_ma": 100}]
    result.write_text("((B:33.3333,A:33.3333):66.6667,(D:60,C:60):40);")
    before = result.read_bytes()
    dated, ages, height = validate_dated(result, original, bounds)
    assert height == pytest.approx(100)
    assert dated["B"].dist == 33.3333 and result.read_bytes() == before
    for text, error in [("((A:50,C:50):50,(B:50,D:50):50);", "topology"),
                        ("((A:30,B:40):60,(C:25,D:25):75);", "ultrametric"),
                        ("((A:1,B:1):1,(C:1,D:1):1);", "calibration")]:
        result.write_text(text)
        with pytest.raises(ValueError, match=error):
            validate_dated(result, original, bounds)


@pytest.mark.parametrize("variance", [0, 1, 2])
def test_native_retains_basal_tips_and_estimates_root_split(tmp_path, lsd2, variance):
    tree, qc, cal = inputs(tmp_path)
    before = tree.read_bytes()
    out = tmp_path / "dated tree"
    date(tree, qc, cal, out, lsd2, {"variance": variance})
    result = read_tree(out / "species_tree.dated.nwk", ["O1", "O2", "A", "B", "C"])
    validate_root(result, ["O1", "O2"])
    basal = result.common_ancestor(["O1", "O2"])
    assert basal.dist == pytest.approx(66.6667)
    assert result.common_ancestor(["A", "B"]).dist == pytest.approx(33.3333)
    assert all(result.get_distance(result, n) == pytest.approx(100, abs=0.00016) for n in result.leaves())
    report = json.loads((out / "provenance.json").read_text())
    assert report["method"] == "LSD2 least squares" and report["root_position_reestimated"]
    assert report["outgroups_retained"] and report["unique_time_scale"]
    assert report["confidence_intervals"] is False and report["numsites"] == 32000
    assert tree.read_bytes() == before
    # Output time edges are exactly the native NEXUS edges, with node names restored.
    native = read_tree(out / "lsd2.dated.nwk")
    edges = lambda t: {tuple(sorted(n.leaf_names())): n.dist for n in t.traverse() if not n.is_root}
    assert edges(result) == edges(native)
    assert (out / "lsd2.outgroups.txt").read_text() == "2\nO1\nO2\n"
    # The substitution output is deliberately not the published dating tree.
    su = read_tree(out / "lsd2.result.nwk")
    assert su.get_distance("O1", "A") < 1


def test_only_interval_calibrations_report_nonunique_scale(tmp_path, lsd2):
    tree, qc, cal = inputs(tmp_path, bounds=(90, 110))
    out = tmp_path / "dated"
    date(tree, qc, cal, out, lsd2)
    report = json.loads((out / "provenance.json").read_text())
    assert report["root_age_ma"] == pytest.approx(100, abs=.001)
    assert report["root_age_solution_bounds_ma"] == [90, 110]
    assert report["unique_time_scale"] is False
    assert report["solution_bounds_are_confidence_intervals"] is False
    assert report["review_status"] == "needs_review"


@pytest.mark.parametrize("text", [
    "((O1:0,O2:0):.05,(A:.04,(B:.02,C:.02):0):.03);",
    "((O1:.02,O2:.02):.03,(A:.04,B:.04,C:.04):.03);",
    "(O1:.01,(O2:.02,(A:.04,(B:.02,C:.02):.02):.06):.01);",
])
def test_zero_edges_polytomies_and_anchor_root_are_supported(tmp_path, lsd2, text):
    tree, qc, cal = inputs(tmp_path, text)
    date(tree, qc, cal, tmp_path / "out", lsd2)
    original = read_tree(tree)
    result = read_tree(tmp_path / "out/species_tree.dated.nwk")
    assert len(list(original.traverse())) == len(list(result.traverse()))
    validate_root(result, ["O1", "O2"])


def test_nonmonophyletic_basal_group_rejected_without_running(tmp_path):
    tree, qc, cal = inputs(tmp_path, "((O1:.02,A:.02):.03,(O2:.04,(B:.02,C:.02):.02):.03);")
    with pytest.raises(ValueError, match="not monophyletic"):
        date(tree, qc, cal, tmp_path / "out", "must-not-run")


def test_rooting_preserves_pairwise_distances_and_all_tips():
    tree = Tree("(O1:.01,(O2:.02,(A:.04,(B:.02,C:.02):.02):.06):.01);", parser=1)
    names = list(tree.leaf_names())
    pairs = {(a, b): tree.get_distance(a, b) for a in names for b in names}
    root_on_outgroup(tree, ["O1", "O2"])
    validate_root(tree, ["O1", "O2"])
    assert set(tree.leaf_names()) == set(names)
    assert all(tree.get_distance(a, b) == pytest.approx(d) for (a, b), d in pairs.items())


def test_explicit_multiple_basal_species_expand_all_samples(tmp_path):
    names = ["O1_run1", "O1_run2", "O2_run1", "A", "B", "C"]
    samples, meta, out, qc = [tmp_path / p for p in ("samples", "meta", "out", "qc")]
    write_tsv(samples, ["species", "species_id", "taxid", "cds"], [
        {"species": n, "species_id": n.split("_")[0], "taxid": i, "cds": "unused"} for i, n in enumerate(names)])
    prepare_root(samples, meta, out, qc, ["O1", "O2"])
    assert out.read_text().splitlines() == names[:3]
    assert json.loads(qc.read_text())["outgroup_species"] == ["O1", "O2"]


def test_real_astral_to_lsd2_multiple_basal_samples(tmp_path, lsd2):
    binary = os.environ.get("ASTRAL_BIN")
    if not binary:
        pytest.skip("set ASTRAL_BIN for the real ASTRAL-IV/LSD2 connection")
    trees, merge_qc, manifest = [tmp_path / n for n in ("genes.nwk", "merge.json", "manifest.tsv")]
    trees.write_text("((O1:.02,O2:.02):.03,(A:.04,(B:.02,C:.02):.02):.03);\n" * 12)
    write_json(merge_qc, {"status": "retained", "mean_gene_length": 250, "retained": [{"sites": 250}] * 12})
    write_tsv(manifest, ["species"], [{"species": n} for n in ["O1", "O2", "A", "B", "C"]])
    group = tmp_path / "outgroup.txt"
    group.write_text("O1\nO2\n")
    tree, qc = tmp_path / "species_tree.nwk", tmp_path / "species_tree.json"
    astral(trees, merge_qc, manifest, tree, qc, binary, None, 1, 17, outgroup_file=group)
    validate_root(read_tree(tree), ["O1", "O2"])
    assert json.loads(qc.read_text())["astral_root_anchor"] == "O1"
    cal = tmp_path / "cal.tsv"
    write_tsv(cal, ["taxa", "min_age_ma", "max_age_ma", "source"], [
        {"taxa": "O1,A", "min_age_ma": 100, "max_age_ma": 100, "source": "synthetic"}])
    date(tree, qc, cal, tmp_path / "dated", lsd2)
    validate_root(read_tree(tmp_path / "dated/species_tree.dated.nwk"), ["O1", "O2"])


def test_only_floating_point_negative_residuals_are_clamped():
    text = "((A:-5e-15,B:0):100,(C:50,D:50):50);"
    normalized, records = normalize_numerical_zeros(text, 100)
    assert Tree(normalized, parser=1)["A"].dist == 0
    assert len(records) == 1 and records[0]["native_age_length_ma"] == -5e-15
    with pytest.raises(ValueError, match="beyond floating-point"):
        normalize_numerical_zeros(text.replace("-5e-15", "-0.0001"), 100)


def test_root_split_input_ratio_does_not_fix_dated_split(tmp_path, lsd2):
    outputs = []
    for index, (left, right) in enumerate([(.03, .03), (.054, .006)]):
        work = tmp_path / str(index)
        work.mkdir()
        tree, qc, cal = inputs(work, f"((O1:.02,O2:.02):{left},(A:.04,(B:.02,C:.02):.02):{right});")
        date(tree, qc, cal, work / "dated", lsd2)
        outputs.append((work / "dated/species_tree.dated.nwk").read_text())
    assert outputs[0] == outputs[1]
