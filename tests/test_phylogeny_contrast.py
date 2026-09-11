"""Pairs on completed molecular trees, including exclusion-only recomputation."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest
import yaml

from common import read_tsv, sha256, write_json, write_tsv
from contrast_pairs import from_tree
from filter_species import discover, export
from plot_contrast_tree import plot

ROOT = Path(__file__).resolve().parents[1]
SPECIES = ["O", "A", "A2", "B", "C", "D", "E", "U"]
KNOWN = SPECIES[1:-1]
TREE = "(O:1,(((A:1,A2:1)0.9:1,B:1)0.9:1,(C:1,(D:1,E:1)0.9:1)0.9:1,U:1)0.9:1);"


@pytest.fixture(autouse=True)
def pinned_nwkit():
    nwkit = pytest.importorskip("nwkit")
    if nwkit.__version__ != "0.27.0":
        pytest.skip("requires workflow/envs/timetree.yaml")


@pytest.fixture
def molecular_snapshot(tmp_path):
    source = tmp_path / "completed"
    rows = [dict(species=s, scientific_name=s, taxid=i+1, odb_species=s,
                 run=f"R{i}", cds="unavailable", busco_percent=70 if s == "A2" else 95)
            for i, s in enumerate(SPECIES)]
    rows.append({**rows[1], "run": "replicate_A"})
    write_tsv(source / "metadata/samples.tsv", list(rows[0]), rows)
    write_tsv(source / "metadata/metadata_high_busco.tsv", list(rows[0]), rows)
    traits = tmp_path / "traits.tsv"
    write_tsv(traits, ["species", "C4", "other"],
              [dict(species=s, C4="" if s in {"O", "U"} else int(s in {"B", "D"}),
                    other="" if s == "B" else int(s == "D")) for s in SPECIES])
    full = source / "phylogeny/all"
    full.mkdir(parents=True)
    (full / "species_tree.nwk").write_text(TREE + "\n")
    write_json(full / "species_tree.json", {"species": len(SPECIES), "outgroup": "O"})
    observed = source / "phylogeny/phenotyped"
    write_tsv(observed / "selection/samples.tsv", list(rows[0]), [r for r in rows if r["species"] in KNOWN])
    (observed / "species_tree.nwk").write_text("(A:1,((A2:1,B:1)0.9:1,(C:1,(D:1,E:1)0.9:1)0.9:1)0.9:1);\n")
    write_json(observed / "species_tree.json", {"species": len(KNOWN), "outgroup": "A"})
    # This legacy representative result must never be read, pruned or rerun.
    (source / "phylogeny/representatives/contrast").mkdir(parents=True)
    (source / "phylogeny/representatives/contrast/contrast_pairs.tsv").write_text("historical NCBI representative result\n")
    return source, traits


def inputs(source, traits, branch="phylogeny/all"):
    samples = source / ("metadata/samples.tsv" if branch == "phylogeny/all" else f"{branch}/selection/samples.tsv")
    return dict(tree=source / branch / "species_tree.nwk", tree_qc=source / branch / "species_tree.json",
                samples=samples, metadata=source / "metadata/metadata_high_busco.tsv", traits=traits)


def state(folder):
    return {str(p.relative_to(folder)): (sha256(p), p.stat().st_mtime_ns)
            for p in folder.rglob("*") if p.is_file()}


def pair_members(out):
    return {frozenset([r["representative_a"], r["representative_b"]])
            for r in read_tsv(out / "contrast_pairs.tsv")}


def test_direct_membership_uses_molecular_tree_and_inherits_root(molecular_snapshot, tmp_path):
    from ete4 import Tree
    source, traits = molecular_snapshot
    before = state(source)
    out = tmp_path / "pairs"
    report = from_tree(**inputs(source, traits), outdir=out)
    assert pair_members(out) == {frozenset(["A", "B"]), frozenset(["D", "E"])}
    meta = {r["species"]: r for r in read_tsv(out / "species_metadata.tsv")}
    assert set(meta) == set(SPECIES)
    assert meta["A2"]["representative"] == "A" and meta["A"]["n_species_in_group"] == "2"
    assert meta["A2"]["contrast_pair_id"] == meta["B"]["contrast_pair_id"]
    assert meta["C"]["contrast_pair_id"] == ""
    assert all(meta[n]["contrast_pair_id"] == meta[n]["group"] == "" for n in ["O", "U"])
    assert report["outgroup"] == "O" and report["outgroup_in_observed_tree"] is False
    assert report["rooting"] == "source_root_inherited" and report["reestimated"] is False
    old, new = Tree(TREE, parser=0), Tree((out / "observed_tree.nwk").read_text(), parser=0)
    assert set(new.leaf_names()) == set(KNOWN)
    for a in KNOWN:
        for b in KNOWN:
            assert old.get_distance(a, b) == pytest.approx(new.get_distance(a, b))
    assert state(source) == before


def test_exclusion_reforms_pairs_and_does_not_just_remove_rows(molecular_snapshot, tmp_path):
    source, traits = molecular_snapshot
    out = tmp_path / "pairs"
    from_tree(**inputs(source, traits), outdir=out, exclude_species=["E"])
    assert pair_members(out) == {frozenset(["A", "B"]), frozenset(["C", "D"])}
    assert "E" not in {r["species"] for r in read_tsv(out / "species_metadata.tsv")}
    from_tree(**inputs(source, traits), outdir=out)
    assert frozenset(["D", "E"]) in pair_members(out)


@pytest.mark.parametrize("keep", [{"A", "A2", "C", "E", "O", "U"}, {"A", "O"}, {"O", "U"}, set()])
def test_degenerate_subsets_clear_pairs_and_render(molecular_snapshot, tmp_path, keep):
    source, traits = molecular_snapshot
    out = tmp_path / "pairs"
    from_tree(**inputs(source, traits), outdir=out)
    report = from_tree(**inputs(source, traits), outdir=out, exclude_species=sorted(set(SPECIES) - keep))
    assert report["contrast_pairs"] == 0 and pair_members(out) == set()
    meta = read_tsv(out / "species_metadata.tsv")
    assert {r["species"] for r in meta} == keep
    assert all(r["contrast_pair_id"] == "" for r in meta)
    plot(out / "summary_tree.nwk", out / "species_metadata.tsv", out / "summary.json", out)
    assert (out / "summary_tree.pdf").read_bytes().startswith(b"%PDF")
    assert (out / "summary_tree.svg").is_file()


def test_phenotyped_scope_and_different_pair_trait(molecular_snapshot, tmp_path):
    source, traits = molecular_snapshot
    out = tmp_path / "pairs"
    report = from_tree(**inputs(source, traits, "phylogeny/phenotyped"), outdir=out, trait="other", exclude_species=["A", "O"])
    meta = {r["species"]: r for r in read_tsv(out / "species_metadata.tsv")}
    assert set(meta) == set(KNOWN) - {"A"}
    assert meta["B"]["group"] == meta["B"]["contrast_pair_id"] == ""
    assert report["excluded_species"] == ["A"] and report["outgroup_in_observed_tree"] is False


@pytest.mark.parametrize("damage", ["duplicate_tip", "missing_tip", "wrong_root", "wrong_qc", "missing_score", "conflicting_score"])
def test_invalid_source_is_rejected(molecular_snapshot, tmp_path, damage):
    source, traits = molecular_snapshot
    args = inputs(source, traits)
    if damage == "duplicate_tip":
        args["tree"].write_text(TREE.replace("A2:", "A:") + "\n")
    elif damage == "missing_tip":
        args["tree"].write_text(TREE.replace(",U:1", "") + "\n")
    elif damage == "wrong_root":
        write_json(args["tree_qc"], {"species": len(SPECIES), "outgroup": "A"})
    elif damage == "wrong_qc":
        write_json(args["tree_qc"], {"species": 99, "outgroup": "O"})
    else:
        rows = read_tsv(args["metadata"])
        if damage == "missing_score":
            rows = [r for r in rows if r["species"] != "A"]
        else:
            rows[-1]["busco_percent"] = "80"
        write_tsv(args["metadata"], list(rows[0]), rows)
    with pytest.raises(ValueError):
        from_tree(**args, outdir=tmp_path / "pairs")


def test_filter_recomputes_both_completed_trees_without_inference(molecular_snapshot, tmp_path):
    source, traits = molecular_snapshot
    before = state(source)
    out = tmp_path / "filtered"
    report = export(source, ["E"], out, traits)
    for branch in ["phylogeny/all", "phylogeny/phenotyped"]:
        assert report["sections"][f"{branch}_contrast"]["status"] == "ready"
        assert frozenset(["C", "D"]) in pair_members(out / branch / "contrast")
        assert not (out / branch / "species_tree.nwk").exists()
        assert not (out / branch / "gene_trees.nwk").exists()
    assert not (out / "phylogeny/representatives").exists()
    assert state(source) == before
    assert all(sha256(out / r["path"]) == r["sha256"] for r in report["outputs"])
    export(source, [], out, traits)
    assert frozenset(["D", "E"]) in pair_members(out / "phylogeny/all/contrast")
    export(source, KNOWN, out, traits)
    assert read_tsv(out / "phylogeny/phenotyped/contrast/species_metadata.tsv") == []
    assert read_tsv(out / "phylogeny/all/contrast/contrast_pairs.tsv") == []
    assert state(source) == before


def test_filter_inventory_tracks_new_trees_but_not_legacy_pairs(molecular_snapshot):
    source, traits = molecular_snapshot
    original = discover(source, traits)
    assert not any("/contrast/" in p for p in original["files"])
    (source / "phylogeny/representatives/contrast/contrast_pairs.tsv").write_text("changed history\n")
    assert discover(source, traits) == original
    (source / "phylogeny/phenotyped/species_tree.nwk").unlink()
    changed = discover(source, traits)
    assert changed != original and changed["sections"]["phylogeny/phenotyped_contrast"]["status"] == "incomplete"


def test_full_snakefile_exclusions_never_request_inference(molecular_snapshot, workflow_project, command_environment):
    snakemake = os.environ.get("SNAKEMAKE_BIN") or shutil.which("snakemake")
    if not snakemake:
        pytest.skip("Snakemake required")
    source, traits = molecular_snapshot
    project = workflow_project
    target = project / "results/test"
    target.parent.mkdir()
    shutil.copytree(source, target)
    before = state(target)
    cfg = project / "override.yaml"
    environment = command_environment({"python": sys.executable})
    argv = [snakemake, "--snakefile", str(ROOT / "workflow/Snakefile"), "--configfile", str(cfg), "--cores", "1", "--"]
    def run(excluded, target_name="phylogeny_contrast_pairs"):
        cfg.write_text(yaml.safe_dump({"analysis": "test", "exclude_species": excluded,
                                      "phylogeny": {"species_sets": ["all", "phenotyped"]},
                                      "inputs": {"species_trait": str(traits)}}))
        process = subprocess.run(argv + [target_name], cwd=project, env=environment, text=True,
                                 stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        assert process.returncode == 0, process.stdout + "\n" + "\n".join(p.read_text() for p in (project / "logs").rglob("*.log"))
        return process.stdout
    run(["E"])
    out = target / "filtered"
    assert frozenset(["C", "D"]) in pair_members(out / "phylogeny/all/contrast")
    assert "Nothing to be done" in run(["E"])
    run(["A"], "filter_species")
    assert frozenset(["A2", "B"]) in pair_members(out / "phylogeny/all/contrast")
    assert "Nothing to be done" in run(["A"])
    assert {p: v for p, v in state(target).items() if not p.startswith("filtered/")} == before
    assert not list(target.rglob("gene_trees.nwk"))
