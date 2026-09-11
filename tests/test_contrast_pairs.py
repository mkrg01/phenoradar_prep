"""Root selection and contrast membership, using the pinned nwkit implementation."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest
import yaml

from common import read_tsv, write_tsv
from species_traits import read_species_traits
from phylogeny_root import ncbi_tree, prepare_root, resolve_outgroup
from contrast_pairs import prepare, summarize, skim
from plot_contrast_tree import plot

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(autouse=True)
def pinned_nwkit():
    nwkit = pytest.importorskip("nwkit")
    if nwkit.__version__ != "0.27.0":
        pytest.skip("requires workflow/envs/timetree.yaml")


def test_traits_normalize_names_and_reject_collisions(tmp_path):
    path = tmp_path / "traits.tsv"
    path.write_text("species\tC4\nPlant alpha\t0\nPlant_beta\t1\nPlant gamma\t\n")
    assert read_species_traits(path) == {"Plant_alpha": "0", "Plant_beta": "1", "Plant_gamma": ""}
    path.write_text("species\tC4\nPlant alpha\t0\nPlant_alpha\t0\n")
    with pytest.raises(ValueError, match="duplicate"):
        read_species_traits(path)
    path.write_text("species\tC4\nPlant alpha\t2\n")
    with pytest.raises(ValueError, match="C4 must"):
        read_species_traits(path)


def test_ncbi_guide_uses_only_frozen_sqlite(tmp_path, tiny_inputs, monkeypatch):
    import nwkit.util
    monkeypatch.setattr(nwkit.util, "_download_ete_taxdump", lambda *a: pytest.fail("unexpected download"))
    samples = tmp_path / "samples.tsv"
    write_tsv(samples, ["species", "taxid", "cds"],
              [{"species": n, "taxid": t, "cds": "unused"} for n, t in [("Alpha_plant", 42), ("Beta_sp-X", 43), ("Gamma_plant", 44)]])
    original = Path(tiny_inputs["taxonomy_db"]).read_bytes()
    ncbi_tree(samples, tiny_inputs["taxonomy_db"], tmp_path / "tree.nwk", tmp_path / "taxids.tsv")
    from ete4 import Tree
    assert set(Tree((tmp_path / "tree.nwk").read_text(), parser=9).leaf_names()) == {"Alpha_plant", "Beta_sp-X", "Gamma_plant"}
    assert Path(tiny_inputs["taxonomy_db"]).read_bytes() == original


def test_rooting_uses_ncbi_split_and_rejects_multi_species_side(tmp_path):
    from ete4 import Tree
    scores = dict.fromkeys("OABCD", 90)
    name, record = resolve_outgroup(Tree("(O,((A,B),(C,D)));", parser=9), scores,
                                   reference=lambda *a: pytest.fail("unnecessary reference lookup"))
    assert name == "O" and record["source"] == "ncbi"
    with pytest.raises(ValueError, match="multi-species"):
        resolve_outgroup(Tree("((A,B),(C,D));", parser=9), scores)


def test_rooting_reference_preserves_every_basal_lineage():
    from ete4 import Tree
    scores = {"Out_species": 50, "A_species": 90, "A_other": 80, "B_species": 85}
    tree = Tree("(Out_species,(A_species,A_other),B_species);", parser=9)
    calls = []
    def reference(names):
        calls.append(names)
        return Tree("(Out_species,(A_species,B_species));", parser=9), {"test": True}
    name, record = resolve_outgroup(tree, scores, reference)
    assert name == "Out_species" and record["source"] == "nwkit_apgiv"
    assert calls == [["A_species", "B_species", "Out_species"]]


@pytest.mark.parametrize("guide", ["(A_species,(O_species,B_species));", "(O_species,A_species);"])
def test_rooting_never_selects_partial_or_unresolved_basal_clade(guide):
    from ete4 import Tree
    tree = Tree("(O_species,(A_species,A_other),B_species);", parser=9)
    with pytest.raises(ValueError):
        resolve_outgroup(tree, {n: 90 for n in tree.leaf_names()},
                         reference=lambda names: (Tree(guide, parser=9), {}))


def contrast_inputs(tmp_path):
    source = tmp_path / "input"
    source.mkdir()
    species = ["O", "A", "A2", "B", "C", "D", "E", "U"]
    write_tsv(source / "samples.tsv", ["species", "cds", "taxid"],
              [{"species": n, "cds": "unused", "taxid": i + 1} for i, n in enumerate(species)])
    # Deliberately wrong metadata traits must never enter the grouping.
    write_tsv(source / "metadata.tsv", ["species", "busco_percent", "C4"],
              [{"species": n, "busco_percent": 90 if n != "A2" else 70, "C4": 1} for n in species])
    write_tsv(source / "traits.tsv", ["species", "C4"],
              [{"species": n, "C4": "" if n == "U" else int(n in {"B", "D"})} for n in species])
    (source / "ncbi.nwk").write_text("(O,(((A,A2),B),(C,D),E,U));\n")
    (source / "outgroup.txt").write_text("O\n")
    prepare(source / "samples.tsv", source / "metadata.tsv", source / "traits.tsv", source / "ncbi.nwk",
            tmp_path / "selection")
    return source


def test_root_is_selected_within_skim_representatives_without_adding_species(tmp_path):
    source = contrast_inputs(tmp_path)
    traits = read_tsv(source / "traits.tsv")
    next(r for r in traits if r["species"] == "O")["C4"] = ""
    write_tsv(source / "traits.tsv", ["species", "C4"], traits)
    (source / "ncbi.nwk").write_text("(O,(E,(((A,A2),B),(C,D),U)));\n")
    selection = tmp_path / "selection"
    prepare(source / "samples.tsv", source / "metadata.tsv", source / "traits.tsv", source / "ncbi.nwk", selection)
    original = (selection / "samples.tsv").read_bytes()
    selected = {r["species"] for r in read_tsv(selection / "samples.tsv")}
    assert selected == {r["leaf_name"] for r in read_tsv(selection / "ncbi_skim.sampled.tsv")}
    assert "O" not in selected
    prepare_root(selection / "samples.tsv", source / "metadata.tsv", source / "outgroup.txt", source / "root.json",
                 tree=selection / "ncbi_skim.nwk")
    assert (source / "outgroup.txt").read_text().strip() == "E"
    assert json.loads((source / "root.json").read_text())["candidate_species_count"] == len(selected)
    assert (selection / "samples.tsv").read_bytes() == original
    with pytest.raises(ValueError, match="absent from inference species"):
        prepare_root(selection / "samples.tsv", source / "metadata.tsv", source / "outgroup.txt", source / "root.json",
                     outgroup="O")


def test_two_skims_compose_all_species_membership_and_plot(tmp_path):
    source = contrast_inputs(tmp_path)
    selected = {r["species"] for r in read_tsv(tmp_path / "selection/samples.tsv")}
    assert selected == {"O", "A", "B", "C", "D", "E"}
    tree = source / "inferred.nwk"
    tree.write_text("(O:1,(((A:1,C:1)0.95:1,B:1)0.98:1,(D:1,E:1)0.97:1)1:1);\n")
    out = tmp_path / "out"
    summarize(tree, tmp_path / "selection", source / "outgroup.txt", out)
    pairs = read_tsv(out / "contrast_pairs.tsv")
    assert len(pairs) == 2
    meta = {r["species"]: r for r in read_tsv(out / "species_metadata.tsv")}
    assert {meta[n]["group"] for n in ["A", "A2", "C"]} == {meta["A"]["group"]}
    assert meta["A"]["n_species_in_group"] == "3"
    assert meta["A"]["contrast_pair_id"] == meta["B"]["contrast_pair_id"]
    assert meta["D"]["contrast_pair_id"] == meta["E"]["contrast_pair_id"]
    assert meta["A"]["contrast_pair_id"] != meta["D"]["contrast_pair_id"]
    assert meta["U"]["C4"] == meta["U"]["contrast_pair_id"] == ""
    assert meta["O"]["role"] == "outgroup" and meta["O"]["contrast_pair_id"] == ""
    assert meta["O"]["C4"] == "0" and meta["O"]["group"] != ""
    assert meta["O"]["is_representative"] == "1"
    plot(out / "summary_tree.nwk", out / "species_metadata.tsv", out / "summary.json", out)
    assert (out / "summary_tree.pdf").read_bytes().startswith(b"%PDF")
    assert "Pair" in (out / "summary_tree.svg").read_text()


def test_polytomy_is_not_arbitrarily_split_into_pairs(tmp_path):
    source = contrast_inputs(tmp_path)
    tree = source / "inferred.nwk"
    tree.write_text("(O:1,(A:1,B:1,C:1,D:1,E:1):1);\n")
    out = tmp_path / "out"
    summarize(tree, tmp_path / "selection", source / "outgroup.txt", out)
    assert read_tsv(out / "contrast_pairs.tsv") == []
    assert len(json.loads((out / "summary.json").read_text())["unresolved_contrastive_clades"]) == 1
    assert all(r["contrast_pair_id"] == "" for r in read_tsv(out / "species_metadata.tsv"))


def test_empty_contrastive_result_and_deterministic_ties(tmp_path):
    from ete4 import Tree
    rows = [{"leaf_name": n, "trait": "0", "busco_percent": 90} for n in "ABCD"]
    _, a = skim(Tree("((A,B),(C,D));", parser=9), rows, tmp_path / "a", 13, topology_only=True)
    _, b = skim(Tree("((A,B),(C,D));", parser=9), rows[::-1], tmp_path / "b", 13, topology_only=True)
    assert a == b
    _, c = skim(Tree("((A,B),(C,D));", parser=9), rows, tmp_path / "c", 13, contrastive=True)
    assert not c and (tmp_path / "c.nwk").read_text() == ""
    assert read_tsv(tmp_path / "c.sampled.tsv") == []


def test_both_workflow_branches_infer_with_automatic_root(tmp_path, command_environment, workflow_project):
    from test_phylogeny import phylogeny_inputs, trimal_binary
    snakemake = os.environ.get("SNAKEMAKE_BIN") or shutil.which("snakemake")
    famsa = os.environ.get("FAMSA_BIN") or shutil.which("famsa")
    vft = os.environ.get("VERYFASTTREE_BIN") or shutil.which("VeryFastTree")
    if not all([snakemake, famsa, vft]) or not (ROOT / "resources/phylogeny_tools/bin/astral4_int128").exists():
        pytest.skip("real inference tools required")
    source, species = phylogeny_inputs(tmp_path)
    import sqlite3
    with sqlite3.connect(source / "taxa.sqlite") as db:
        db.execute("INSERT INTO species VALUES (3,2,'Ingroup','','no rank','3,2,1')")
        db.execute("INSERT INTO species VALUES (4,3,'Inner clade','','no rank','4,3,2,1')")
        for i in range(1, len(species)):
            parent, track = (3, "3,2,1") if i == 1 else (4, "4,3,2,1")
            db.execute("UPDATE species SET parent=?, track=? WHERE taxid=?", (parent, f"{42+i},{track}",42+i))
    traits = source / "traits.tsv"
    write_tsv(traits, ["species", "C4"], [{"species": n.replace("_", " "), "C4": "" if i == 0 else i % 2}
                                         for i, n in enumerate(species)])
    cfg = {"analysis": "test", "inputs": {"metadata": str(source / "metadata.tsv"), "species_trait": str(traits),
           "busco": str(source / "busco.tsv"), "cds_dir": str(source / "cds"), "quant_dir": str(source / "quant")},
           "taxonomy": {"source": str(source / "taxa.sqlite")},
           "phylogeny": {"busco_full_dir": str(source / "busco"), "outgroup": "auto", "max_markers": 3,
                         "align_threads": 1, "tree_threads": 1, "astral_threads": 2, "astral_mem_gb": 4}}
    config = tmp_path / "config.yaml"
    config.write_text(yaml.safe_dump(cfg))
    env = command_environment({"python": sys.executable, "famsa": famsa, "trimal": trimal_binary(), "VeryFastTree": vft})
    argv = [snakemake, "--snakefile", str(ROOT / "workflow/Snakefile"), "--configfile", str(config),
            "--cores", "2", "--resources", "mem_mb=8000"]
    conda_prefix = os.environ.get("PHYLOGENY_CONDA_PREFIX")
    if conda_prefix:
        argv += ["--use-conda", "--conda-prefix", conda_prefix]
        env.pop("PYTHONPATH", None)
    def run(targets):
        result = subprocess.run(argv + ["--", *targets], cwd=workflow_project, env=env, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        assert result.returncode == 0, result.stdout + "\n" + "\n".join(
            p.read_text()[-3000:] for p in (tmp_path / "logs").rglob("*.log"))
        return result.stdout
    run(["contrast_pairs"])
    result = tmp_path / "results/test"
    assert not (result / "phylogeny/all/species_tree.nwk").exists()
    assert not (result / "phylogeny/all/gene_trees.nwk").exists()
    assert (result / "phylogeny/representatives/rooting/outgroup.txt").read_text().strip() == species[1]
    assert json.loads((result / "phylogeny/representatives/species_tree.json").read_text())["outgroup"] == species[1]
    selected = {r["species"] for r in read_tsv(result / "phylogeny/representatives/selection/samples.tsv")}
    assert selected == {r["leaf_name"] for r in read_tsv(result / "phylogeny/representatives/selection/ncbi_skim.sampled.tsv")}
    assert species[0] not in selected and species[1] in selected
    assert "Nothing to be done" in run(["contrast_pairs"])
    run(["phylogeny"])
    assert (result / "phylogeny/all/rooting/outgroup.txt").read_text().strip() == species[0]
    assert json.loads((result / "phylogeny/all/species_tree.json").read_text())["outgroup"] == species[0]
    full_tree = result / "phylogeny/all/species_tree.nwk"
    full_mtime = full_tree.stat().st_mtime_ns
    contrast_tree = result / "phylogeny/representatives/species_tree.nwk"
    contrast_mtime = contrast_tree.stat().st_mtime_ns
    # The post-inference target uses the completed full tree, with no NCBI
    # representative selection. Selecting additional runs retains old outputs.
    run(["phylogeny_contrast_pairs"])
    full_pairs = result / "phylogeny/all/contrast/contrast_pairs.tsv"
    full_pair_mtime = full_pairs.stat().st_mtime_ns
    assert {r["species"] for r in read_tsv(full_pairs.parent / "species_metadata.tsv")} == set(species)
    assert "Nothing to be done" in run(["phylogeny_contrast_pairs"])
    cfg["phylogeny"]["species_sets"] = ["phenotyped"]
    config.write_text(yaml.safe_dump(cfg))
    run(["phylogeny_contrast_pairs"])
    observed = result / "phylogeny/phenotyped"
    assert {r["species"] for r in read_tsv(observed / "contrast/species_metadata.tsv")} == set(species[1:])
    observed_mtime = (observed / "species_tree.nwk").stat().st_mtime_ns
    cfg["phylogeny"]["species_sets"] = ["all", "phenotyped"]
    config.write_text(yaml.safe_dump(cfg))
    assert "Nothing to be done" in run(["phylogeny_contrast_pairs"])
    assert full_pairs.stat().st_mtime_ns == full_pair_mtime
    (full_pairs.parent / "summary_tree.pdf").unlink()
    run(["phylogeny_contrast_pairs"])
    assert full_pairs.stat().st_mtime_ns == full_pair_mtime
    assert full_tree.stat().st_mtime_ns == full_mtime
    assert contrast_tree.stat().st_mtime_ns == contrast_mtime
    assert (observed / "species_tree.nwk").stat().st_mtime_ns == observed_mtime
    cfg["phylogeny"]["species_sets"] = ["all"]
    config.write_text(yaml.safe_dump(cfg))
    (result / "phylogeny/representatives/contrast/summary_tree.pdf").unlink()
    run(["contrast_pairs"])
    assert contrast_tree.stat().st_mtime_ns == contrast_mtime
    assert full_tree.stat().st_mtime_ns == full_mtime
    # Change trait labels while keeping the representative set: rebuild contrast
    # membership, and leave the full inference completely unchanged.
    rows = read_tsv(traits)
    rows[1]["C4"] = "0"
    write_tsv(traits, ["species", "C4"], rows)
    run(["contrast_pairs", "phylogeny"])
    assert full_tree.stat().st_mtime_ns == full_mtime
    assert not (result / "orthogroups/expression").exists() and not (result / "orthogroups/mapping").exists()
