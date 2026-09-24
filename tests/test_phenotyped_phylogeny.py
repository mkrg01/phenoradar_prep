"""Observed-trait selection and independent, reusable BUSCO inference runs."""
import json
import os
from pathlib import Path
import shutil
import sqlite3
import subprocess
import sys

import pytest
import yaml

from common import read_tsv, write_tsv
from infer_phylogeny import read_tree
from species_traits import select_phenotyped

ROOT = Path(__file__).resolve().parents[1]


def test_selection_keeps_zero_and_all_observed_species_after_busco_filtering(tmp_path):
    samples, traits = tmp_path / "samples.tsv", tmp_path / "traits.tsv"
    rows = [{"species": f"Plant_{n}", "run": n, "cds": "unused", "C4": "1"} for n in "ABCDEFGH"]
    rows.append(dict(rows[0], run="replicate"))
    write_tsv(samples, list(rows[0]), rows)
    write_tsv(traits, ["species", "C4"], [
        {"species": f"Plant {n}", "C4": v}
        for n, v in zip("ABCDEFGX", ["0", "1", "0", "1", "", "NA", "NaN", "1"])])
    out = tmp_path / "selection"
    select_phenotyped(samples, traits, out)
    assert read_tsv(out / "samples.tsv") == [r for r in rows if r["species"] in {f"Plant_{n}" for n in "ABCD"}]
    report = json.loads((out / "selection.json").read_text())
    assert report["dataset_species"] == 8 and report["inference_species"] == 4
    assert report["missing_trait_species"] == [f"Plant_{n}" for n in "EFGH"]
    assert report["missing_trait_rows"] == ["Plant_H"]
    assert report["trait"] == "C4" and report["species_trait"]["sha256"]


@pytest.mark.parametrize("trait,values", [("C4", ["0"] * 4), ("height", ["0", "1.2", "3", "4.5"])])
def test_selection_has_no_two_state_constraint_and_counts_species_not_runs(tmp_path, trait, values):
    samples, traits = tmp_path / "samples.tsv", tmp_path / "traits.tsv"
    rows = [{"species": n, "cds": "unused"} for n in "ABCD"]
    write_tsv(samples, list(rows[0]), rows * 2)
    write_tsv(traits, ["species", trait], [{"species": n, trait: v} for n, v in zip("ABCD", values)])
    select_phenotyped(samples, traits, tmp_path / "out", trait=trait)
    assert len(read_tsv(tmp_path / "out/samples.tsv")) == 8
    with pytest.raises(ValueError, match="at least 5 phenotyped species; found 4"):
        select_phenotyped(samples, traits, tmp_path / "too_small", trait=trait, min_taxa=5)
    with pytest.raises(ValueError, match="requires species and absent"):
        select_phenotyped(samples, traits, tmp_path / "invalid", trait="absent")


@pytest.mark.parametrize("sets", ["phenotyped", ["all", "all"], ["unknown"], [["all"]]])
def test_invalid_trees_fail_at_configuration(tmp_path, workflow_project, sets):
    snakemake = os.environ.get("SNAKEMAKE_BIN") or shutil.which("snakemake")
    if not snakemake:
        pytest.skip("Snakemake required")
    config = tmp_path / "override.yaml"
    config.write_text(yaml.safe_dump({"phylogeny": {"trees": sets}}))
    result = subprocess.run([snakemake, "--snakefile", str(ROOT / "workflow/Snakefile"),
                             "--configfile", str(config), "--list-rules"],
                            cwd=workflow_project, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            env={**os.environ, "XDG_CACHE_HOME": str(tmp_path / "cache")})
    assert result.returncode != 0 and "phylogeny.trees must be a list" in result.stdout


def test_species_set_switching_preserves_inference_and_dating(tmp_path, workflow_project, command_environment, seed_taxonomy):
    from test_phylogeny import phylogeny_inputs, trimal_binary
    nwkit = pytest.importorskip("nwkit")
    if nwkit.__version__ != "0.27.0":
        pytest.skip("requires workflow/envs/timetree.yaml")
    snakemake = os.environ.get("SNAKEMAKE_BIN") or shutil.which("snakemake")
    famsa = os.environ.get("FAMSA_BIN") or shutil.which("famsa")
    vft = os.environ.get("VERYFASTTREE_BIN") or shutil.which("VeryFastTree")
    if not all([snakemake, famsa, vft]) or not (ROOT / "resources/phylogeny_tools/bin/astral4_int128").exists():
        pytest.skip("real inference tools required")
    treepl = os.environ.get("TREEPL_BIN") or shutil.which("treePL")
    commands = {"python": sys.executable, "famsa": famsa, "trimal": trimal_binary(), "VeryFastTree": vft}
    if treepl:
        commands["treePL"] = treepl
    env = command_environment(commands)
    source, species = phylogeny_inputs(tmp_path)
    # Removing the unobserved basal species must select a different outgroup.
    with sqlite3.connect(source / "taxa.sqlite") as db:
        db.execute("INSERT INTO species VALUES (3,2,'Ingroup','','no rank','3,2,1')")
        db.execute("INSERT INTO species VALUES (4,3,'Inner clade','','no rank','4,3,2,1')")
        for i in range(1, len(species)):
            parent, track = (3, "3,2,1") if i == 1 else (4, "4,3,2,1")
            db.execute("UPDATE species SET parent=?, track=? WHERE taxid=?", (parent, f"{42+i},{track}", 42+i))
    traits = source / "species_trait.tsv"
    trait_rows = [{"species": n, "C4": "" if i == 0 else "0"} for i, n in enumerate(species)]
    write_tsv(traits, ["species", "C4"], trait_rows)
    seed_taxonomy(source / "taxa.sqlite")
    cfg = {"run_name": "test",
        "trait": "C4",
        "phylogeny": {"outgroup": "auto", "max_markers": 3}}
    config = tmp_path / "override.yaml"
    argv = [snakemake, "--snakefile", str(ROOT / "workflow/Snakefile"), "--configfile", str(config),
            "--cores", "2", "--resources", "mem_mb=8000",
            "--set-threads", "align_busco_marker=1", "infer_busco_gene_tree=1", "infer_busco_species_tree=2",
            "--set-resources", "infer_busco_species_tree:mem_mb=4000"]
    if os.environ.get("PHYLOGENY_CONDA_PREFIX"):
        argv += ["--use-conda", "--conda-prefix", os.environ["PHYLOGENY_CONDA_PREFIX"]]

    def run(sets, target="phylogeny"):
        cfg["phylogeny"]["trees"] = sets
        config.write_text(yaml.safe_dump(cfg))
        result = subprocess.run(argv + ["--", target], cwd=workflow_project, env=env, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        assert result.returncode == 0, result.stdout + "\n" + "\n".join(
            p.read_text()[-3000:] for p in (tmp_path / "logs").rglob("*.log"))
        return result.stdout

    def timestamps(*folders):
        return {p: p.stat().st_mtime_ns for folder in folders for p in folder.rglob("*") if p.is_file()}

    result = tmp_path / "results/test"
    full, observed = result / "phylogeny/all", result / "phylogeny/phenotyped"
    run(["phenotyped"], "phylogeny_prepare")
    assert not full.exists() and not (observed / "species_tree.nwk").exists()
    assert {r["species"] for r in read_tsv(observed / "plan/species.tsv")} == set(species[1:])
    assert (observed / "rooting/outgroup.txt").read_text().strip() == species[1]
    assert all(int(r["species"]) == 5 for r in read_tsv(observed / "plan/marker_stats.tsv"))
    run(["phenotyped"])
    read_tree(observed / "species_tree.nwk", species[1:])
    observed_times = timestamps(observed)
    assert not full.exists() and not (result / "phylogeny/all/rooting").exists()
    # The all-species branch must not even require the phenotype input file.
    unavailable = traits.with_suffix(".unavailable")
    traits.rename(unavailable)
    try:
        run(["all"])
    finally:
        unavailable.rename(traits)
    read_tree(full / "species_tree.nwk", species)
    assert (full / "rooting/outgroup.txt").read_text().strip() == species[0]
    assert timestamps(observed) == observed_times
    inference_times = timestamps(full, observed)
    run(["all", "phenotyped"])
    assert timestamps(full, observed) == inference_times
    assert "Nothing to be done" in run(["all", "phenotyped"])
    run(["phenotyped"])
    assert timestamps(full, observed) == inference_times

    # Retrieve separate calibrations from recorded responses, without HTTP.
    from timetree_calibrations import cached_response
    from test_timetree_calibrations import cache_missing_estimates, fake_fetch, payload
    cache = workflow_project / "resources/timetree_cache"
    for ids in [range(42, 48), range(43, 48)]:
        cached_response(ids, cache, delay=0,
                        backend=(fake_fetch(payload(ids, mrca=str(list(ids)))), {"synthetic": True}))
    for branch in [full, observed]:
        cache_missing_estimates(branch / "species_tree.nwk", dict(zip(species, range(42, 48))), cache)
    cfg["phylogeny"]["dating"] = {"enabled": True, "calibration_source": "timetree"}
    run(["phenotyped"], "phylogeny_calibrations")
    assert not (full / "timetree").exists()
    assert json.loads((observed / "timetree/provenance.json").read_text())["status"] == "ready"
    calibration_times = timestamps(observed / "timetree")
    run(["all", "phenotyped"], "phylogeny_calibrations")
    assert timestamps(observed / "timetree") == calibration_times
    assert {r["species"] for r in read_tsv(full / "timetree/taxa.tsv")} == set(species)
    assert {r["species"] for r in read_tsv(observed / "timetree/taxa.tsv")} == set(species[1:])

    if treepl:
        run(["phenotyped"], "timetree")
        read_tree(observed / "dating/species_tree.dated.nwk", species[1:])
        assert not (full / "dating").exists()
        dated_times = timestamps(observed / "dating")
        # The phylogeny target stops at inference even with dating enabled.
        run(["all", "phenotyped"])
        assert not (full / "dating").exists()
        run(["all", "phenotyped"], "timetree")
        read_tree(full / "dating/species_tree.dated.nwk", species)
        assert timestamps(observed / "dating") == dated_times
        assert "Nothing to be done" in run(["all", "phenotyped"])
    assert all(p.stat().st_mtime_ns == stamp for p, stamp in inference_times.items())

    # Changed membership rebuilds only the observed inference, keeping all output.
    full_times = timestamps(full)
    cfg["phylogeny"]["dating"]["enabled"] = False
    trait_rows[-1]["C4"] = "NA"
    write_tsv(traits, ["species", "C4"], trait_rows)
    run(["all", "phenotyped"])
    read_tree(observed / "species_tree.nwk", species[1:-1])
    assert timestamps(full) == full_times
    assert "Nothing to be done" in run(["all", "phenotyped"])
    assert not (result / "phylogeny/representatives").exists() and not (result / "orthogroups/expression").exists()
