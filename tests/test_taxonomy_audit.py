"""Real MonoPhy integration: native roles, missing ranks and species-tree reuse."""
import json
import os
from pathlib import Path
import shutil
import sqlite3
import subprocess
import sys

import pytest
import yaml

from common import read_tsv, sha256, write_json, write_tsv
from taxonomy_audit import Taxonomy, audit, parse_tree, validate_settings

ROOT = Path(__file__).resolve().parents[1]
NAMES = ["A_one", "A_query", "A_three", "A_two", "B_five", "B_four", "B_one", "B_three", "B_two", "Outgroup"]
WRONG = "(Outgroup:1,(((A_one:1,A_two:1):1,A_three:1):1,((B_one:1,A_query:1):1,((B_two:1,B_three:1):1,(B_four:1,B_five:1):1):1):1):1);"
RIGHT = "(Outgroup:1,(((A_one:1,A_two:1):1,(A_three:1,A_query:1):1):1,((B_one:1,B_two:1):1,(B_three:1,(B_four:1,B_five:1):1):1):1):1);"


def require_monophy():
    if not shutil.which("Rscript"):
        pytest.skip("Rscript with real MonoPhy 1.3.2 required")


@pytest.fixture
def audit_inputs(tmp_path):
    folder = tmp_path / "source"
    folder.mkdir(parents=True)
    taxonomy = folder / "taxa.sqlite"
    with sqlite3.connect(taxonomy) as db:
        db.executescript("""
            CREATE TABLE species(taxid INTEGER PRIMARY KEY, spname TEXT, rank TEXT, track TEXT);
            CREATE TABLE merged(taxid_old INTEGER, taxid_new INTEGER);
            INSERT INTO species VALUES (1,'root','no rank','1');
            INSERT INTO species VALUES (5,'OrderX','order','5,1');
            INSERT INTO species VALUES (10,'FamilyA','family','10,5,1');
            INSERT INTO species VALUES (20,'FamilyB','family','20,5,1');
            INSERT INTO species VALUES (30,'FamilyO','family','30,5,1');
            INSERT INTO merged VALUES (999,101);
        """)
        for i, name in enumerate(NAMES):
            family = 10 if name.startswith("A_") else 20 if name.startswith("B_") else 30
            db.execute("INSERT INTO species VALUES (?,?,?,?)", (100+i, name, "species", f"{100+i},{family},5,1"))
    samples = folder / "samples.tsv"
    rows = [dict(species=name, scientific_name=name.replace("_", " "), taxid=100+i, run=f"SRR{i+1}") for i, name in enumerate(NAMES)]
    rows.append({**rows[1], "run": "SRR99"})
    write_tsv(samples, list(rows[0]), rows)
    tree = folder / "species_tree.nwk"
    tree.write_text(WRONG + "\n")
    # Its recorded inference input deliberately does not exist: no gene-tree dependency.
    write_json(folder / "species_tree.json", {"species": len(NAMES), "outgroup": "Outgroup",
               "input": {"path": "/absent/gene_trees.nwk", "sha256": "0" * 64}})
    return dict(tree=tree, tree_qc=folder/"species_tree.json", samples=samples, taxonomy=taxonomy)


def run_audit(inputs, out, **settings):
    require_monophy()
    return audit(**inputs, outdir=out, settings={"ranks": ["family"], **settings})


def test_native_roles_run_linkage_and_legacy_report_replacement(audit_inputs, tmp_path):
    before = {k: sha256(v) for k, v in audit_inputs.items()}
    out = tmp_path / "audit"
    out.mkdir()
    write_json(out / "summary.json", {"report_type": "taxonomy_audit", "schema_version": 1})
    (out / "gene_evidence.tsv").write_text("obsolete detector output")
    (out / "contexts").mkdir()
    (out / "contexts/obsolete.svg").write_text("obsolete plot")
    result = run_audit(audit_inputs, out, ranks=["family", "subfamily", "species", "order"])
    assert result["method"] == "MonoPhy" and result["engine"]["version"] == "1.3.2"
    assert result["candidate_species"] == 1
    candidates = read_tsv(out / "candidates.tsv")
    assert {(r["species"], r["role"], r["focal_taxon"]) for r in candidates} == {
        ("A_query", "outlier", "FamilyA"), ("A_query", "intruder", "FamilyB")}
    assert all(r["run_ids"] == "SRR2;SRR99" for r in candidates)
    assert "gene_trees" not in result and not (out / "gene_evidence.tsv").exists()
    assert not (out / "contexts").exists()
    flagged = [r["run"] for r in read_tsv(out / "samples.tsv") if r["status"] == "review_flag"]
    assert flagged == ["SRR2", "SRR99"]
    rank_states = read_tsv(out / "rank_status.tsv")
    assert {r["status"] for r in rank_states if r["rank"] == "subfamily"} == {"missing_taxonomy_rank"}
    assert {r["status"] for r in rank_states if r["rank"] == "species"} == {"monotypic"}
    # Keep both roles, although native TipStates gives Intruder priority.
    assert next(r for r in rank_states if r["rank"] == "family" and r["species"] == "A_query")["outlier"] == "TRUE"
    assert (out / "ranks/family/monophy.rds").stat().st_size > 0
    for rank in ["family", "subfamily", "species", "order"]:
        assert (out / f"ranks/{rank}/tree.pdf").stat().st_size > 1000
        assert (out / f"ranks/{rank}/tree.svg").stat().st_size > 1000
    assert {k: sha256(v) for k, v in audit_inputs.items()} == before
    assert set(result["inputs"]) == {"tree", "tree_qc", "samples", "taxonomy"}
    # Verify the stored native object agrees with a fresh direct MonoPhy call.
    code = '''suppressPackageStartupMessages(library(MonoPhy)); a<-commandArgs(TRUE);
    t<-ape::read.tree(a[1]); d<-read.delim(a[2],check.names=FALSE);
    native<-AssessMonophyly(t,taxonomy=d[,c("tip","family")],outliercheck=TRUE,outlierlevel=0.5);
    stopifnot(identical(native,readRDS(a[3])))'''
    subprocess.run([shutil.which("Rscript"), "--vanilla", "-e", code, str(audit_inputs["tree"]),
                    str(out / "taxonomy.tsv"), str(out / "ranks/family/monophy.rds")], check=True)
    # Changing settings affects native roles, not an old reference-count filter.
    run_audit(audit_inputs, out, outlierlevel=0.4)
    assert len(read_tsv(out / "candidates.tsv")) > 2
    audit_inputs["tree"].write_text(RIGHT + "\n")
    result = run_audit(audit_inputs, out)
    assert result["candidate_species"] == 0 and read_tsv(out / "candidates.tsv") == []
    assert not (out / "ranks/subfamily").exists()


def test_missing_rank_is_not_an_intruder_and_all_flags_stay_in_plot(audit_inputs, tmp_path):
    with sqlite3.connect(audit_inputs["taxonomy"]) as db:
        db.execute("UPDATE species SET track='101,5,1' WHERE taxid=101")
    out = tmp_path / "audit"
    summary = run_audit(audit_inputs, out)
    assert summary["candidate_species"] == 0
    status = next(r for r in read_tsv(out / "rank_status.tsv") if r["species"] == "A_query")
    assert status["status"] == "missing_taxonomy_rank"
    tree = parse_tree((out / "ranks/family/assessment_tree.nwk").read_text(), NAMES)
    assert "A_query" not in set(tree.leaf_names())
    plots = read_tsv(out / "plot_members.tsv")
    assert len(plots) == len(NAMES) and next(r for r in plots if r["species"] == "A_query")["color"] == "#777777"


def test_grouped_intruders_are_detected_by_monophy(audit_inputs, tmp_path):
    # Two A tips travel together; this was excluded by the retired single-intruder rule.
    text = WRONG.replace("(A_one:1,A_two:1):1,A_three:1", "A_one:1,A_two:1")
    text = text.replace("A_query:1", "(A_query:1,A_three:1):1")
    audit_inputs["tree"].write_text(text)
    out = tmp_path / "audit"
    run_audit(audit_inputs, out)
    b_intruders = {r["species"] for r in read_tsv(out / "candidates.tsv")
                   if r["focal_taxon"] == "FamilyB" and r["role"] == "intruder"}
    assert b_intruders == {"A_query", "A_three"}
    shown = read_tsv(out / "plot_members.tsv")
    assert all(r["species"] == r["display_tip"] for r in shown if r["symbol"] != "none")


def test_rank_with_one_annotated_tip_is_reported_without_native_assessment(audit_inputs, tmp_path):
    with sqlite3.connect(audit_inputs["taxonomy"]) as db:
        db.execute("INSERT INTO species VALUES (11,'SubfamilyA','subfamily','11,10,5,1')")
        db.execute("UPDATE species SET track='101,11,10,5,1' WHERE taxid=101")
    out = tmp_path / "audit"
    summary = run_audit(audit_inputs, out, ranks=["subfamily"])
    assert summary["candidate_species"] == 0
    groups = read_tsv(out / "taxon_results.tsv")
    assert len(groups) == 1 and groups[0]["monophyly"] == "Monotypic"
    assert (out / "ranks/subfamily/not_assessed.txt").exists()
    assert not (out / "ranks/subfamily/monophy.rds").exists()
    assert {r["status"] for r in read_tsv(out / "samples.tsv")} == {"not_assessable"}


def test_failed_validation_preserves_report(audit_inputs, tmp_path):
    out = tmp_path / "audit"
    run_audit(audit_inputs, out)
    before = {p.relative_to(out): sha256(p) for p in out.rglob("*") if p.is_file()}
    report = json.loads(audit_inputs["tree_qc"].read_text())
    report["outgroup"] = "A_one"
    write_json(audit_inputs["tree_qc"], report)
    with pytest.raises(ValueError, match="root does not match"):
        run_audit(audit_inputs, out)
    assert {p.relative_to(out): sha256(p) for p in out.rglob("*") if p.is_file()} == before
    with pytest.raises(ValueError, match="overlaps"):
        audit(**audit_inputs, outdir=audit_inputs["tree"].parent)


def test_manifest_conflicts_fail(audit_inputs, tmp_path):
    rows = read_tsv(audit_inputs["samples"])
    rows[-1]["taxid"] = "9999"
    write_tsv(audit_inputs["samples"], list(rows[0]), rows)
    with pytest.raises(ValueError, match="conflicting species"):
        audit(**audit_inputs, outdir=tmp_path / "audit")


def test_frozen_taxonomy_merged_ids(audit_inputs):
    taxonomy = Taxonomy(audit_inputs["taxonomy"])
    try:
        assert taxonomy.lineage(999)["family"] == (10, "FamilyA")
        assert taxonomy.lineage(999999) == {}
    finally:
        taxonomy.db.close()


@pytest.mark.parametrize("text", [WRONG.replace("A_two", "A_one"), WRONG.replace("A_two", "Unknown"),
                                   WRONG + WRONG, WRONG.replace("A_two:1", "A_two:-1")])
def test_invalid_tree_inputs(text):
    with pytest.raises(ValueError):
        parse_tree(text, NAMES, exact=True)


@pytest.mark.parametrize("settings", [{"min_reference_species": 2}, {"max_plot_species": 24},
    {"ranks": ["family", "family"]}, {"enabled": "false"}, {"outlierlevel": 1.1},
    {"outlierlevel": True}, {"outlierlevel": float("nan")}, {"ranks": ["no rank"]},
    {"collapse_monophyletic": "true"}, {"exclude_species": ["A_query"]}, None, []])
def test_invalid_or_obsolete_settings(settings):
    with pytest.raises(ValueError):
        validate_settings(settings)


def test_snakemake_reuses_species_tree_without_gene_inputs(audit_inputs, workflow_project, command_environment):
    require_monophy()
    snakemake = os.environ.get("SNAKEMAKE_BIN") or shutil.which("snakemake")
    if not snakemake:
        pytest.skip("Snakemake required")
    project = workflow_project
    metadata = project / "results/test/metadata"
    metadata.mkdir(parents=True)
    shutil.copyfile(audit_inputs["samples"], metadata / "samples.tsv")
    taxdir = project / "resources/taxonomy"
    taxdir.mkdir(parents=True)
    shutil.copyfile(audit_inputs["taxonomy"], taxdir / "taxa.sqlite")
    originals = [metadata / "samples.tsv", taxdir / "taxa.sqlite"]
    (project / "results/test/run.json").write_text("{}\n")
    for branch in ["phylogeny/all", "phylogeny/phenotyped"]:
        folder = project / "results/test" / branch
        folder.mkdir(parents=True)
        for key in ["tree", "tree_qc"]:
            dest = folder / audit_inputs[key].name
            shutil.copyfile(audit_inputs[key], dest)
            originals.append(dest)
        (folder / "species_coverage.tsv").write_text("species\tgene_trees\n")
        if branch == "phylogeny/phenotyped":
            (folder / "selection").mkdir()
            shutil.copyfile(audit_inputs["samples"], folder / "selection/samples.tsv")
    before = {p: (sha256(p), p.stat().st_mtime_ns) for p in originals}
    override = project / "override.yaml"
    cfg = {"analysis": "test", "phylogeny": {"species_sets": ["all"]}, "taxonomy_audit": {"ranks": ["family"]}}
    environment = command_environment({"python": sys.executable})
    wrapper = project / "workflow/AuditSnakefile"
    wrapper.parent.mkdir()
    for name in ["scripts", "rules", "envs"]:
        (wrapper.parent / name).symlink_to(ROOT / "workflow" / name, target_is_directory=True)
    wrapper.write_text("\n".join(line for line in (ROOT / "workflow/Snakefile").read_text().splitlines()
                                 if not line.startswith("include:") or '"rules/taxonomy_audit.smk"' in line) + "\n")
    argv = [snakemake, "--snakefile", str(wrapper), "--configfile", str(override), "--cores", "1", "--"]
    def run(target="taxonomy_audit"):
        override.write_text(yaml.safe_dump(cfg))
        process = subprocess.run(argv + [target], cwd=project, env=environment, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        assert process.returncode == 0, process.stdout + "\n" + "\n".join(p.read_text() for p in (project/"logs").rglob("*.log"))
        return process.stdout
    full = project / "results/test/phylogeny/all/taxonomy_audit"
    phenotyped = project / "results/test/phylogeny/phenotyped/taxonomy_audit"
    assert "Nothing to be done" in run("phylogeny")
    run()
    assert json.loads((full/"summary.json").read_text())["candidate_species"] == 1
    assert not phenotyped.exists()
    assert "Nothing to be done" in run()
    cfg["taxonomy_audit"]["outlierlevel"] = 0.4
    run()
    assert json.loads((full/"summary.json").read_text())["candidate_species"] > 1
    full_times = {p: p.stat().st_mtime_ns for p in full.rglob("*") if p.is_file()}
    cfg["phylogeny"]["species_sets"] = ["phenotyped"]
    run()
    assert phenotyped.is_dir()
    assert {p: p.stat().st_mtime_ns for p in full.rglob("*") if p.is_file()} == full_times
    assert {p: (sha256(p), p.stat().st_mtime_ns) for p in originals} == before
    assert not list((project / "results").rglob("gene_trees.nwk"))
    assert not (project / "results/test/phylogeny/representatives").exists()
