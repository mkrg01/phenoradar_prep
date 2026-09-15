"""Automatically collect completed results without upstream jobs."""
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest
import yaml

from common import file_record, read_tsv, sha256, write_json, write_tsv
from phenoradar_inputs import discover, export


ROOT = Path(__file__).resolve().parents[1]
SPECIES = ["Plant_A", "Plant_B-x", "Plant_C"]
META_FIELDS = ["species", "C4", "contrast_pair_id", "family"]
EXPR_FIELDS = ["species", "run", "orthogroup", "tpm"]
KO_FIELDS = ["species", "run", "ko", "tpm_sum"]


@pytest.fixture
def snapshot(tmp_path):
    source = tmp_path / "completed"
    samples = [dict(species=species, run=f"R{i}", odb_species=species.replace("-", "_"),
                    scientific_name=species.replace("_", " "), taxid=str(10 + i))
               for i, species in enumerate(SPECIES)]
    write_tsv(source / "metadata/samples.tsv", list(samples[0]), samples)
    write_tsv(source / "metadata/species_metadata.tsv", META_FIELDS,
              [dict(species=species, C4=trait, contrast_pair_id="", family="Plantaceae")
               for species, trait in zip(SPECIES, ["0", "1", ""])])
    write_tsv(source / "orthogroups/expression/tpm.tsv", EXPR_FIELDS,
              [dict(species=r["species"], run=r["run"], orthogroup="OG1", tpm="0.000001")
               for r in samples])
    write_tsv(source / "orthogroups/expression/mapping_qc.tsv", ["species", "run", "mapped_genes"],
              [dict(species=r["species"], run=r["run"], mapped_genes="1") for r in samples])
    write_tsv(source / "kegg/ko_tpm_sum.tsv", KO_FIELDS,
              [dict(species=r["species"], run=r["run"], ko="K00001", tpm_sum="123.45000")
               for r in samples])
    write_tsv(source / "kegg/mapping_qc.tsv", ["species", "run", "quantified_genes"],
              [dict(species=r["species"], run=r["run"], quantified_genes="1") for r in samples])
    write_tsv(source / "kegg/ko_support.tsv",
              ["species", "run", "ko", "annotated_genes", "quantified_genes", "tpm_sum"],
              [dict(species=r["species"], run=r["run"], ko="K00001", annotated_genes="1",
                    quantified_genes="1", tpm_sum="123.45000") for r in samples])
    # A KO can belong to several modules; exporting must preserve that overlap.
    write_tsv(source / "kegg/ko_modules.tsv", ["ko", "module"],
              [dict(ko="K00001", module="M00001"), dict(ko="K00001", module="M00002")])
    write_tsv(source / "kegg/ko_pathways.tsv", ["ko", "pathway"],
              [dict(ko="K00001", pathway="map00010")])
    write_json(source / "kegg/reference_qc.json", {
        "results": [file_record(source / f"kegg/ko_{name}s.tsv") for name in ["module", "pathway"]]})
    alignment = source / "orthogroups/alignments/OG1.faa"
    alignment.parent.mkdir()
    alignment.write_text("".join(f">{species}_g1 original description\naAC--D\n" for species in SPECIES))
    write_json(source / "orthogroups/alignments/provenance.json",
               {"alignments": [dict(orthogroup="OG1", alignment=file_record(alignment))]})
    return source


def state(folder):
    return {str(path.relative_to(folder)): (sha256(path), path.stat().st_mtime_ns)
            for path in folder.rglob("*") if path.is_file()}


def test_default_export_collects_available_results_and_preserves_source_values(snapshot, tmp_path):
    source = snapshot
    original = state(source)
    out = tmp_path / "inputs"
    export(source, out)
    assert state(source) == original
    assert read_tsv(out / "species_metadata.tsv") == read_tsv(source / "metadata/species_metadata.tsv")
    for relative, original_path in [
        ("tpm.tsv", "orthogroups/expression/tpm.tsv"), ("kegg/ko_tpm_sum.tsv", "kegg/ko_tpm_sum.tsv"),
        ("kegg/ko_modules.tsv", "kegg/ko_modules.tsv"), ("alignments/OG1.faa", "orthogroups/alignments/OG1.faa"),
    ]:
        assert (out / relative).is_symlink()
        assert (out / relative).resolve() == (source / original_path).resolve()
        assert (out / relative).read_bytes() == (source / original_path).read_bytes()
    files = {str(path.relative_to(out)) for path in out.rglob("*") if path.is_file()}
    assert files == {"species_metadata.tsv", "metadata/species_metadata.tsv", "metadata/samples.tsv",
                     "tpm.tsv", "orthogroups/expression/tpm.tsv", "orthogroups/expression/mapping_qc.tsv",
                     "kegg/ko_tpm_sum.tsv", "kegg/ko_support.tsv", "kegg/mapping_qc.tsv",
                     "kegg/ko_modules.tsv", "kegg/ko_pathways.tsv", "kegg/reference_qc.json",
                     "alignments/OG1.faa", "orthogroups/alignments/OG1.faa",
                     "orthogroups/alignments/provenance.json"}
    assert not (out / "manifest.json").exists()
    assert not (out / "alignments/members.tsv").exists()


def test_kegg_inputs_need_no_orthogroup_or_alignment_results(snapshot, tmp_path):
    shutil.rmtree(snapshot / "orthogroups/expression")
    shutil.rmtree(snapshot / "orthogroups/alignments")
    out = tmp_path / "inputs"
    export(snapshot, out)
    assert (out / "kegg/ko_tpm_sum.tsv").is_file()
    assert (out / "kegg/ko_modules.tsv").is_file()
    assert (out / "kegg/ko_pathways.tsv").is_file()
    assert not (out / "tpm.tsv").exists()
    assert not (out / "alignments").exists()
    assert not (snapshot / "orthogroups/mapping").exists()


def test_missing_group_map_does_not_hide_other_maps_or_expression(snapshot, tmp_path):
    (snapshot / "kegg/ko_modules.tsv").unlink()
    out = tmp_path / "inputs"
    export(snapshot, out)
    assert (out / "kegg/ko_tpm_sum.tsv").is_file()
    assert not (out / "kegg/ko_modules.tsv").exists()
    assert (out / "kegg/ko_pathways.tsv").is_file()


def test_discovery_uses_only_existing_absolute_dependencies(snapshot):
    inventory = discover(snapshot)
    assert all(Path(path).is_absolute() and Path(path).is_file() for path in inventory["files"])
    assert all(inventory["sections"][key] == "ready" for key in ["orthogroups", "kegg", "alignments"])
    assert not any("/orthogroups/mapping/" in path or "/proteins/" in path for path in inventory["files"])


@pytest.mark.parametrize("missing", ["kegg/ko_support.tsv", "kegg/mapping_qc.tsv"])
def test_incomplete_expression_does_not_hide_completed_reference_maps(snapshot, tmp_path, missing):
    (snapshot / missing).unlink()
    out = tmp_path / "inputs"
    export(snapshot, out)
    assert (out / "tpm.tsv").is_file()
    assert not (out / "kegg/ko_tpm_sum.tsv").exists()
    assert (out / "kegg/ko_modules.tsv").is_file()


def test_unverified_group_maps_do_not_hide_completed_expression(snapshot, tmp_path):
    (snapshot / "kegg/reference_qc.json").unlink()
    out = tmp_path / "inputs"
    export(snapshot, out)
    assert (out / "kegg/ko_tpm_sum.tsv").is_file()
    assert not (out / "kegg/ko_modules.tsv").exists()


def test_incomplete_alignment_is_skipped(snapshot, tmp_path):
    (snapshot / "orthogroups/alignments/provenance.json").unlink()
    out = tmp_path / "inputs"
    export(snapshot, out)
    assert not (out / "alignments").exists()
    assert discover(snapshot)["sections"]["alignments"] == "incomplete"


def test_collection_does_not_require_expression(snapshot, tmp_path):
    shutil.rmtree(snapshot / "orthogroups/expression")
    shutil.rmtree(snapshot / "kegg")
    out = tmp_path / "inputs"
    export(snapshot, out)
    assert (out / "species_metadata.tsv").is_file()
    assert (out / "alignments/OG1.faa").is_file()
    assert not (out / "tpm.tsv").exists()


@pytest.mark.parametrize("branch,fields,value_column", [
    ("orthogroups/expression/tpm.tsv", EXPR_FIELDS, "tpm"), ("kegg/ko_tpm_sum.tsv", KO_FIELDS, "tpm_sum"),
])
@pytest.mark.parametrize("bad_value", ["", "NA", "NaN", "inf", "-0.1", "oops"])
def test_expression_requires_finite_nonnegative_numbers(snapshot, tmp_path, branch, fields, value_column, bad_value):
    rows = read_tsv(snapshot / branch)
    rows[0][value_column] = bad_value
    write_tsv(snapshot / branch, fields, rows)
    with pytest.raises(ValueError):
        export(snapshot, tmp_path / "inputs")


@pytest.mark.parametrize("branch,fields", [("orthogroups/expression/tpm.tsv", EXPR_FIELDS), ("kegg/ko_tpm_sum.tsv", KO_FIELDS)])
@pytest.mark.parametrize("problem", ["duplicate", "missing_species", "wrong_run", "unknown_species"])
def test_expression_coordinates_must_match_selected_samples(snapshot, tmp_path, branch, fields, problem):
    rows = read_tsv(snapshot / branch)
    if problem == "duplicate":
        rows.append(dict(rows[0]))
    elif problem == "missing_species":
        rows.pop()
    elif problem == "wrong_run":
        rows[0]["run"] = rows[1]["run"]
    else:
        rows[0]["species"] = "Unknown_plant"
    write_tsv(snapshot / branch, fields, rows)
    with pytest.raises(ValueError):
        export(snapshot, tmp_path / "inputs")


def test_multiple_runs_per_species_are_collected_without_aggregation(snapshot, tmp_path):
    for relative in ["metadata/samples.tsv", "orthogroups/expression/tpm.tsv",
                     "orthogroups/expression/mapping_qc.tsv", "kegg/ko_tpm_sum.tsv",
                     "kegg/mapping_qc.tsv", "kegg/ko_support.tsv"]:
        path = snapshot / relative
        rows = read_tsv(path)
        rows.append({**rows[0], "run": "Replicate"})
        write_tsv(path, list(rows[0]), rows)
    out = tmp_path / "inputs"
    export(snapshot, out)
    assert read_tsv(out / "tpm.tsv") == read_tsv(snapshot / "orthogroups/expression/tpm.tsv")
    assert len(read_tsv(out / "species_metadata.tsv")) == 3
    assert len(read_tsv(out / "metadata/samples.tsv")) == 4


def test_metadata_must_cover_selected_species_without_duplicates(snapshot, tmp_path):
    path = snapshot / "metadata/species_metadata.tsv"
    rows = read_tsv(path)
    write_tsv(path, META_FIELDS, rows[:-1] + [rows[0]])
    with pytest.raises(ValueError):
        export(snapshot, tmp_path / "inputs")


def test_all_contrast_branches_keep_independent_pair_ids(snapshot, tmp_path):
    for name in ["all", "phenotyped", "representatives"]:
        branch = snapshot / f"phylogeny/{name}/contrast"
        write_tsv(branch / "species_metadata.tsv", ["species", "C4", "contrast_pair_id"],
                  [dict(species=species, C4=trait, contrast_pair_id=f"{name}_pair")
                   for species, trait in zip(SPECIES[:2], ["0", "1"])])
        write_json(branch / "summary.json", {"pairs": 1})
        write_tsv(branch / "contrast_pairs.tsv", ["species1", "species2"],
                  [dict(species1=SPECIES[0], species2=SPECIES[1])])
    out = tmp_path / "inputs"
    export(snapshot, out)
    assert read_tsv(out / "species_metadata.tsv") == read_tsv(snapshot / "metadata/species_metadata.tsv")
    for name in ["all", "phenotyped", "representatives"]:
        branch = out / f"phylogeny/{name}/contrast"
        assert {r["contrast_pair_id"] for r in read_tsv(branch / "species_metadata.tsv")} == {f"{name}_pair"}
        assert (branch / "contrast_pairs.tsv").is_symlink()


def test_incomplete_contrast_is_skipped(snapshot, tmp_path):
    write_tsv(snapshot / "phylogeny/representatives/contrast/species_metadata.tsv", META_FIELDS,
              read_tsv(snapshot / "metadata/species_metadata.tsv"))
    out = tmp_path / "inputs"
    export(snapshot, out)
    assert not (out / "phylogeny/representatives/contrast").exists()


def completed_tree(snapshot, newick, branch="all"):
    tree = snapshot / f"phylogeny/{branch}/species_tree.nwk"
    tree.parent.mkdir(parents=True, exist_ok=True)
    tree.write_text(newick)
    write_json(tree.with_suffix(".json"), {"species": 3})
    return tree


def test_completed_newick_and_available_gzipped_annotations_are_linked(snapshot, tmp_path):
    tree = completed_tree(snapshot, "((Plant_A,Plant_B-x),Plant_C);\n")
    annotations = snapshot / "orthogroup_annotations.tsv.gz"
    with gzip.open(annotations, "wt") as handle:
        handle.write("OG1\t3193\tA useful annotation\n")
    out = tmp_path / "inputs"
    export(snapshot, out)
    assert (out / "phylogeny/all/species_tree.nwk").resolve() == tree
    assert (out / "orthogroup_annotations.tsv.gz").resolve() == annotations
    assert not (out / "orthogroup_annotations.tsv").exists()


@pytest.mark.parametrize("newick", [
    "(Plant_A,Plant_B-x);", "((Plant_A,Plant_B-x),Unknown_plant);",
    "((Plant_A,Plant_B-x),(Plant_C,Extra_plant));", "((Plant_A,Plant_B-x),(Plant_C,Plant_C));",
])
def test_full_species_tree_requires_exact_unique_species_tips(snapshot, tmp_path, newick):
    completed_tree(snapshot, newick)
    with pytest.raises(ValueError):
        export(snapshot, tmp_path / "inputs")


def test_alignment_completion_checksum_is_verified(snapshot, tmp_path):
    path = snapshot / "orthogroups/alignments/OG1.faa"
    path.write_text(path.read_text().replace("aAC", "GAC"))
    with pytest.raises(ValueError):
        export(snapshot, tmp_path / "inputs")


@pytest.mark.parametrize("group,old,new", [("module", "M00001", "M00003"),
                                          ("pathway", "map00010", "map00020")])
def test_wellformed_kegg_map_changes_fail_completion_checksum(snapshot, tmp_path, group, old, new):
    path = snapshot / f"kegg/ko_{group}s.tsv"
    path.write_text(path.read_text().replace(old, new))
    with pytest.raises(ValueError, match="checksum"):
        export(snapshot, tmp_path / "inputs")


def test_refresh_removes_stale_links_and_failure_preserves_previous_export(snapshot, tmp_path):
    out = tmp_path / "inputs"
    export(snapshot, out)
    shutil.rmtree(snapshot / "kegg")
    shutil.rmtree(snapshot / "orthogroups/alignments")
    export(snapshot, out)
    assert (out / "tpm.tsv").is_file()
    assert not (out / "kegg").exists()
    assert not (out / "alignments").exists()
    before = state(out)
    completed_tree(snapshot, "(Unknown_plant,Plant_A);")
    with pytest.raises(ValueError):
        export(snapshot, out)
    assert state(out) == before


def make_filtered(source, exclusions):
    filtered = source / "filtered"
    keep = set(SPECIES) - set(exclusions)
    for relative in ["metadata/samples.tsv", "metadata/species_metadata.tsv", "orthogroups/expression/tpm.tsv", "orthogroups/expression/mapping_qc.tsv"]:
        rows = read_tsv(source / relative)
        write_tsv(filtered / relative, list(rows[0]), [row for row in rows if row["species"] in keep])
    records = [{**file_record(path), "path": str(path.relative_to(filtered)), "storage": "file"}
               for path in sorted(filtered.rglob("*")) if path.is_file()]
    report = dict(report_type="species_filter", schema_version=2, source=str(source.resolve()),
                  exclude_species=exclusions, retained_species=sorted(keep),
                  retained_runs=[f"R{i}" for i, species in enumerate(SPECIES) if species in keep], outputs=records)
    write_json(filtered / "manifest.json", report)
    return filtered


def test_exclusions_select_only_matching_completed_filtered_snapshot(snapshot, tmp_path):
    filtered = make_filtered(snapshot, ["Plant_C"])
    out = tmp_path / "inputs"
    export(snapshot, out, exclusions=["Plant_C"])
    assert (out / "tpm.tsv").resolve() == (filtered / "orthogroups/expression/tpm.tsv").resolve()
    assert {row["species"] for row in read_tsv(out / "species_metadata.tsv")} == set(SPECIES[:2])
    # Clearing exclusions returns to the original completed dataset.
    export(snapshot, out, exclusions=[])
    assert (out / "tpm.tsv").resolve() == (snapshot / "orthogroups/expression/tpm.tsv").resolve()
    assert {row["species"] for row in read_tsv(out / "species_metadata.tsv")} == set(SPECIES)


@pytest.mark.parametrize("problem", ["missing", "exclusions", "source", "retained_species",
                                     "retained_runs", "original_run", "report_type", "hash"])
def test_stale_or_invalid_filtered_snapshots_are_rejected(snapshot, tmp_path, problem):
    filtered = make_filtered(snapshot, ["Plant_C"])
    path = filtered / "manifest.json"
    report = json.loads(path.read_text())
    if problem == "missing":
        path.unlink()
    elif problem == "hash":
        table = filtered / "orthogroups/expression/tpm.tsv"
        table.write_text(table.read_text().replace("0.000001", "42"))
    elif problem == "original_run":
        samples = snapshot / "metadata/samples.tsv"
        rows = read_tsv(samples)
        rows[0]["run"] = "New_R0"
        write_tsv(samples, list(rows[0]), rows)
    else:
        field = "exclude_species" if problem == "exclusions" else problem
        report[field] = {"exclusions": ["Plant_A"], "source": str(tmp_path / "elsewhere"),
                         "retained_species": SPECIES, "retained_runs": ["R0", "Unknown_run"],
                         "report_type": "other"}[problem]
        write_json(path, report)
    with pytest.raises(ValueError):
        export(snapshot, tmp_path / "inputs", exclusions=["Plant_C"])


def test_full_snakefile_collects_frozen_results_without_upstream_inputs(snapshot, workflow_project, command_environment):
    snakemake = os.environ.get("SNAKEMAKE_BIN") or shutil.which("snakemake")
    if not snakemake:
        pytest.skip("Snakemake required")
    project = workflow_project
    target = project / "results/test"
    target.parent.mkdir()
    shutil.copytree(snapshot, target)
    before = state(target)
    cfg = project / "override.yaml"
    cfg.write_text(yaml.safe_dump({"run_name": "test"}))
    environment = command_environment({"python": sys.executable})
    argv = [snakemake, "--snakefile", str(ROOT / "workflow/Snakefile"), "--configfile", str(cfg),
            "--cores", "1", "--", "phenoradar_inputs"]

    def run(success=True):
        process = subprocess.run(argv, cwd=project, env=environment, text=True,
                                 stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        logs = "\n".join(path.read_text() for path in (project / "logs").rglob("*.log"))
        assert (process.returncode == 0) == success, process.stdout + "\n" + logs
        return process.stdout

    run()
    out = target / "phenoradar_inputs"

    def published_state():
        links = {str(path.relative_to(out)): (os.readlink(path), path.lstat().st_ino)
                 for path in out.rglob("*") if path.is_symlink()}
        return out.stat().st_ino, state(out), links

    assert (out / "tpm.tsv").is_symlink()
    assert (out / "kegg/ko_modules.tsv").is_file()
    # This explicit target always validates its completed inputs before refresh.
    assert "Nothing to be done" not in run()
    changed = target / "orthogroups/expression/tpm.tsv"
    changed.write_text(changed.read_text().replace("0.000001", "0.000002"))
    assert "Nothing to be done" not in run()
    assert read_tsv(out / "tpm.tsv")[0]["tpm"] == "0.000002"
    # Snakemake must not remove the previous publication before validation or
    # while cleaning a failed job. Source edits naturally affect existing links.
    valid_expression = changed.read_text()
    changed.write_text(valid_expression.replace("0.000002", "NaN"))
    previous = published_state()
    run(success=False)
    assert published_state() == previous
    changed.write_text(valid_expression)
    before["orthogroups/expression/tpm.tsv"] = (sha256(changed), changed.stat().st_mtime_ns)
    user_file = out / "keep.txt"
    user_file.write_text("user-owned file must survive a rejected refresh\n")
    previous = published_state()
    run(success=False)
    assert published_state() == previous
    assert user_file.read_text() == "user-owned file must survive a rejected refresh\n"
    user_file.unlink()
    # Newly completed outputs are discovered; removed outputs leave no stale links.
    completed_tree(target, "((Plant_A,Plant_B-x),Plant_C);")
    run()
    assert (out / "phylogeny/all/species_tree.nwk").is_symlink()
    for relative in ["kegg", "orthogroups/alignments"]:
        shutil.rmtree(target / relative)
    before = {name: value for name, value in state(target).items() if not name.startswith("phenoradar_inputs/")}
    run()
    assert not (out / "kegg").exists()
    assert not (out / "alignments").exists()
    after = {name: value for name, value in state(target).items() if not name.startswith("phenoradar_inputs/")}
    assert after == before
    assert not (project / "work").exists()


def test_collects_all_tree_variants_and_future_inputs(snapshot, tmp_path):
    expected = set()
    for branch, newick in [("all", "((Plant_A,Plant_B-x),Plant_C);"),
                           ("phenotyped", "(Plant_A,Plant_B-x);"),
                           ("representatives", "(Plant_A,Plant_C);")]:
        tree = completed_tree(snapshot, newick, branch)
        expected.add(str(tree.relative_to(snapshot)))
        base = f"phylogeny/{branch}"
        write_json(snapshot / base / "gene_trees.json", {"retained": [{"marker": "marker1"}]})
        for relative, contents in {
            "gene_trees.nwk": newick + "\n",
            "alignments/marker1.faa": ">Plant_A\nMK--\n",
            "alignments/raw/marker1.faa": ">Plant_A\nMK--\n",
            "markers/marker1.faa": ">Plant_A\nMK\n",
            "gene_trees/marker1.nwk": "(Plant_A,Plant_C);\n",
            "dating/species_tree.dated.nwk": newick + "\n",
            "dating/node_ages.tsv": "node\tage_ma\nroot\t100\n",
            "timetree/calibrations.tsv": "taxa\tmin_age_ma\tmax_age_ma\tsource\n",
            "rooting/outgroup.txt": "Plant_A\n",
            "taxonomy_check/ranks/family.nwk": newick + "\n",
            "taxonomy_check/samples.tsv": "species\tstatus\nPlant_A\tno_flag\n",
        }.items():
            path = snapshot / base / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(contents)
            expected.add(f"{base}/{relative}")
        for directory, marker in [("dating", "provenance.json"), ("timetree", "provenance.json"),
                                  ("rooting", "outgroup.json"), ("taxonomy_check", "summary.json")]:
            write_json(snapshot / base / directory / marker, {"completed": True})
        # Native solver work and stale loci are not completed downstream results.
        work = snapshot / base / "dating/treepl_runs/attempt/input.nwk"
        work.parent.mkdir(parents=True)
        work.write_text(newick)
        (snapshot / base / "alignments/stale.faa").write_text(">Unknown_species\nAA\n")
    mapping = snapshot / "orthogroups/mapping"
    mapping.mkdir()
    (mapping / "mappings.sqlite").write_bytes(b"completed mapping database")
    write_tsv(mapping / "gene_orthogroups.tsv", ["gene_id", "orthogroup"],
              [dict(gene_id="Plant_A_g1", orthogroup="OG1")])
    write_json(mapping / "merge_qc.json", {"pairs": 1})
    expected.update(f"orthogroups/mapping/{name}" for name in
                    ["mappings.sqlite", "gene_orthogroups.tsv", "merge_qc.json"])
    for name in SPECIES:
        stem = name.replace("-", "_") + "_protein"
        path = snapshot / "proteins" / f"{stem}.fa"
        path.parent.mkdir(exist_ok=True)
        path.write_text(f">{name}_g1\nMK\n")
        write_json(path.with_suffix(".json"), {"protein": file_record(path)})
        expected.add(str(path.relative_to(snapshot)))
    (snapshot / "proteins/Stale_protein.fa").write_text(">Stale_g1\nMK\n")
    before = state(snapshot)
    out = tmp_path / "inputs"
    export(snapshot, out)
    assert state(snapshot) == before
    for relative in expected:
        assert (out / relative).is_symlink(), relative
        assert (out / relative).resolve() == snapshot / relative
    assert not list(out.rglob("treepl_runs"))
    assert not list(out.rglob("stale.faa"))
    assert not (out / "proteins/Stale_protein.fa").exists()
    assert not (out / "species_tree.nwk").exists()  # No arbitrary choice among tree types.


def test_incomplete_tree_and_dating_are_skipped_independently(snapshot, tmp_path):
    tree = completed_tree(snapshot, "((Plant_A,Plant_B-x),Plant_C);")
    tree.with_suffix(".json").unlink()
    dated = tree.parent / "dating/species_tree.dated.nwk"
    dated.parent.mkdir()
    dated.write_text(tree.read_text())
    out = tmp_path / "inputs"
    export(snapshot, out)
    assert not (out / "phylogeny/all").exists()
    write_json(tree.with_suffix(".json"), {"species": 3})
    export(snapshot, out)
    assert (out / "phylogeny/all/species_tree.nwk").is_symlink()
    assert not (out / "phylogeny/all/dating").exists()


def test_unknown_files_inside_collected_directories_survive_failed_refresh(snapshot, tmp_path):
    out = tmp_path / "inputs"
    export(snapshot, out)
    user_file = out / "metadata/notes.txt"
    user_file.write_text("Do not delete my notes\n")
    before = state(out)
    with pytest.raises(ValueError, match="unrecognized file"):
        export(snapshot, out)
    assert state(out) == before


def test_completed_busco_sequences_and_alignments_collect_before_tree_inference(snapshot, tmp_path):
    base = snapshot / "phylogeny/all"
    write_json(base / "plan/provenance.json", {"complete": True})
    write_tsv(base / "plan/markers.tsv", ["marker"], [{"marker": "locus1"}, {"marker": "locus2"}])
    write_tsv(base / "plan/species.tsv", ["species"], [{"species": "Plant_A"}])
    for relative, text in [("species/Plant_A.faa", ">locus1\nMK\n"),
                           ("alignments/raw/locus1.faa", ">Plant_A\nMK-\n")]:
        path = base / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        write_json(path.with_suffix(".json"), {"completed": True})
    # An incomplete second alignment and a stale locus must stay outside the collection.
    (base / "alignments/raw/locus2.faa").write_text(">Plant_A\nMK\n")
    (base / "alignments/raw/stale.faa").write_text(">Plant_A\nMK\n")
    write_json(base / "alignments/raw/stale.json", {"completed": True})
    out = tmp_path / "inputs"
    export(snapshot, out)
    assert (out / "phylogeny/all/species/Plant_A.faa").is_symlink()
    assert (out / "phylogeny/all/alignments/raw/locus1.faa").is_symlink()
    assert not (out / "phylogeny/all/alignments/raw/locus2.faa").exists()
    assert not (out / "phylogeny/all/alignments/raw/stale.faa").exists()
    assert not (out / "phylogeny/all/species_tree.nwk").exists()
