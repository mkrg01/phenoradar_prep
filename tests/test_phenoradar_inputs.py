"""Publish only completed, coherent PhenoRadar inputs without upstream jobs."""
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
from phenoradar_inputs import discover, export, validate_settings


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


def expression_only(**settings):
    return dict(orthogroups=True, kegg=False, alignments=False, **settings)


@pytest.mark.parametrize("settings", [
    [], "auto", {"unknown": True}, {"orthogroups": "yes"}, {"kegg": 1},
    {"alignments": None}, {"kegg_groups": "module"}, {"kegg_groups": ["other"]},
    {"kegg_groups": ["module", "module"]}, {"contrast": "../contrast"},
])
def test_invalid_settings_are_rejected(settings):
    with pytest.raises(ValueError):
        validate_settings(settings)


def test_default_export_is_minimal_and_preserves_source_values(snapshot, tmp_path):
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
    assert files == {"species_metadata.tsv", "tpm.tsv", "kegg/ko_tpm_sum.tsv",
                     "kegg/ko_modules.tsv", "alignments/OG1.faa"}
    assert not (out / "manifest.json").exists()
    assert not (out / "alignments/members.tsv").exists()


def test_kegg_inputs_need_no_orthogroup_or_alignment_results(snapshot, tmp_path):
    shutil.rmtree(snapshot / "orthogroups/expression")
    shutil.rmtree(snapshot / "orthogroups/alignments")
    out = tmp_path / "inputs"
    export(snapshot, out, settings={"kegg": True, "kegg_groups": ["module", "pathway"]})
    assert (out / "kegg/ko_tpm_sum.tsv").is_file()
    assert (out / "kegg/ko_modules.tsv").is_file()
    assert (out / "kegg/ko_pathways.tsv").is_file()
    assert not (out / "tpm.tsv").exists()
    assert not (out / "alignments").exists()
    assert not (snapshot / "orthogroups/mapping").exists()


def test_kegg_group_maps_are_opt_in(snapshot, tmp_path):
    out = tmp_path / "inputs"
    export(snapshot, out, settings={"kegg_groups": []})
    assert (out / "kegg/ko_tpm_sum.tsv").is_file()
    assert not (out / "kegg/ko_modules.tsv").exists()
    assert not (out / "kegg/ko_pathways.tsv").exists()


def test_discovery_uses_only_existing_absolute_dependencies(snapshot):
    inventory = discover(snapshot)
    assert all(Path(path).is_absolute() and Path(path).is_file() for path in inventory["files"])
    assert inventory["sections"] == dict(orthogroups="ready", kegg="ready", alignments="ready")
    assert not any("/orthogroups/mapping/" in path or "/proteins/" in path for path in inventory["files"])


@pytest.mark.parametrize("missing", ["kegg/ko_support.tsv", "kegg/reference_qc.json",
                                     "kegg/ko_modules.tsv", "kegg/mapping_qc.tsv"])
def test_auto_skips_incomplete_kegg_but_explicit_request_fails(snapshot, tmp_path, missing):
    (snapshot / missing).unlink()
    out = tmp_path / "inputs"
    export(snapshot, out)
    assert (out / "tpm.tsv").is_file()
    assert not (out / "kegg").exists()
    before = state(out)
    with pytest.raises(ValueError):
        export(snapshot, out, settings={"kegg": True})
    assert state(out) == before


def test_auto_skips_incomplete_alignment_and_required_mode_reports_it(snapshot, tmp_path):
    (snapshot / "orthogroups/alignments/provenance.json").unlink()
    out = tmp_path / "inputs"
    export(snapshot, out)
    assert not (out / "alignments").exists()
    with pytest.raises(ValueError):
        export(snapshot, tmp_path / "required", settings={"alignments": True})


def test_export_requires_at_least_one_completed_expression_branch(snapshot, tmp_path):
    shutil.rmtree(snapshot / "orthogroups/expression")
    shutil.rmtree(snapshot / "kegg")
    with pytest.raises(ValueError):
        export(snapshot, tmp_path / "inputs")


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


def test_multiple_runs_per_species_are_not_silently_aggregated(snapshot, tmp_path):
    path = snapshot / "metadata/samples.tsv"
    rows = read_tsv(path)
    rows.append({**rows[0], "run": "Replicate"})
    write_tsv(path, list(rows[0]), rows)
    with pytest.raises(ValueError):
        export(snapshot, tmp_path / "inputs")


def test_metadata_must_cover_selected_species_without_duplicates(snapshot, tmp_path):
    path = snapshot / "metadata/species_metadata.tsv"
    rows = read_tsv(path)
    write_tsv(path, META_FIELDS, rows[:-1] + [rows[0]])
    with pytest.raises(ValueError):
        export(snapshot, tmp_path / "inputs")


def test_contrast_pairs_are_left_joined_onto_complete_metadata(snapshot, tmp_path):
    branch = snapshot / "phylogeny/phenotyped/contrast"
    write_tsv(branch / "species_metadata.tsv", ["species", "C4", "contrast_pair_id"],
              [dict(species=species, C4=trait, contrast_pair_id="pair_1")
               for species, trait in zip(SPECIES[:2], ["0", "1"])])
    write_json(branch / "summary.json", {"pairs": 1})
    out = tmp_path / "inputs"
    export(snapshot, out, settings={"contrast": "phylogeny/phenotyped/contrast"})
    rows = {row["species"]: row for row in read_tsv(out / "species_metadata.tsv")}
    assert set(rows) == set(SPECIES)
    assert rows["Plant_A"]["contrast_pair_id"] == "pair_1"
    assert rows["Plant_B-x"]["contrast_pair_id"] == "pair_1"
    assert rows["Plant_C"]["contrast_pair_id"] == ""
    assert rows["Plant_C"]["C4"] == ""
    assert all(row["family"] == "Plantaceae" for row in rows.values())


def test_selected_contrast_requires_completion(snapshot, tmp_path):
    write_tsv(snapshot / "phylogeny/representatives/contrast/species_metadata.tsv", META_FIELDS,
              read_tsv(snapshot / "metadata/species_metadata.tsv"))
    with pytest.raises(ValueError):
        export(snapshot, tmp_path / "inputs", settings={"contrast": "phylogeny/representatives/contrast"})


def test_external_newick_and_gzipped_annotations_are_linked(snapshot, tmp_path, monkeypatch):
    tree = tmp_path / "tree.nwk"
    tree.write_text("((Plant_A,Plant_B-x),Plant_C);\n")
    annotations = tmp_path / "ogs.tsv.gz"
    with gzip.open(annotations, "wt") as handle:
        handle.write("OG1\t3193\tA useful annotation\n")
    monkeypatch.chdir(tmp_path)
    out = tmp_path / "inputs"
    export(snapshot, out, settings={"tree": "tree.nwk", "orthogroup_annotations": "ogs.tsv.gz"})
    assert (out / "species_tree.nwk").resolve() == tree
    assert (out / "orthogroup_annotations.tsv.gz").resolve() == annotations
    assert not (out / "orthogroup_annotations.tsv").exists()


@pytest.mark.parametrize("newick", [
    "(Plant_A,Plant_B-x);", "((Plant_A,Plant_B-x),Unknown_plant);",
    "((Plant_A,Plant_B-x),(Plant_C,Extra_plant));", "((Plant_A,Plant_B-x),(Plant_C,Plant_C));",
])
def test_selected_tree_requires_exact_unique_species_tips(snapshot, tmp_path, newick):
    tree = tmp_path / "tree.nwk"
    tree.write_text(newick)
    with pytest.raises(ValueError):
        export(snapshot, tmp_path / "inputs", settings={"tree": str(tree)})


def test_alignment_completion_checksum_is_verified(snapshot, tmp_path):
    path = snapshot / "orthogroups/alignments/OG1.faa"
    path.write_text(path.read_text().replace("aAC", "GAC"))
    with pytest.raises(ValueError):
        export(snapshot, tmp_path / "inputs", settings={"alignments": True})


@pytest.mark.parametrize("group,old,new", [("module", "M00001", "M00003"),
                                          ("pathway", "map00010", "map00020")])
def test_wellformed_kegg_map_changes_fail_completion_checksum(snapshot, tmp_path, group, old, new):
    path = snapshot / f"kegg/ko_{group}s.tsv"
    path.write_text(path.read_text().replace(old, new))
    with pytest.raises(ValueError, match="checksum"):
        export(snapshot, tmp_path / "inputs", settings={"kegg_groups": [group]})


def test_refresh_removes_stale_optional_links_and_failure_preserves_previous_export(snapshot, tmp_path):
    out = tmp_path / "inputs"
    export(snapshot, out, settings={"kegg_groups": ["module", "pathway"]})
    export(snapshot, out, settings=expression_only())
    assert (out / "tpm.tsv").is_file()
    assert not (out / "kegg").exists()
    assert not (out / "alignments").exists()
    before = state(out)
    with pytest.raises(ValueError):
        export(snapshot, out, settings=expression_only(tree=str(tmp_path / "missing.nwk")))
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
    export(snapshot, out, settings=expression_only(), exclusions=["Plant_C"])
    assert (out / "tpm.tsv").resolve() == (filtered / "orthogroups/expression/tpm.tsv").resolve()
    assert {row["species"] for row in read_tsv(out / "species_metadata.tsv")} == set(SPECIES[:2])
    # Clearing exclusions returns to the original completed dataset.
    export(snapshot, out, settings=expression_only(), exclusions=[])
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
        export(snapshot, tmp_path / "inputs", settings=expression_only(), exclusions=["Plant_C"])


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
    cfg.write_text(yaml.safe_dump({"analysis": "test", "phenoradar": {"kegg_groups": ["module"]}}))
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
    # Changing selected branches removes stale links without requesting producers.
    cfg.write_text(yaml.safe_dump({"analysis": "test", "phenoradar": {"kegg": False, "alignments": False}}))
    run()
    assert not (out / "kegg").exists()
    assert not (out / "alignments").exists()
    after = {name: value for name, value in state(target).items() if not name.startswith("phenoradar_inputs/")}
    assert after == before
    assert not (project / "work").exists()
