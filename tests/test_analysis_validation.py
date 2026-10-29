"""Native output checks preserve data and never create an input bundle."""
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
from analysis_validation import discover, validate


ROOT = Path(__file__).resolve().parents[1]
SPECIES = ["Plant_A", "Plant_B-x", "Plant_C"]
META_FIELDS = ["species", "C4", "contrast_pair_id", "family"]
EXPR_FIELDS = ["species", "run", "orthogroup", "tpm"]
KO_FIELDS = ["species", "run", "ko", "tpm_sum"]
SEQ_FIELDS = ["species", "sequence_id", "tpm"]


@pytest.fixture
def snapshot(tmp_path):
    source = tmp_path / "completed"
    samples = [dict(species=species, run=f"R{i}", odb_species=species.replace("-", "_"),
                    scientific_name=species.replace("_", " "), taxid=str(10 + i))
               for i, species in enumerate(SPECIES)]
    write_tsv(source / "metadata/samples.tsv", list(samples[0]), samples)
    write_tsv(source / "expression/sequence_tpm.tsv.gz", SEQ_FIELDS,
              [dict(species=r["species"], sequence_id=r["species"] + "_g1", tpm="0.00000")
               for r in samples])
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
    # A KO can belong to several modules; validation must allow that overlap.
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


def test_discovery_uses_only_existing_absolute_dependencies(snapshot):
    inventory = discover(snapshot)
    assert all(Path(path).is_absolute() and Path(path).is_file() for path in inventory["files"])
    assert all(inventory["sections"][key] == "ready" for key in ["orthogroups", "kegg", "alignments", "sequence_expression"])
    assert not any("/orthogroups/mapping/" in path or "/proteins/" in path for path in inventory["files"])


@pytest.mark.parametrize("branch,fields,value_column", [
    ("orthogroups/expression/tpm.tsv", EXPR_FIELDS, "tpm"), ("kegg/ko_tpm_sum.tsv", KO_FIELDS, "tpm_sum"),
    ("expression/sequence_tpm.tsv.gz", SEQ_FIELDS, "tpm"),
])
@pytest.mark.parametrize("bad_value", ["", "NA", "NaN", "inf", "-0.1", "oops"])
def test_expression_requires_finite_nonnegative_numbers(snapshot, tmp_path, branch, fields, value_column, bad_value):
    rows = read_tsv(snapshot / branch)
    rows[0][value_column] = bad_value
    write_tsv(snapshot / branch, fields, rows)
    with pytest.raises(ValueError):
        validate(snapshot)


@pytest.mark.parametrize("branch,fields", [("orthogroups/expression/tpm.tsv", EXPR_FIELDS), ("kegg/ko_tpm_sum.tsv", KO_FIELDS),
                                         ("expression/sequence_tpm.tsv.gz", SEQ_FIELDS)])
@pytest.mark.parametrize("problem", ["duplicate", "missing_species", "wrong_run", "unknown_species"])
def test_expression_coordinates_must_match_selected_samples(snapshot, tmp_path, branch, fields, problem):
    rows = read_tsv(snapshot / branch)
    if problem == "duplicate":
        rows.append(dict(rows[0]))
    elif problem == "missing_species":
        rows.pop()
    elif problem == "wrong_run":
        if "run" in fields:
            rows[0]["run"] = rows[1]["run"]
        else:
            rows[0]["sequence_id"] = rows[1]["sequence_id"]
    else:
        rows[0]["species"] = "Unknown_plant"
    write_tsv(snapshot / branch, fields, rows)
    with pytest.raises(ValueError):
        validate(snapshot)


def test_metadata_must_cover_selected_species_without_duplicates(snapshot, tmp_path):
    path = snapshot / "metadata/species_metadata.tsv"
    rows = read_tsv(path)
    write_tsv(path, META_FIELDS, rows[:-1] + [rows[0]])
    with pytest.raises(ValueError):
        validate(snapshot)


def test_older_snapshot_without_sequence_expression_remains_valid(snapshot):
    (snapshot / "expression/sequence_tpm.tsv.gz").unlink()
    assert validate(snapshot)["sections"]["sequence_expression"] == "absent"


@pytest.mark.parametrize("layout", ["native", "mixed", "legacy"])
def test_expression_layouts_validate_without_modifying_inputs(snapshot, layout):
    for index, branch in enumerate(["orthogroups/expression/tpm.tsv", "kegg/ko_tpm_sum.tsv"]):
        path = snapshot / branch
        rows = read_tsv(path)
        if layout == "native" or (layout == "mixed" and index == 0):
            fields = [key for key in rows[0] if key != "run"]
            write_tsv(path, fields, [{key: row[key] for key in fields} for row in rows])
    before = state(snapshot)
    assert validate(snapshot)["samples"] == len(SPECIES)
    assert state(snapshot) == before


@pytest.mark.parametrize("problem", ["unknown_sample", "missing_sample", "duplicate", "discontiguous", "ambiguous_manifest"])
def test_native_expression_resolves_explicit_metadata_without_parsing_ids(snapshot, problem):
    from analysis_validation import validate_expression
    path = snapshot / "orthogroups/expression/tpm.tsv"
    fields = ["species", "orthogroup", "tpm"]
    rows = [{key: row[key] for key in fields} for row in read_tsv(path)]
    runs = {f"private_run_with_underscores_{i}": name for i, name in enumerate(SPECIES)}
    if problem == "unknown_sample":
        rows[0]["species"] = "Unknown_R0"
    elif problem == "missing_sample":
        rows.pop()
    elif problem in {"duplicate", "discontiguous"}:
        extra = dict(rows[0])
        if problem == "discontiguous":
            extra["orthogroup"] = "OG2"
        rows.append(extra)
    else:
        runs["extra_run"] = SPECIES[0]
    write_tsv(path, fields, rows)
    with pytest.raises(ValueError):
        validate_expression(path, runs, "orthogroup", "tpm")


def completed_tree(snapshot, newick, branch="all"):
    tree = snapshot / f"phylogeny/{branch}/species_tree.nwk"
    tree.parent.mkdir(parents=True, exist_ok=True)
    tree.write_text(newick)
    write_json(tree.with_suffix(".json"), {"species": 3})
    return tree


@pytest.mark.parametrize("newick", [
    "(Plant_A,Plant_B-x);", "((Plant_A,Plant_B-x),Unknown_plant);",
    "((Plant_A,Plant_B-x),(Plant_C,Extra_plant));", "((Plant_A,Plant_B-x),(Plant_C,Plant_C));",
])
def test_full_species_tree_requires_exact_unique_species_tips(snapshot, tmp_path, newick):
    completed_tree(snapshot, newick)
    with pytest.raises(ValueError):
        validate(snapshot)


def test_alignment_completion_checksum_is_verified(snapshot, tmp_path):
    path = snapshot / "orthogroups/alignments/OG1.faa"
    path.write_text(path.read_text().replace("aAC", "GAC"))
    with pytest.raises(ValueError):
        validate(snapshot)


@pytest.mark.parametrize("group,old,new", [("module", "M00001", "M00003"),
                                          ("pathway", "map00010", "map00020")])
def test_wellformed_kegg_map_changes_fail_completion_checksum(snapshot, tmp_path, group, old, new):
    path = snapshot / f"kegg/ko_{group}s.tsv"
    path.write_text(path.read_text().replace(old, new))
    with pytest.raises(ValueError, match="checksum"):
        validate(snapshot)


def make_filtered(source, exclusions):
    filtered = source / "filtered"
    keep = set(SPECIES) - set(exclusions)
    for relative in ["metadata/samples.tsv", "metadata/species_metadata.tsv", "orthogroups/expression/tpm.tsv", "orthogroups/expression/mapping_qc.tsv", "expression/sequence_tpm.tsv.gz"]:
        rows = read_tsv(source / relative)
        write_tsv(filtered / relative, list(rows[0]), [row for row in rows if row["species"] in keep])
    records = [{**file_record(path), "path": str(path.relative_to(filtered)), "storage": "file"}
               for path in sorted(filtered.rglob("*")) if path.is_file()]
    report = dict(report_type="species_filter", schema_version=2, source=str(source.resolve()),
                  exclude_species=exclusions, retained_species=sorted(keep),
                  retained_runs=[f"R{i}" for i, species in enumerate(SPECIES) if species in keep], outputs=records)
    write_json(filtered / "manifest.json", report)
    return filtered


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
        validate(snapshot, exclusions=["Plant_C"])


def test_validation_is_read_only_and_preserves_legacy_bundle(snapshot):
    legacy = snapshot / "phenoradar_inputs"
    legacy.mkdir()
    (legacy / "notes.txt").write_text("Existing user notes\n")
    (legacy / "tpm.tsv").write_text("Existing exported values\n")
    (legacy / "species_metadata.tsv").symlink_to(snapshot / "metadata/species_metadata.tsv")
    before = state(snapshot)
    report = validate(snapshot)
    assert report["source"] == str(snapshot)
    assert report["sections"]["orthogroups"] == "ready"
    assert not any("phenoradar_inputs" in name for name in report["files"])
    assert state(snapshot) == before
    assert list(read_tsv(snapshot / "orthogroups/expression/tpm.tsv")[0]) == EXPR_FIELDS
    assert (legacy / "species_metadata.tsv").is_symlink()


def test_verification_does_not_generate_optional_results(snapshot):
    shutil.rmtree(snapshot / "orthogroups/expression")
    shutil.rmtree(snapshot / "orthogroups/alignments")
    before = state(snapshot)
    report = validate(snapshot)
    assert report["sections"]["orthogroups"] == "absent"
    assert report["sections"]["kegg"] == "ready"
    assert state(snapshot) == before
    assert not (snapshot / "phenoradar_inputs").exists()


def test_multiple_runs_per_sample_key_are_rejected(snapshot):
    for relative in ["metadata/samples.tsv", "orthogroups/expression/tpm.tsv",
                     "orthogroups/expression/mapping_qc.tsv", "kegg/ko_tpm_sum.tsv",
                     "kegg/mapping_qc.tsv", "kegg/ko_support.tsv"]:
        path = snapshot / relative
        rows = read_tsv(path)
        rows.append({**rows[0], "run": "Replicate"})
        write_tsv(path, list(rows[0]), rows)
    before = state(snapshot)
    with pytest.raises(ValueError, match="multiple runs per sample key"):
        validate(snapshot)
    assert state(snapshot) == before


def test_completed_tree_and_gzipped_annotations_are_verified_in_place(snapshot):
    tree = completed_tree(snapshot, "((Plant_A,Plant_B-x),Plant_C);\n")
    with gzip.open(snapshot / "orthogroup_annotations.tsv.gz", "wt") as handle:
        handle.write("OG1\t3193\tA useful annotation\n")
    before = state(snapshot)
    report = validate(snapshot)
    assert str(tree) in report["files"]
    assert state(snapshot) == before


def test_filtered_outputs_are_selected_without_copying(snapshot):
    filtered = make_filtered(snapshot, ["Plant_C"])
    before = state(snapshot)
    report = validate(snapshot, exclusions=["Plant_C"])
    assert report["source"] == str(filtered)
    assert report["samples"] == 2
    assert validate(snapshot)["samples"] == 3
    assert state(snapshot) == before


def test_completion_requires_native_outputs_even_with_legacy_exports(snapshot):
    import analysis
    write_json(snapshot / 'run.json', {})
    write_json(snapshot / 'metadata/selection.json', {})
    pipeline = {'alignment': {'enabled': True}}
    analysis.validate_outputs(snapshot, {'pipeline': pipeline})
    legacy = snapshot / 'phenoradar_inputs/tpm.tsv'
    legacy.parent.mkdir()
    legacy.write_text('Previously collected TPM\n')
    (snapshot / 'orthogroups/expression/tpm.tsv').unlink()
    before = state(snapshot)
    with pytest.raises(ValueError, match='required analysis outputs.*orthogroups/expression/tpm.tsv'):
        analysis.validate_outputs(snapshot, {'pipeline': pipeline})
    assert state(snapshot) == before


def test_retired_target_does_not_touch_existing_outputs(snapshot, workflow_project, command_environment):
    snakemake = os.environ.get("SNAKEMAKE_BIN") or shutil.which("snakemake")
    if not snakemake:
        pytest.skip("Snakemake required")
    target = workflow_project / "results/test"
    target.parent.mkdir()
    shutil.copytree(snapshot, target)
    legacy = target / "phenoradar_inputs"
    legacy.mkdir()
    (legacy / "keep.txt").write_text("Keep existing export\n")
    before = state(target)
    cfg = workflow_project / "override.yaml"
    cfg.write_text(yaml.safe_dump({"run_name": "test"}))
    process = subprocess.run([snakemake, "--snakefile", str(ROOT / "workflow/Snakefile"),
                              "--configfile", str(cfg), "--cores", "1", "--", "phenoradar_inputs"],
                             cwd=workflow_project, env=command_environment({"python": sys.executable}),
                             text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    assert process.returncode != 0 and "MissingRuleException" in process.stdout
    assert state(target) == before
