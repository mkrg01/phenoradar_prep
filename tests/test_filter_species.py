"""Manual species exports preserve numeric values, sites, sources and branches."""
import json
import os
from pathlib import Path
import shutil
import sqlite3
import subprocess
import sys

import pytest
import yaml

from common import file_record, read_tsv, sha256, write_json, write_tsv
from filter_species import discover, export, fasta_records, validate_exclusions
from filter_species_phylogeny import read_tree

ROOT = Path(__file__).resolve().parents[1]
SPECIES = ["Plant_A", "Plant_B-x", "Plant_C", "Plant_D", "Plant_E"]
TREE = "((Plant_A:1,Plant_B-x:2)0.9:3,(Plant_C:4,(Plant_D:5,Plant_E:6)0.8:7)0.7:8);"


@pytest.fixture
def snapshot(tmp_path):
    source = tmp_path / "completed"
    rows = [dict(species=s, scientific_name=s.replace("_", " "), taxid=10+i,
                 odb_species=s.replace("-", "_"), run=f"R{i}") for i, s in enumerate(SPECIES)]
    rows.append({**rows[0], "run": "R5"})
    write_tsv(source / "metadata/samples.tsv", list(rows[0]), rows)
    write_tsv(source / "metadata/metadata_all.tsv", list(rows[0]), rows + [{**rows[0], "species": "Not_selected", "run": "Old"}])
    write_tsv(source / "metadata/metadata_high_busco.tsv", list(rows[0]), rows)
    write_json(source / "metadata/selection.json", {"selected_species": 5})
    write_json(source / "run.json", {"original_configuration": True})
    traits = tmp_path / "traits.tsv"
    write_tsv(traits, ["species", "C4", "other"], [{"species": s.replace("_", " "), "C4": str(i % 2), "other": ""}
                                                  for i, s in enumerate(SPECIES)])
    genes = [("Plant_A_g99", SPECIES[0]), ("Plant_B-x_g2", SPECIES[1]), ("Plant_B-x_g3", SPECIES[1]),
             ("Plant_C_g4", SPECIES[2]), ("Plant_D_g5", SPECIES[3]), ("Plant_E_g6", SPECIES[4])]
    pairs = [(g, "OGshared") for g, s in genes] + [("Plant_A_g99", "OGempty"), ("Plant_B-x_g2", "OGamb")]
    folder = source / "orthogroups/mapping"
    folder.mkdir(parents=True)
    with sqlite3.connect(folder / "mappings.sqlite") as db:
        db.executescript("CREATE TABLE genes(query TEXT PRIMARY KEY,species TEXT); CREATE TABLE mappings(query TEXT,og TEXT,PRIMARY KEY(query,og));")
        db.executemany("INSERT INTO genes VALUES (?,?)", genes)
        db.executemany("INSERT INTO mappings VALUES (?,?)", pairs)
    write_tsv(folder / "gene_orthogroups.tsv", ["#query", "ODB_OG"], [{"#query": g, "ODB_OG": og} for g, og in pairs])
    write_json(folder / "merge_qc.json", {"unique_gene_og_pairs": len(pairs)})
    (source / "proteins").mkdir()
    for row in rows[:5]:
        (source / "proteins" / f"{row['odb_species']}_protein.fa").write_text(
            "".join(f">{gene} original description\nAACD\n" for gene, s in genes if s == row["species"]))
    for name, value in [("tpm", "1000000.000000"), ("tpm_sum", "0.0000001")]:
        write_tsv(source / f"orthogroups/expression/{name}.tsv", ["species", "run", "orthogroup", name],
                  [{"species": r["species"], "run": r["run"], "orthogroup": "OGshared", name: value} for r in rows])
        write_tsv(source / f"orthogroups/expression/{name}_wide.tsv", ["species", "run", "OGshared", "OGzero"],
                  [{"species": r["species"], "run": r["run"], "OGshared": value, "OGzero": "0"} for r in rows])
    write_tsv(source / "orthogroups/expression/mapping_qc.tsv", ["species", "run", "fraction"],
              [{"species": r["species"], "run": r["run"], "fraction": ".000001"} for r in rows])
    alignments = []
    for og in sorted({og for g, og in pairs}):
        path = source / "orthogroups/alignments" / f"{og}.faa"
        path.parent.mkdir(exist_ok=True)
        records = [(g, "aA----" if s == SPECIES[0] else "--C--D") for g, s in genes if (g, og) in pairs]
        path.write_text("".join(f">{g} full header\n{seq}\n" for g, seq in records))
        alignments.append(dict(orthogroup=og, alignment=file_record(path)))
    write_json(source / "orthogroups/alignments/provenance.json", {"alignments": alignments})
    write_tsv(source / "kegg/genes.tsv", ["species", "gene_id", "assignment_status"],
              [dict(species=s, gene_id=g, assignment_status="unique") for g, s in genes])
    write_tsv(source / "kegg/gene_kos.tsv", ["species", "gene_id", "ko", "accepted"],
              [dict(species=s, gene_id=g, ko="K00001", accepted="1") for g, s in genes])
    for name in ["ko_tpm_sum", "ko_support", "mapping_qc"]:
        write_tsv(source / f"kegg/{name}.tsv", ["species", "run", "ko", "tpm_sum"],
                  [dict(species=r["species"], run=r["run"], ko="K00001", tpm_sum="0.00000") for r in rows])
    write_tsv(source / "kegg/ko_tpm_sum_wide.tsv", ["species", "run", "K00001", "K00002"],
              [dict(species=r["species"], run=r["run"], K00001="0", K00002="") for r in rows])
    write_tsv(source / "kegg/ko_modules.tsv", ["ko", "module"], [dict(ko="K00001", module="M00001")])
    folder = source / "phylogeny/all"
    folder.mkdir(parents=True)
    (folder / "species_tree.nwk").write_text(TREE + "\n")
    (folder / "gene_trees.nwk").write_text(TREE + "\n((Plant_A:1,Plant_B-x:2):3,(Plant_C:4,Plant_D:5):6);\n")
    write_json(folder / "species_tree.json", {"species": 5, "outgroup": "Plant_A", "input": file_record(folder / "gene_trees.nwk"),
               "branch_length_unit": "substitutions_per_site"})
    write_json(folder / "gene_trees.json", {"status": "retained", "retained": [{"marker": "marker1"}, {"marker": "marker2"}]})
    write_tsv(folder / "species_coverage.tsv", ["species", "gene_trees"], [dict(species=s, gene_trees=1) for s in SPECIES])
    (folder / "alignments").mkdir()
    (folder / "alignments/marker1.faa").write_text("".join(f">{s}\n{'AA----' if s == SPECIES[0] else '--CC--'}\n" for s in SPECIES))
    (folder / "alignments/marker1.columns.tsv").write_text("trimmed_column_1based\tfamsa_column_1based\n1\t7\n")
    (folder / "dating").mkdir()
    (folder / "dating/species_tree.dated.nwk").write_text("(Plant_A:4,(Plant_B-x:3,(Plant_C:2,(Plant_D:1,Plant_E:1):1):1):1);\n")
    for relative in ["phylogeny/phenotyped", "phylogeny/representatives", "phylogeny/all/taxonomy_audit"]:
        path = source / relative
        path.mkdir(parents=True)
        (path / "original.txt").write_text("Plant_A must stay in this historical result\n")
    return source, traits


def state(folder):
    return {str(p.relative_to(folder)): (sha256(p), p.stat().st_mtime_ns)
            for p in folder.rglob("*") if p.is_file()}


def test_filtered_metadata_and_alignments_feed_phenoradar(snapshot, tmp_path):
    from phenoradar_inputs import export as collect_inputs
    source, traits = snapshot
    path = source / "metadata/metadata_high_busco.tsv"
    rows = [{**row, "family": "Plantaceae"} for row in read_tsv(path)]
    write_tsv(path, list(rows[0]), rows)
    filtered = source / "filtered"
    export(source, ["Plant_A"], filtered, traits)
    metadata = read_tsv(filtered / "metadata/species_metadata.tsv")
    assert {r["species"] for r in metadata} == set(SPECIES[1:])
    assert all(r["family"] == "Plantaceae" and r["contrast_pair_id"] == "" for r in metadata)
    out = tmp_path / "phenoradar_inputs"
    collect_inputs(source, out, exclusions=["Plant_A"],
                   settings=dict(orthogroups=True, kegg=False, alignments=True,
                                 tree=str(filtered / "phylogeny/all/species_tree.pruned.nwk")))
    assert (out / "tpm.tsv").resolve() == filtered / "orthogroups/expression/tpm.tsv"
    assert (out / "alignments/OGshared.faa").is_symlink()
    assert not (out / "alignments/OGempty.faa").exists()
    assert (out / "species_tree.nwk").is_symlink()


@pytest.mark.parametrize("value", [None, "Plant_A", {"species": "Plant_A"}, [1], ["Plant A"], ["../Plant_A"], ["Plant_A", "Plant_A"]])
def test_exclusions_require_exact_list(value):
    with pytest.raises(ValueError):
        validate_exclusions(value)


def test_export_all_outputs_preserves_values_and_sources(snapshot, tmp_path):
    source, traits = snapshot
    before = state(source)
    out = tmp_path / "filtered"
    summary = export(source, ["Plant_A"], out, traits)
    assert state(source) == before
    assert summary["counts"]["species"] == dict(before=5, after=4)
    assert summary["counts"]["runs"] == dict(before=6, after=4)
    assert {r["run"] for r in read_tsv(out / "excluded_samples.tsv")} == {"R0", "R5"}
    assert {r["species"] for r in read_tsv(out / "metadata/samples.tsv")} == set(SPECIES[1:])
    assert "Not_selected" not in (out / "metadata/metadata_all.tsv").read_text()
    for branch in ["orthogroups/expression", "kegg"]:
        for path in (out / branch).glob("*.tsv"):
            original = read_tsv(source / branch / path.name)
            if original and "species" in original[0]:
                assert read_tsv(path) == [r for r in original if r["species"] != "Plant_A"]
    wide = read_tsv(out / "kegg/ko_tpm_sum_wide.tsv")
    assert all(row["K00001"] == "0" and row["K00002"] == "" for row in wide)
    with sqlite3.connect(out / "orthogroups/mapping/mappings.sqlite") as db:
        assert db.execute("SELECT query FROM genes WHERE species='Plant_A'").fetchall() == []
        assert db.execute("SELECT og FROM mappings WHERE query='Plant_B-x_g2' ORDER BY og").fetchall() == [("OGamb",), ("OGshared",)]
        assert db.execute("PRAGMA foreign_key_check").fetchall() == []
    assert not (out / "orthogroups/alignments/OGempty.faa").exists()
    assert not (out / "orthogroups/alignments/members.tsv").exists()
    assert list(fasta_records(out / "orthogroups/alignments/OGshared.faa")) == [r for r in fasta_records(source / "orthogroups/alignments/OGshared.faa") if not r[0].startswith("Plant_A_g99 ")]
    assert all(len(seq) == 6 for header, seq in fasta_records(out / "orthogroups/alignments/OGshared.faa"))
    assert (out / "proteins/Plant_B_x_protein.fa").is_symlink()
    assert not (out / "proteins/Plant_A_protein.fa").exists()
    old = read_tree(TREE, set(SPECIES))
    new = read_tree((out / "phylogeny/all/species_tree.pruned.nwk").read_text(), set(SPECIES[1:]), exact=True)
    for a in SPECIES[1:]:
        for b in SPECIES[1:]:
            assert old.get_distance(a, b) == pytest.approx(new.get_distance(a, b))
    assert all(not n.name for n in new.traverse() if not n.is_leaf)
    report = json.loads((out / "phylogeny/all/pruning.json").read_text())
    assert report["reestimated"] is False and report["rooting"] == "original_outgroup_removed_requires_review"
    assert not (out / "phylogeny/all/species_tree.json").exists()
    assert not (out / "phylogeny/all/dating/node_ages.tsv").exists()
    assert sha256(out / "phylogeny/all/alignments/marker1.columns.tsv") == sha256(source / "phylogeny/all/alignments/marker1.columns.tsv")
    assert all(len(seq) == 6 for header, seq in fasta_records(out / "phylogeny/all/alignments/marker1.faa"))
    for branch in ["phylogeny/phenotyped", "phylogeny/representatives", "phylogeny/all/taxonomy_audit"]:
        assert not (out / branch).exists()
    assert all(sha256(out / record["path"]) == record["sha256"] for record in summary["outputs"])


def test_rerun_always_uses_originals_and_removes_stale_exports(snapshot, tmp_path):
    source, traits = snapshot
    out = tmp_path / "filtered"
    export(source, ["Plant_A"], out, traits)
    export(source, [], out, traits)
    assert (out / "orthogroups/alignments/OGempty.faa").is_file()
    assert (out / "proteins/Plant_A_protein.fa").is_symlink()
    assert len(read_tsv(out / "metadata/samples.tsv")) == 6
    export(source, SPECIES[:4], out, traits)
    assert not (out / "phylogeny/all/species_tree.pruned.nwk").exists()
    assert not (out / "phylogeny/all/dating").exists()
    assert (out / "phylogeny/all/gene_trees.pruned.nwk").read_text() == ""
    assert all(r["status"] == "fewer_than_two_tips" for r in read_tsv(out / "phylogeny/all/gene_trees.tsv"))


@pytest.mark.parametrize("excluded", [["Unknown"], SPECIES])
def test_unknown_species_and_empty_selection_preserve_prior_export(snapshot, tmp_path, excluded):
    source, traits = snapshot
    out = tmp_path / "filtered"
    export(source, [], out)
    before = state(out)
    with pytest.raises(ValueError):
        export(source, excluded, out)
    assert state(out) == before


def test_incomplete_branch_is_reported_and_never_built(snapshot, tmp_path):
    source, traits = snapshot
    (source / "phylogeny/all/species_tree.nwk").unlink()
    out = tmp_path / "filtered"
    result = export(source, ["Plant_A"], out)
    assert result["sections"]["phylogeny"]["status"] == "incomplete"
    assert not (out / "phylogeny/all").exists()
    assert (out / "orthogroups/expression/tpm.tsv").is_file()


def test_alignments_without_odb_or_membership_table_use_gene_ids(snapshot, tmp_path):
    source, traits = snapshot
    shutil.rmtree(source / "orthogroups/mapping")
    out = tmp_path / "filtered"
    result = export(source, ["Plant_B-x"], out)
    assert result["sections"]["alignments"]["status"] == "ready"
    assert not (out / "orthogroups/alignments/members.tsv").exists()
    assert not (out / "orthogroups/alignments/OGamb.faa").exists()
    assert list(fasta_records(out / "orthogroups/alignments/OGshared.faa")) == [
        r for r in fasta_records(source / "orthogroups/alignments/OGshared.faa") if not r[0].startswith("Plant_B-x_g")]


def test_alignment_without_completion_record_is_incomplete(snapshot, tmp_path):
    source, traits = snapshot
    (source / "orthogroups/alignments/provenance.json").unlink()
    result = export(source, [], tmp_path / "filtered")
    assert result["sections"]["alignments"]["status"] == "incomplete"
    assert not (tmp_path / "filtered/orthogroups/alignments").exists()


@pytest.mark.parametrize("damage", ["format", "unknown_species", "duplicate", "empty", "unequal_lengths", "missing_gene", "missing_og", "duplicate_og"])
def test_alignment_id_and_membership_checks_use_fasta(snapshot, tmp_path, damage):
    source, traits = snapshot
    path = source / "orthogroups/alignments/OGshared.faa"
    records = list(fasta_records(path))
    if damage == "format":
        records[0] = ("arbitrary_id", records[0][1])
    elif damage == "unknown_species":
        records[0] = ("Unknown_g1", records[0][1])
    elif damage == "duplicate":
        records.append(records[0])
    elif damage == "empty":
        records = []
    elif damage == "unequal_lengths":
        records[0] = (records[0][0], "AA")
    elif damage == "missing_gene":
        records.pop()
    path.write_text("".join(f">{header}\n{seq}\n" for header, seq in records))
    report_path = source / "orthogroups/alignments/provenance.json"
    report = json.loads(report_path.read_text())
    for record in report["alignments"]:
        if record["orthogroup"] == "OGshared":
            record["alignment"] = file_record(path)
    if damage == "missing_og":
        report["alignments"] = [r for r in report["alignments"] if r["orthogroup"] != "OGamb"]
    elif damage == "duplicate_og":
        report["alignments"].append(report["alignments"][0])
    write_json(report_path, report)
    # The recorded hash is consistent: validation must inspect IDs/content.
    out = tmp_path / "filtered"
    with pytest.raises(ValueError):
        export(source, ["Plant_A"], out)
    assert not out.exists()


@pytest.mark.parametrize("damage", ["run_identity", "missing_run", "alignment_hash", "ko_owner", "odb_owner", "odb_pair", "gene_tree_hash"])
def test_inconsistent_inputs_fail_without_changing_originals(snapshot, tmp_path, damage):
    source, traits = snapshot
    if damage in {"run_identity", "missing_run"}:
        path = source / "orthogroups/expression/tpm_wide.tsv"
        rows = read_tsv(path)
        if damage == "run_identity":
            rows[1]["species"] = "Plant_A"
        else:
            rows.pop()
        write_tsv(path, list(rows[0]), rows)
    elif damage == "alignment_hash":
        with (source / "orthogroups/alignments/OGshared.faa").open("a") as handle:
            handle.write(">new_gene\nAAAAAA\n")
    elif damage == "ko_owner":
        path = source / "kegg/gene_kos.tsv"
        rows = read_tsv(path)
        rows[1]["species"] = "Plant_A"
        write_tsv(path, list(rows[0]), rows)
    elif damage == "odb_owner":
        with sqlite3.connect(source / "orthogroups/mapping/mappings.sqlite") as db:
            db.execute("UPDATE genes SET species='Plant_A' WHERE query='Plant_B-x_g2'")
    elif damage == "odb_pair":
        with sqlite3.connect(source / "orthogroups/mapping/mappings.sqlite") as db:
            db.execute("UPDATE mappings SET og='changed_OG' WHERE query='Plant_B-x_g2' AND og='OGamb'")
    else:
        path = source / "phylogeny/all/gene_trees.nwk"
        path.write_text(path.read_text().replace("Plant_A:1", "Plant_A:2"))
    before = state(source)
    with pytest.raises(ValueError):
        export(source, ["Plant_A"], tmp_path / "filtered")
    assert state(source) == before
    assert not (tmp_path / "filtered").exists()


def test_output_overlap_and_unowned_directory_rejected(snapshot, tmp_path):
    source, traits = snapshot
    for out in [source, source / "orthogroups/expression", tmp_path]:
        with pytest.raises(ValueError, match="overlap"):
            export(source, [], out)
    out = tmp_path / "unrelated"
    out.mkdir()
    (out / "keep.txt").write_text("retain")
    with pytest.raises(ValueError, match="not owned"):
        export(source, [], out)
    linked = tmp_path / "linked"
    linked.symlink_to(out, target_is_directory=True)
    with pytest.raises(ValueError, match="symlink"):
        export(source, [], linked)
    output = tmp_path / "filtered"
    export(source, ["Plant_A"], output)
    with pytest.raises(ValueError, match="original analysis"):
        export(output, [], tmp_path / "second_filter")


def test_full_snakefile_exports_frozen_results_without_upstream_inputs(snapshot, workflow_project, command_environment):
    snakemake = os.environ.get("SNAKEMAKE_BIN") or shutil.which("snakemake")
    if not snakemake:
        pytest.skip("Snakemake required")
    source, traits = snapshot
    project = workflow_project
    target = project / "results/test"
    target.parent.mkdir()
    shutil.copytree(source, target)
    # Real producer rules are loaded; no raw metadata/CDS/BUSCO/taxonomy exists.
    # Mapping/phylogeny planning would fail if export requested any producer.
    original = state(target)
    cfg = project / "override.yaml"
    environment = command_environment({"python": sys.executable})
    argv = [snakemake, "--snakefile", str(ROOT / "workflow/Snakefile"), "--configfile", str(cfg), "--cores", "1", "--", "filter_species"]
    def run(excluded):
        cfg.write_text(yaml.safe_dump({"analysis": "test", "exclude_species": excluded,
                                      "inputs": {"species_trait": str(traits)}}))
        process = subprocess.run(argv, cwd=project, env=environment, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        assert process.returncode == 0, process.stdout + "\n" + "\n".join(p.read_text() for p in (project / "logs").rglob("*.log"))
        return process.stdout
    run(["Plant_A"])
    out = target / "filtered"
    assert {r["species"] for r in read_tsv(out / "metadata/samples.tsv")} == set(SPECIES[1:])
    assert "Nothing to be done" in run(["Plant_A"])
    run(["Plant_C"])
    assert "Plant_A" in {r["species"] for r in read_tsv(out / "metadata/samples.tsv")}
    # Updating completed source tables must invalidate the subset, independently
    # of the exclusion list. A new optional branch would change input inventory.
    old = state(out)
    changed = target / "orthogroups/expression/tpm_wide.tsv"
    changed.write_text(changed.read_text().replace("\t0\n", "\t0.00000\n"))
    original["orthogroups/expression/tpm_wide.tsv"] = (sha256(changed), changed.stat().st_mtime_ns)
    run(["Plant_C"])
    assert state(out) != old
    after = {p: v for p, v in state(target).items() if not p.startswith("filtered/")}
    assert after == original
    assert not (project / "work").exists()
