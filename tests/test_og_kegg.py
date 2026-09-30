"""Median representatives, OG-level labels, and original-TPM aggregation."""
import json
from pathlib import Path

import pytest

from common import file_record, read_tsv, write_json, write_tsv
from mapping_fixtures import make_mapping
from og_kegg import ANNOTATION_LABEL, aggregate, assign
from run_kofam import run
from select_kegg_representatives import select
from test_kegg_workflow import fake_kofam_command, small_reference


@pytest.fixture
def og_inputs(tmp_path):
    root = tmp_path / "completed"
    proteins = root / "proteins"
    proteins.mkdir(parents=True)
    sequences = {
        "A": {"g1": "AA", "g2": "AXAA", "g3": "MK*K", "g4": "MK*", "g5": "AAAA", "g9": "AAA"},
        "B": {"g1": "AAAA", "g2": "AAAA", "g3": "*"},
        "C": {"g1": "AAAAAA", "g2": "AAAA"},
        "D": {"g1": "AAAAAAAAAA"},
    }
    genes, pairs, samples = [], [], []
    for species, entries in sequences.items():
        path = proteins / f"{species}_protein.fa"
        path.write_text("".join(f">{species}_{gene}\n{sequence}\n" for gene, sequence in entries.items()))
        for gene in entries:
            name = f"{species}_{gene}"
            genes.append((name, species))
            if gene != "g9":
                pairs.append((name, "OG" + gene[1:]))
        samples.append(dict(species=species, run=f"R{species}", odb_species=species,
                            abundance="unused.tsv", scientific_name=species, taxid="1"))
    path = root / "metadata/samples.tsv"
    write_tsv(path, list(samples[0]), samples)
    write_tsv(root / "metadata/metadata_high_busco.tsv", ["species", "run", "family"],
              [dict(species=r["species"], run=r["run"], family="Plants") for r in samples])
    write_tsv(root / "metadata/species_metadata.tsv", ["species", "C4", "contrast_pair_id", "family"],
              [dict(species=s, C4="", contrast_pair_id="", family="Plants") for s in sequences])
    mapping = make_mapping(root / "orthogroups/mapping/snapshot.json", genes, pairs)
    data = json.loads(mapping.read_text())
    for entry in data["proteins"]:
        species = entry["species"]
        digest = file_record(proteins / f"{species}_protein.fa")["sha256"]
        entry["sha256"] = data["tables"][species]["protein_sha256"] = digest
    write_json(mapping, data)
    return dict(samples=path, mapping=mapping, protein_dir=proteins,
                outdir=root / "kegg/representatives", batch_size=2)


@pytest.fixture
def og_annotations(og_inputs, tmp_path, monkeypatch):
    select(**og_inputs)
    reference = small_reference(tmp_path)
    command = fake_kofam_command(tmp_path)
    monkeypatch.setenv("FAKE_KOFAM_LOG", str(tmp_path / "events.txt"))
    monkeypatch.delenv("SLURM_CPUS_PER_TASK", raising=False)
    reps = og_inputs["outdir"]
    out = reps.parent
    plan = json.loads((reps / "provenance.json").read_text())
    for batch in plan["batches"]:
        name = batch["batch"]
        run(reps / "batches" / f"{name}.faa", ANNOTATION_LABEL,
            reference / "reference.json", out / "annotation" / name,
            tmp_path / "work" / name, command=str(command))
    assign(reps, out / "annotation", out)
    return out


def test_selects_median_length_with_deterministic_ties_and_no_extra_representatives(og_inputs):
    select(**og_inputs)
    rows = {r["orthogroup"]: r for r in read_tsv(og_inputs["outdir"] / "representatives.tsv")}
    assert rows["OG1"]["median_length"] == "5.0"
    assert rows["OG1"]["length"] == "6"
    assert rows["OG1"]["representative_gene_id"] == "C_g1"
    assert rows["OG2"]["representative_gene_id"] == "B_g2"  # fewer unknowns, then ID
    assert rows["OG3"]["selection_status"] == "no_valid_representative"
    assert rows["OG3"]["invalid_member_count"] == "2"
    assert rows["OG4"]["length"] == "2"  # one terminal stop removed
    plan = json.loads((og_inputs["outdir"] / "provenance.json").read_text())
    assert plan["orthogroups"] == 5 and plan["representatives"] == 4
    assert plan["unmapped_proteins"] == 1
    assert [b["orthogroups"] for b in plan["batches"]] == [["OG1", "OG2"], ["OG4", "OG5"]]
    before = (og_inputs["outdir"] / "representatives.tsv").read_bytes()
    select(**og_inputs)
    assert (og_inputs["outdir"] / "representatives.tsv").read_bytes() == before


def test_changed_protein_is_rejected_before_selection(og_inputs):
    (og_inputs["protein_dir"] / "A_protein.fa").write_text(">A_g1\nAAAA\n")
    with pytest.raises(ValueError, match="protein differs from ODB"):
        select(**og_inputs)


def test_empty_mapping_produces_empty_annotation_without_search(og_inputs):
    data = json.loads(og_inputs["mapping"].read_text())
    for name, entry in data["tables"].items():
        import gzip
        path = og_inputs["mapping"].parent / entry["table"]["path"]
        rows = []
        with gzip.open(path, "rt") as handle:
            for line in list(handle)[1:]:
                rows.append(line.split("\t")[0])
        with gzip.open(path, "wt") as handle:
            handle.write("gene_id\torthogroup\n")
            handle.writelines(gene + "\t\n" for gene in rows)
        entry["table"] = {**file_record(path), "path": entry["table"]["path"]}
    write_json(og_inputs["mapping"], data)
    select(**og_inputs)
    reps = og_inputs["outdir"]
    assert json.loads((reps / "provenance.json").read_text())["batches"] == []
    assign(reps, reps.parent / "no_annotation_directory", reps.parent)
    assert read_tsv(reps.parent / "orthogroups.tsv") == []
    assert read_tsv(reps.parent / "og_kos.tsv") == []


def write_expression(root, species, values, total=None):
    run = "R" + species
    path = root / f"{run}.tsv"
    fields = ["species", "run", "orthogroup", "tpm_sum", "tpm"]
    write_tsv(path, fields, [dict(species=species, run=run, orthogroup=og, tpm_sum=value, tpm=999999)
                             for og, value in values.items()])
    write_json(root / f"{run}.qc.json", dict(species=species, run=run, expression=file_record(path),
               retained_tpm=sum(values.values()), total_tpm=sum(values.values()) if total is None else total))


@pytest.fixture
def og_expression(tmp_path):
    root = tmp_path / "og_expression"
    write_expression(root, "A", {"OG1": 5, "OG2": 0, "OG3": 4, "OG4": 2, "OG5": 3}, total=24)
    write_expression(root, "B", {"OG1": 7, "OG2": 11, "OG3": 0})
    write_expression(root, "C", {"OG1": 0, "OG2": 20})
    write_expression(root, "D", {"OG1": 15})
    return root


def test_og_labels_keep_real_donors_and_all_statuses(og_annotations):
    groups = {r["orthogroup"]: r for r in read_tsv(og_annotations / "orthogroups.tsv")}
    assert groups["OG1"]["representative_gene_id"] == "C_g1"
    assert groups["OG1"]["selected_ko"] == "K00001"
    assert groups["OG2"]["assignment_status"] == "ambiguous"
    assert groups["OG3"]["assignment_status"] == "no_valid_representative"
    assert groups["OG4"]["assignment_status"] == "below_threshold"
    assert not (og_annotations / "genes.tsv").exists()
    assert not (og_annotations / "gene_kos.tsv").exists()


@pytest.mark.parametrize("policy", ["duplicate", "drop"])
def test_aggregates_original_og_tpm_preserving_zero_and_missing(og_inputs, og_annotations, og_expression, policy):
    aggregate(og_inputs["samples"], og_annotations, og_expression, og_annotations, ambiguity=policy)
    values = {(r["run"], r["ko"]): float(r["tpm_sum"]) for r in read_tsv(og_annotations / "ko_tpm_sum.tsv")}
    assert values["RA", "K00001"] == 5  # never the normalized value 999999
    assert values["RC", "K00001"] == 0
    if policy == "duplicate":
        assert values["RA", "K00002"] == 0
        assert values["RB", "K00002"] == values["RB", "K00003"] == 11
        assert ("RD", "K00002") not in values
        support = {(r["run"], r["ko"]): r for r in read_tsv(og_annotations / "ko_support.tsv")}
        assert support["RD", "K00002"]["tpm_sum"] == ""
        assert support["RA", "K00002"]["quantified_orthogroups"] == "1"
    else:
        assert {ko for run, ko in values} == {"K00001"}
    qc = {r["run"]: r for r in read_tsv(og_annotations / "mapping_qc.tsv")}
    assert float(qc["RA"]["odb_mapped_tpm"]) == 14
    assert float(qc["RA"]["retained_tpm_fraction"]) == pytest.approx(5 / 24)
    assert qc["RA"]["no_valid_representative_orthogroups"] == "1"
    assert json.loads((og_annotations / "provenance.json").read_text())["additional_representatives"] is False


def test_error_policy_rejects_multi_ko_og_even_with_zero_tpm(og_inputs, og_annotations, og_expression):
    with pytest.raises(ValueError, match="multiple accepted KOs"):
        aggregate(og_inputs["samples"], og_annotations, og_expression, og_annotations, ambiguity="error")


def test_all_zero_run_is_valid(og_inputs, og_annotations, og_expression):
    write_expression(og_expression, "A", {"OG1": 0, "OG2": 0})
    aggregate(og_inputs["samples"], og_annotations, og_expression, og_annotations)
    qc = read_tsv(og_annotations / "mapping_qc.tsv")[0]
    assert qc["retained_tpm_fraction"] == ""


def test_expression_corruption_is_rejected(og_inputs, og_annotations, og_expression):
    path = og_expression / "RA.tsv"
    path.write_text(path.read_text().replace("\t5\t", "\t6\t"))
    with pytest.raises(ValueError, match="file changed"):
        aggregate(og_inputs["samples"], og_annotations, og_expression, og_annotations)


def test_species_filter_preserves_shared_og_annotation_when_donor_is_removed(
        og_inputs, og_annotations, og_expression, tmp_path):
    from filter_species import export
    from phenoradar_inputs import export as collect
    aggregate(og_inputs["samples"], og_annotations, og_expression, og_annotations)
    source = og_annotations.parent
    out = source / "filtered"
    export(source, ["C"], out)
    assert (out / "kegg/orthogroups.tsv").read_bytes() == (source / "kegg/orthogroups.tsv").read_bytes()
    assert "C" not in {r["species"] for r in read_tsv(out / "kegg/ko_tpm_sum.tsv")}
    # The removed sample remains the recorded donor; no annotation is recomputed.
    assert read_tsv(out / "kegg/orthogroups.tsv")[0]["representative_species"] == "C"
    assert json.loads((out / "kegg/filter_qc.json").read_text())["representatives_reselected"] is False
    # A collector includes OG evidence and ignores stale direct-gene output.
    (source / "kegg/genes.tsv").write_text("stale\n")
    (source / "kegg/gene_kos.tsv").write_text("stale\n")
    destination = tmp_path / "inputs"
    collect(source, destination)
    assert (destination / "kegg/og_kos.tsv").is_file()
    assert (destination / "kegg/annotation_provenance.json").is_file()
    assert (destination / "kegg/representatives/representatives.tsv").is_file()
    assert not (destination / "kegg/genes.tsv").exists()
    assert not (destination / "kegg/gene_kos.tsv").exists()

    filtered_inputs = tmp_path / "filtered_inputs"
    collect(source, filtered_inputs, exclusions=["C"])
    assert {r["species"] for r in read_tsv(filtered_inputs / "kegg/ko_tpm_sum.tsv")} == {"A", "B", "D"}
    assert read_tsv(filtered_inputs / "kegg/orthogroups.tsv")[0]["representative_species"] == "C"
    assert (filtered_inputs / "kegg/source_provenance.json").is_file()
