import json
from pathlib import Path

import pytest

from aggregate_ko_tpm import GENE_FIELDS, HIT_FIELDS, aggregate
from common import file_record, read_tsv, write_json, write_tsv
from merge_kegg import merge


def annotate(root, species, assignments):
    """Small contract fixture: each gene gets (KO, hit status) pairs."""
    genes, hits = [], []
    for gene, pairs in assignments.items():
        accepted = [ko for ko, status in pairs if status == "accepted"]
        states = {status for _, status in pairs}
        status = ("unique" if len(accepted) == 1 else "ambiguous" if accepted else
                  "threshold_missing" if "threshold_missing" in states else
                  "below_threshold" if pairs else "unannotated")
        genes.append(dict(species=species, gene_id=gene, assignment_status=status,
                          accepted_ko_count=len(accepted), selected_ko=accepted[0] if len(accepted) == 1 else "",
                          terminal_stop_stripped=0))
        for ko, hit_status in pairs:
            hits.append(dict(species=species, gene_id=gene, ko=ko, score=10,
                             threshold="" if hit_status == "threshold_missing" else 5 if hit_status == "accepted" else 20,
                             evalue="1e-10", assignment_status=hit_status, accepted=int(hit_status == "accepted")))
    write_tsv(root / "genes.tsv", GENE_FIELDS, genes)
    write_tsv(root / "gene_kos.tsv", HIT_FIELDS, hits)
    (root / "detail.tsv").write_text("# fixture detail\n")
    provenance(root, species)


def provenance(root, species):
    write_json(root / "provenance.json", dict(schema_version=1, fingerprint="fixture", identity={"species": species},
                                             results=[file_record(root / name) for name in
                                                      ("genes.tsv", "gene_kos.tsv", "detail.tsv")]))


@pytest.fixture
def ko_inputs(tmp_path):
    annotation = tmp_path / "annotations"
    alpha = {f"Alpha_plant_g{i}": hits for i, hits in enumerate([
        [("K00001", "accepted"), ("K00009", "below_threshold")], [("K00001", "accepted")],
        [("K00002", "accepted"), ("K00003", "accepted")], [("K00004", "below_threshold")],
        [("K00005", "threshold_missing")], [], [("K00006", "accepted")], [("K00007", "accepted")],
    ], 1)}
    annotate(annotation / "Alpha_plant", "Alpha_plant", alpha)
    annotate(annotation / "Beta_sp_X", "Beta_sp-X", {"Beta_sp-X_g1": [("K00002", "accepted")], "Beta_sp-X_g2": []})
    samples = tmp_path / "samples.tsv"
    manifest = []
    for name, directory, run in [("Alpha_plant", "Alpha_plant", "A1"), ("Alpha_plant", "Alpha_plant", "A2"),
                                 ("Beta_sp-X", "Beta_sp_X", "B1")]:
        abundance = tmp_path / f"{run}.abundance.tsv"
        values = {1: 20, 2: 30, 3: 10, 4: 10, 5: 10, 6: 10, 8: 0, 9: 10}
        if run == "A2":
            values.update({1: 5, 2: 15})
        elif run == "B1":
            values = {1: 7, 2: 1}
        write_tsv(abundance, ["target_id", "tpm"],
                  [dict(target_id=f"{name}_g{i}", tpm=value) for i, value in values.items()])
        manifest.append(dict(species=name, odb_species=directory, run=run, abundance=str(abundance)))
    write_tsv(samples, list(manifest[0]), manifest)
    return dict(samples=samples, annotation_dir=annotation, run_dir=tmp_path / "runs", outdir=tmp_path / "merged")


def aggregate_runs(inputs, **options):
    for row in read_tsv(inputs["samples"]):
        run = row["run"]
        aggregate(inputs["samples"], run, inputs["annotation_dir"] / row["odb_species"],
                  inputs["run_dir"] / f"{run}.tsv", inputs["run_dir"] / f"{run}.qc.json", **options)


def test_explicit_drop_preserves_unique_assignment_and_support(ko_inputs):
    aggregate_runs(ko_inputs, ambiguity="drop")
    rows = {row["ko"]: row for row in read_tsv(ko_inputs["run_dir"] / "A1.tsv")}
    assert set(rows) == {"K00001", "K00006", "K00007"}
    assert float(rows["K00001"]["tpm_sum"]) == 50
    assert rows["K00001"]["annotated_genes"] == rows["K00001"]["quantified_genes"] == "2"
    assert rows["K00006"]["tpm_sum"] == ""
    assert rows["K00006"]["quantified_genes"] == "0"
    assert rows["K00007"]["tpm_sum"] == "0.0"
    assert rows["K00007"]["quantified_genes"] == "1"
    qc = json.loads((ko_inputs["run_dir"] / "A1.qc.json").read_text())
    assert qc["total_tpm"] == 100
    assert qc["retained_tpm"] == 50
    assert qc["retained_tpm_fraction"] == 0.5
    assert qc["unique_genes"] == 4
    assert qc["ambiguous_tpm"] == qc["below_threshold_tpm"] == qc["threshold_missing_tpm"] == 10
    assert qc["unannotated_tpm"] == qc["no_protein_tpm"] == 10
    assert qc["retained_targets"] == 3
    assert sum(float(row["tpm_sum"]) for row in rows.values() if row["tpm_sum"]) == 50
    assert qc["ko_tpm_sum"] == 50
    assert qc["quantified_assignments"] == 3
    merge(**ko_inputs)


def test_default_adds_full_tpm_to_each_accepted_ko_and_counts_coverage_once(ko_inputs):
    aggregate_runs(ko_inputs)
    rows = {row["ko"]: row for row in read_tsv(ko_inputs["run_dir"] / "A1.tsv")}
    assert set(rows) == {"K00001", "K00002", "K00003", "K00006", "K00007"}
    assert float(rows["K00001"]["tpm_sum"]) == 50
    for ko in ("K00002", "K00003"):
        assert float(rows[ko]["tpm_sum"]) == 10
        assert rows[ko]["annotated_genes"] == rows[ko]["quantified_genes"] == "1"
    assert rows["K00006"]["tpm_sum"] == ""  # no observation
    assert rows["K00007"]["tpm_sum"] == "0.0"  # measured zero
    qc = json.loads((ko_inputs["run_dir"] / "A1.qc.json").read_text())
    assert qc["ambiguity"] == "duplicate"
    assert qc["total_tpm"] == 100
    assert qc["retained_tpm"] == 60
    assert qc["retained_tpm_fraction"] == 0.6
    assert qc["retained_targets"] == 4
    assert qc["quantified_assignments"] == 5
    assert qc["ko_tpm_sum"] == 70
    assert qc["ambiguous_tpm"] == 10


def test_merge_runs_missing_evidence_and_selection(ko_inputs):
    aggregate_runs(ko_inputs)
    merge(**ko_inputs)
    wide = read_tsv(ko_inputs["outdir"] / "ko_tpm_sum_wide.tsv")
    assert [row["run"] for row in wide] == ["A1", "A2", "B1"]
    assert float(wide[0]["K00001"]) == 50
    assert float(wide[1]["K00001"]) == 20
    assert wide[0]["K00002"] == wide[0]["K00003"] == "10.0"
    assert wide[0]["K00006"] == wide[2]["K00001"] == ""
    assert wide[0]["K00007"] == "0.0"
    numeric = read_tsv(ko_inputs["outdir"] / "ko_tpm_sum.tsv")
    assert len(numeric) == 9
    assert all(row["tpm_sum"] != "" for row in numeric)
    assert len(read_tsv(ko_inputs["outdir"] / "ko_support.tsv")) == 11
    qc = read_tsv(ko_inputs["outdir"] / "mapping_qc.tsv")[0]
    assert float(qc["retained_tpm"]) == 60
    assert float(qc["ko_tpm_sum"]) == 70
    assert len(read_tsv(ko_inputs["outdir"] / "genes.tsv")) == 10  # species annotation not duplicated per run
    assert any(row["assignment_status"] == "accepted" and row["gene_id"] == "Alpha_plant_g3"
               for row in read_tsv(ko_inputs["outdir"] / "gene_kos.tsv"))
    manifest = read_tsv(ko_inputs["samples"])
    write_tsv(ko_inputs["samples"], list(manifest[0]), [manifest[-1]])
    # Stale, even corrupted unselected products are ignored.
    (ko_inputs["run_dir"] / "A1.tsv").write_text("broken")
    (ko_inputs["annotation_dir"] / "Alpha_plant/genes.tsv").write_text("broken")
    merge(**ko_inputs)
    assert [row["run"] for row in read_tsv(ko_inputs["outdir"] / "ko_tpm_sum_wide.tsv")] == ["B1"]
    assert {row["species"] for row in read_tsv(ko_inputs["outdir"] / "genes.tsv")} == {"Beta_sp-X"}


@pytest.mark.parametrize("ambiguity", ["duplicate", "drop", "error"])
@pytest.mark.parametrize("kind", ["ambiguous", "unannotated", "zero"])
def test_multi_ko_only_no_retained_or_all_zero_runs(ko_inputs, kind, ambiguity):
    rows = read_tsv(ko_inputs["samples"])
    sample = rows[-1]
    write_tsv(ko_inputs["samples"], list(sample), [sample])
    hits = [] if kind == "unannotated" else [("K00001", "accepted"), ("K00002", "accepted")]
    annotate(ko_inputs["annotation_dir"] / "Beta_sp_X", "Beta_sp-X",
             {"Beta_sp-X_g1": hits, "Beta_sp-X_g2": hits})
    if kind == "zero":
        values = read_tsv(sample["abundance"])
        for row in values:
            row["tpm"] = 0
        write_tsv(sample["abundance"], list(values[0]), values)
    if kind != "unannotated" and ambiguity == "error":
        with pytest.raises(ValueError, match="multiple accepted KOs"):
            aggregate_runs(ko_inputs, ambiguity=ambiguity)
        return
    aggregate_runs(ko_inputs, ambiguity=ambiguity)
    merge(**ko_inputs)
    report = json.loads((ko_inputs["run_dir"] / "B1.qc.json").read_text())
    if kind != "unannotated" and ambiguity == "duplicate":
        total = 0 if kind == "zero" else 8
        assert report["total_tpm"] == report["retained_tpm"] == total
        assert report["ko_tpm_sum"] == 2 * total
        assert report["retained_tpm_fraction"] == (None if kind == "zero" else 1)
        assert report["retained_targets"] == 2
        assert report["quantified_assignments"] == 4
        assert not report["no_retained_kos"]  # observed zero remains evidence
        numeric = read_tsv(ko_inputs["outdir"] / "ko_tpm_sum.tsv")
        assert {r["ko"] for r in numeric} == {"K00001", "K00002"}
        assert all(float(r["tpm_sum"]) == total for r in numeric)
    else:
        assert report["retained_tpm"] == report["ko_tpm_sum"] == 0
        assert report["no_retained_kos"]
        assert read_tsv(ko_inputs["outdir"] / "ko_tpm_sum.tsv") == []
        assert read_tsv(ko_inputs["outdir"] / "ko_tpm_sum_wide.tsv") == [{"species": "Beta_sp-X", "run": "B1"}]


def test_unobserved_multi_ko_has_missing_support_under_default_policy(ko_inputs):
    sample = read_tsv(ko_inputs["samples"])[0]
    values = read_tsv(sample["abundance"])
    write_tsv(sample["abundance"], list(values[0]), [r for r in values if r["target_id"] != "Alpha_plant_g3"])
    aggregate_runs(ko_inputs)
    merge(**ko_inputs)
    support = {r["ko"]: r for r in read_tsv(ko_inputs["run_dir"] / "A1.tsv")}
    for ko in ("K00002", "K00003"):
        assert support[ko]["annotated_genes"] == "1"
        assert support[ko]["quantified_genes"] == "0"
        assert support[ko]["tpm_sum"] == ""


def test_ambiguity_error_does_not_split(ko_inputs):
    with pytest.raises(ValueError, match="multiple accepted KOs"):
        aggregate(ko_inputs["samples"], "A1", ko_inputs["annotation_dir"] / "Alpha_plant",
                  ko_inputs["run_dir"] / "A1.tsv", ko_inputs["run_dir"] / "A1.qc.json", ambiguity="error")
    with pytest.raises(ValueError, match="splitting is unsupported"):
        aggregate(ko_inputs["samples"], "A1", "unused", "unused", "unused", ambiguity="split")


@pytest.mark.parametrize("value", ["nan", "inf", "-1", "abc", ""])
def test_invalid_abundance_numbers_rejected(ko_inputs, value):
    sample = read_tsv(ko_inputs["samples"])[0]
    values = read_tsv(sample["abundance"])
    values[0]["tpm"] = value
    write_tsv(sample["abundance"], list(values[0]), values)
    with pytest.raises(ValueError, match="invalid TPM"):
        aggregate_runs(ko_inputs)


@pytest.mark.parametrize("value", ["", "bad id", "../bad", "Alpha_plant_g2"])
def test_invalid_and_duplicate_target_ids_rejected(ko_inputs, value):
    sample = read_tsv(ko_inputs["samples"])[0]
    values = read_tsv(sample["abundance"])
    values[0]["target_id"] = value
    write_tsv(sample["abundance"], list(values[0]), values)
    with pytest.raises(ValueError, match="target_id"):
        aggregate_runs(ko_inputs)


def test_wrong_annotation_species_rejected(ko_inputs):
    with pytest.raises(ValueError, match="species does not match"):
        aggregate(ko_inputs["samples"], "A1", ko_inputs["annotation_dir"] / "Beta_sp_X", "unused", "unused")


def test_quantified_gene_from_other_species_rejected_at_merge(ko_inputs):
    sample = read_tsv(ko_inputs["samples"])[0]
    values = read_tsv(sample["abundance"])
    values[-1]["target_id"] = "Beta_sp-X_g1"
    write_tsv(sample["abundance"], list(values[0]), values)
    aggregate_runs(ko_inputs)
    with pytest.raises(ValueError, match="different species"):
        merge(**ko_inputs)


@pytest.mark.parametrize("change", ["abundance", "annotation", "run_table", "qc_species"])
def test_provenance_and_manifest_mismatches_rejected(ko_inputs, change):
    aggregate_runs(ko_inputs)
    if change == "abundance":
        path = Path(read_tsv(ko_inputs["samples"])[0]["abundance"])
        path.write_text(path.read_text().replace("20", "21"))
    elif change == "annotation":
        path = ko_inputs["annotation_dir"] / "Alpha_plant/genes.tsv"
        path.write_text(path.read_text().replace("K00001", "K00008"))
    elif change == "run_table":
        path = ko_inputs["run_dir"] / "A1.tsv"
        path.write_text(path.read_text().replace("50.0", "51.0"))
    else:
        path = ko_inputs["run_dir"] / "A1.qc.json"
        report = json.loads(path.read_text())
        report["species"] = "Wrong"
        write_json(path, report)
    with pytest.raises(ValueError, match="changed|does not match"):
        merge(**ko_inputs)


@pytest.mark.parametrize("field,value", [
    ("retained_targets", 5), ("retained_tpm", 70), ("retained_tpm_fraction", 0.7),
    ("quantified_assignments", 4), ("ko_tpm_sum", 60), ("ambiguity", "drop"),
])
def test_merge_rejects_confusing_gene_coverage_with_ko_totals(ko_inputs, field, value):
    aggregate_runs(ko_inputs)
    path = ko_inputs["run_dir"] / "A1.qc.json"
    report = json.loads(path.read_text())
    report[field] = value
    write_json(path, report)
    with pytest.raises(ValueError, match="disagrees"):
        merge(**ko_inputs)


def test_merge_rejects_dividing_multi_ko_tpm_even_with_matching_checksum(ko_inputs):
    aggregate_runs(ko_inputs)
    path = ko_inputs["run_dir"] / "A1.tsv"
    rows = read_tsv(path)
    for row in rows:
        if row["ko"] in {"K00002", "K00003"}:
            row["tpm_sum"] = 5
    write_tsv(path, list(rows[0]), rows)
    qc = path.with_suffix(".qc.json")
    report = json.loads(qc.read_text())
    report.update(result=file_record(path), ko_tpm_sum=60)
    write_json(qc, report)
    with pytest.raises(ValueError, match="TPM disagrees with input abundance"):
        merge(**ko_inputs)


def test_global_gene_ids_must_be_unique(ko_inputs):
    annotate(ko_inputs["annotation_dir"] / "Beta_sp_X", "Beta_sp-X", {"Alpha_plant_g1": [("K00001", "accepted")]})
    aggregate_runs(ko_inputs)
    with pytest.raises(ValueError, match="not globally unique"):
        merge(**ko_inputs)


@pytest.mark.parametrize("change", ["bad_ko", "duplicate_hit", "summary_count", "nan_score"])
def test_invalid_annotation_contract_rejected(ko_inputs, change):
    root = ko_inputs["annotation_dir"] / "Alpha_plant"
    if change == "summary_count":
        rows = read_tsv(root / "genes.tsv")
        rows[0]["accepted_ko_count"] = 2
        write_tsv(root / "genes.tsv", GENE_FIELDS, rows)
    else:
        rows = read_tsv(root / "gene_kos.tsv")
        if change == "bad_ko":
            rows[0]["ko"] = "K01"
        elif change == "duplicate_hit":
            rows.append(rows[0])
        else:
            rows[0]["score"] = "nan"
        write_tsv(root / "gene_kos.tsv", HIT_FIELDS, rows)
    provenance(root, "Alpha_plant")
    with pytest.raises(ValueError):
        aggregate_runs(ko_inputs)


def test_rounded_kofam_accepted_marker_is_authoritative(ko_inputs):
    root = ko_inputs["annotation_dir"] / "Alpha_plant"
    rows = read_tsv(root / "gene_kos.tsv")
    rows[0].update(score="10.0", threshold="10.04")
    write_tsv(root / "gene_kos.tsv", HIT_FIELDS, rows)
    provenance(root, "Alpha_plant")
    aggregate_runs(ko_inputs)
    assert float(read_tsv(ko_inputs["run_dir"] / "A1.tsv")[0]["tpm_sum"]) == 50
