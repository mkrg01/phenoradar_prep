#!/usr/bin/env python3
"""Assign representative KOs to OGs and aggregate original OG TPM without expansion."""
import argparse
from collections import Counter, defaultdict
from contextlib import ExitStack
import csv
import json
import math
from pathlib import Path
import re

from aggregate_ko_tpm import (AMBIGUITY_POLICIES, HIT_FIELDS, load_annotation, number,
                              selected_samples, summed, table, verify_file)
from common import atomic_writer, file_record, now, write_json, write_tsv
from select_kegg_representatives import FIELDS as REPRESENTATIVE_FIELDS

ANNOTATION_LABEL = "orthogroup_representatives"
OG_FIELDS = REPRESENTATIVE_FIELDS + ["assignment_status", "accepted_ko_count", "selected_ko"]
OG_HIT_FIELDS = ["orthogroup", "representative_gene_id"] + HIT_FIELDS[2:]
KO_FIELDS = ["species", "run", "ko", "tpm_sum", "quantified_orthogroups"]
QC_FIELDS = ["species", "run", "annotation_scope", "ambiguity", "reference_orthogroups",
             "observed_orthogroups", "retained_orthogroups", "ambiguous_orthogroups",
             "unannotated_orthogroups", "no_valid_representative_orthogroups",
             "total_tpm", "odb_mapped_tpm", "retained_tpm", "retained_tpm_fraction",
             "quantified_assignments", "ko_tpm_sum"]
ANNOTATION_FILES = ["orthogroups.tsv", "og_kos.tsv"]
EXPRESSION_FILES = ["ko_tpm_sum.tsv", "ko_tpm_sum_wide.tsv", "ko_support.tsv", "mapping_qc.tsv"]


def selection(directory):
    root = Path(directory)
    report = json.loads((root / "provenance.json").read_text())
    if report.get("method") != "median_length_one_representative":
        raise ValueError("unsupported representative selection")
    for entry in report["results"]:
        relative = Path(entry["path"])
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError("invalid representative artifact path")
        verify_file(root / relative, entry)
    rows = table(root / "representatives.tsv", REPRESENTATIVE_FIELDS)
    if len({r["orthogroup"] for r in rows}) != len(rows):
        raise ValueError("duplicate representative orthogroup")
    return report, rows


def assign(representatives, annotation_dir, outdir):
    plan, rows = selection(representatives)
    summaries, hits, receipts = {}, [], []
    reference_ids = set()
    for batch in plan["batches"]:
        name = batch["batch"]
        genes, batch_hits, receipt = load_annotation(Path(annotation_dir) / name, ANNOTATION_LABEL)
        if set(genes) != set(batch["orthogroups"]) or summaries.keys() & genes.keys():
            raise ValueError("annotation batch differs from representative selection")
        report = json.loads(Path(receipt["path"]).read_text())
        expected = Path(representatives) / "batches" / f"{name}.faa"
        verify_file(expected, report["identity"]["protein"])
        reference_ids.add(report["identity"]["reference_id"])
        summaries.update(genes)
        hits.extend(batch_hits)
        receipts.append(receipt)
    if len(reference_ids) > 1:
        raise ValueError("representative annotations use different KEGG references")
    selected = {r["orthogroup"] for r in rows if r["selection_status"] == "selected"}
    if set(summaries) != selected:
        raise ValueError("missing or unexpected representative annotation")
    by_og = {r["orthogroup"]: r for r in rows}
    out = Path(outdir)
    write_tsv(out / "orthogroups.tsv", OG_FIELDS, (
        {**row, **{key: summaries[row["orthogroup"]][key] if row["orthogroup"] in summaries
                    else ("no_valid_representative" if key == "assignment_status"
                          else 0 if key == "accepted_ko_count" else "")
                   for key in ["assignment_status", "accepted_ko_count", "selected_ko"]}}
        for row in rows))
    write_tsv(out / "og_kos.tsv", OG_HIT_FIELDS, (
        {"orthogroup": hit["gene_id"],
         "representative_gene_id": by_og[hit["gene_id"]]["representative_gene_id"],
         **{key: hit[key] for key in HIT_FIELDS[2:]}}
        for hit in sorted(hits, key=lambda r: (r["gene_id"], r["ko"]))))
    write_json(out / "annotation_provenance.json", {
        "schema_version": 1, "annotation_scope": "orthogroup",
        "method": "median_length_one_representative", "created_at": now(),
        "representative_selection": file_record(Path(representatives) / "provenance.json"),
        "annotation_batches": receipts, "reference_ids": sorted(reference_ids),
        "results": [file_record(out / name) for name in ANNOTATION_FILES],
    })


def load_og_annotations(directory):
    root = Path(directory)
    provenance = root / "annotation_provenance.json"
    report = json.loads(provenance.read_text())
    if report.get("annotation_scope") != "orthogroup" or report.get("method") != "median_length_one_representative":
        raise ValueError("unsupported OG annotation scope or method")
    records = {Path(r["path"]).name: r for r in report["results"]}
    if len(records) != len(report["results"]) or set(records) != set(ANNOTATION_FILES):
        raise ValueError("incomplete OG annotation")
    for name, entry in records.items():
        verify_file(root / name, entry)
    rows = table(root / "orthogroups.tsv", OG_FIELDS)
    groups = {r["orthogroup"]: r for r in rows}
    if len(groups) != len(rows):
        raise ValueError("duplicate annotated orthogroup")
    kos, hit_statuses, pairs = defaultdict(set), defaultdict(set), set()
    for hit in table(root / "og_kos.tsv", OG_HIT_FIELDS):
        og, ko = hit["orthogroup"], hit["ko"]
        if og not in groups or (og, ko) in pairs:
            raise ValueError("unknown or duplicate OG/KO pair")
        if hit["representative_gene_id"] != groups[og]["representative_gene_id"]:
            raise ValueError("KO representative differs from selected representative")
        if not re.fullmatch(r"K\d{5}", ko):
            raise ValueError(f"invalid KO identifier: {ko}")
        pairs.add((og, ko))
        status = hit["assignment_status"]
        if status not in {"accepted", "below_threshold", "threshold_missing"} or hit["accepted"] != str(int(status == "accepted")):
            raise ValueError("invalid representative hit status")
        number(hit["score"], "score", nonnegative=False)
        number(hit["evalue"], "evalue")
        if status == "threshold_missing":
            if hit["threshold"]:
                raise ValueError("threshold_missing hit has threshold")
        else:
            number(hit["threshold"], "threshold", nonnegative=False)
        if status == "accepted":
            kos[og].add(ko)
        hit_statuses[og].add(status)
    for og, row in groups.items():
        if row["selection_status"] not in {"selected", "no_valid_representative"}:
            raise ValueError(f"invalid representative selection status: {og}")
        accepted = kos[og]
        expected = ("no_valid_representative" if row["selection_status"] == "no_valid_representative"
                    else "unique" if len(accepted) == 1 else "ambiguous" if accepted
                    else "threshold_missing" if "threshold_missing" in hit_statuses[og]
                    else "below_threshold" if hit_statuses[og] else "unannotated")
        if (row["assignment_status"] != expected or int(row["accepted_ko_count"]) != len(accepted)
                or row["selected_ko"] != (next(iter(accepted)) if len(accepted) == 1 else "")):
            raise ValueError(f"OG summary differs from KO assignments: {og}")
        if expected == "no_valid_representative" and hit_statuses[og]:
            raise ValueError("KO hit without a valid representative")
    return groups, dict(kos), file_record(provenance)


def aggregate(samples, annotation_dir, og_run_dir, outdir, ambiguity="duplicate"):
    if ambiguity not in AMBIGUITY_POLICIES:
        raise ValueError("unsupported KO ambiguity policy")
    manifest = selected_samples(samples)
    groups, assignments, annotation = load_og_annotations(annotation_dir)
    eligible = {og: kos for og, kos in assignments.items()
                if kos and (ambiguity == "duplicate" or len(kos) == 1)}
    all_kos = sorted({ko for kos in eligible.values() for ko in kos})
    out = Path(outdir)
    reports, sources = [], []
    with ExitStack() as stack:
        def writer(name, fields):
            handle = stack.enter_context(atomic_writer(out / name))
            result = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            result.writeheader()
            return result
        long = writer("ko_tpm_sum.tsv", ["species", "run", "ko", "tpm_sum"])
        support = writer("ko_support.tsv", KO_FIELDS)
        wide = writer("ko_tpm_sum_wide.tsv", ["species", "run", *all_kos])
        for sample in manifest:
            species, run = sample["species"], sample["run"]
            path = Path(og_run_dir) / f"{run}.tsv"
            qc_path = Path(og_run_dir) / f"{run}.qc.json"
            qc = json.loads(qc_path.read_text())
            if qc.get("species") != species or qc.get("run") != run:
                raise ValueError(f"OG expression QC differs from sample: {run}")
            verify_file(path, qc.get("expression", {}))
            values = {}
            for row in table(path, ["species", "run", "orthogroup", "tpm_sum"]):
                og = row["orthogroup"]
                if row["species"] != species or row["run"] != run or og not in groups or og in values:
                    raise ValueError(f"invalid or duplicate OG expression row: {run}/{og}")
                values[og] = number(row["tpm_sum"], f"OG TPM: {run}/{og}")
            mapped = summed(values.values())
            total = number(qc["total_tpm"], "total input TPM")
            if not math.isclose(mapped, number(qc["retained_tpm"], "ODB retained TPM"), abs_tol=1e-8):
                raise ValueError(f"OG TPM differs from original-TPM QC: {run}")
            if mapped > total and not math.isclose(mapped, total, abs_tol=1e-8):
                raise ValueError(f"ODB TPM exceeds total input TPM: {run}")
            ambiguous = {og for og in values if len(assignments.get(og, ())) > 1}
            if ambiguity == "error" and ambiguous:
                raise ValueError(f"{run}: quantified orthogroups have multiple accepted KOs")
            contributions = defaultdict(list)
            retained = set()
            for og, value in values.items():
                for ko in eligible.get(og, ()):
                    contributions[ko].append(value)
                    retained.add(og)
            summed_kos = {ko: summed(vals) for ko, vals in contributions.items()}
            for ko in all_kos:
                value = summed_kos.get(ko, "")
                support.writerow(dict(species=species, run=run, ko=ko, tpm_sum=value,
                                      quantified_orthogroups=len(contributions[ko])))
                if ko in summed_kos:
                    long.writerow(dict(species=species, run=run, ko=ko, tpm_sum=value))
            wide.writerow({"species": species, "run": run, **{ko: summed_kos.get(ko, "") for ko in all_kos}})
            retained_tpm = summed(values[og] for og in retained)
            statuses = Counter(groups[og]["assignment_status"] for og in values)
            reports.append(dict(zip(QC_FIELDS, [
                species, run, "orthogroup", ambiguity, len(groups), len(values), len(retained),
                len(ambiguous), sum(not assignments.get(og) for og in values),
                statuses["no_valid_representative"], total, mapped, retained_tpm,
                retained_tpm / total if total else "", sum(map(len, contributions.values())),
                summed(summed_kos.values()),
            ])))
            sources.append({"run": run, "expression": file_record(path), "qc": file_record(qc_path)})
    write_tsv(out / "mapping_qc.tsv", QC_FIELDS, reports)
    write_json(out / "provenance.json", {
        "schema_version": 1, "annotation_scope": "orthogroup", "created_at": now(),
        "method": "median_length_one_representative", "ambiguity": ambiguity,
        "samples": file_record(samples), "annotation": annotation, "og_expression": sources,
        "quantification": "sum_original_OG_tpm_sum_without_renormalization",
        "unmapped_proteins": "not_annotated", "additional_representatives": False,
        "results": [file_record(out / name) for name in EXPRESSION_FILES],
    })


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    assign_parser = sub.add_parser("assign")
    for key in ["representatives", "annotation-dir", "outdir"]:
        assign_parser.add_argument(f"--{key}", required=True)
    aggregate_parser = sub.add_parser("aggregate")
    for key in ["samples", "annotation-dir", "og-run-dir", "outdir"]:
        aggregate_parser.add_argument(f"--{key}", required=True)
    aggregate_parser.add_argument("--ambiguity", choices=AMBIGUITY_POLICIES, default="duplicate")
    args = vars(parser.parse_args())
    action = args.pop("action")
    (assign if action == "assign" else aggregate)(**args)
