#!/usr/bin/env python3
"""Aggregate original transcript TPM for accepted KO assignments."""
import argparse
import csv
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

from common import file_record, now, write_json, write_tsv


HIT_FIELDS = ["species", "gene_id", "ko", "score", "threshold", "evalue", "assignment_status", "accepted"]
GENE_FIELDS = ["species", "gene_id", "assignment_status", "accepted_ko_count", "selected_ko", "terminal_stop_stripped"]
KO_FIELDS = ["species", "run", "ko", "tpm_sum", "annotated_genes", "quantified_genes"]
STATUSES = ("unique", "ambiguous", "below_threshold", "threshold_missing", "unannotated")
AMBIGUITY_POLICIES = ("duplicate", "drop", "error")
QC_FIELDS = ["species", "run", "ambiguity", "targets", "protein_genes", "quantified_proteins",
             "retained_targets", "annotated_kos", "quantified_kos", "total_tpm", "retained_tpm",
             "retained_tpm_fraction", "quantified_assignments", "ko_tpm_sum",
             "no_protein_targets", "no_protein_tpm", "no_retained_kos"]
QC_FIELDS += [f"{status}_{suffix}" for status in STATUSES
              for suffix in ("genes", "targets", "tpm")]


def table(path, fields):
    with open(path, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames or []
        if len(header) != len(set(header)) or not set(fields).issubset(header):
            raise ValueError(f"invalid TSV header: {path}; requires {fields}")
        rows = list(reader)
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise ValueError(f"malformed TSV row: {path}")
    return rows


def identifier(value, label):
    if not value or value in {".", ".."} or re.search(r"\s|[/\\]", value):
        raise ValueError(f"invalid {label}: {value!r}")
    return value


def number(value, label, nonnegative=True):
    try:
        result = float(value)
    except (ValueError, TypeError) as error:
        raise ValueError(f"invalid {label}: {value!r}") from error
    if not math.isfinite(result) or (nonnegative and result < 0):
        raise ValueError(f"invalid {label}: {value!r}")
    return result


def summed(values):
    try:
        value = math.fsum(values)
    except OverflowError as error:
        raise ValueError("TPM total is not finite") from error
    return number(value, "TPM total")


def selected_samples(samples):
    rows = table(samples, ["species", "odb_species", "run", "abundance"])
    runs, species, directories = set(), {}, {}
    for row in rows:
        for key in ("species", "odb_species", "run"):
            identifier(row[key], key)
        if row["run"] in runs:
            raise ValueError(f"duplicate run in manifest: {row['run']}")
        runs.add(row["run"])
        name, directory = row["species"], row["odb_species"]
        if name in species and species[name] != directory:
            raise ValueError(f"inconsistent odb_species for {name}")
        if directory in directories and directories[directory] != name:
            raise ValueError(f"annotation directory collision: {directory}")
        species[name], directories[directory] = directory, name
    if not rows:
        raise ValueError("empty sample manifest")
    return sorted(rows, key=lambda row: (row["species"], row["run"]))


def verify_file(path, record):
    current = file_record(path)
    if current["sha256"] != record.get("sha256") or current["bytes"] != record.get("bytes"):
        raise ValueError(f"file changed after completion: {path}")


def load_annotation(directory, species):
    """Validate completed annotation and the relation between hits and genes."""
    root = Path(directory)
    provenance = root / "provenance.json"
    report = json.loads(provenance.read_text())
    if report.get("identity", {}).get("species") != species:
        raise ValueError(f"annotation species does not match manifest: {species}")
    entries = report.get("results", [])
    records = {Path(entry["path"]).name: entry for entry in entries}
    if len(records) != len(entries) or not {"gene_kos.tsv", "genes.tsv"}.issubset(records):
        raise ValueError(f"incomplete annotation provenance: {provenance}")
    for name, record in records.items():
        verify_file(root / name, record)
    genes = {}
    gene_rows = table(root / "genes.tsv", GENE_FIELDS)
    for row in gene_rows:
        gene = identifier(row["gene_id"], "gene_id")
        if row["species"] != species:
            raise ValueError(f"gene species does not match manifest: {gene}")
        if gene in genes:
            raise ValueError(f"duplicate gene_id: {gene}")
        if row["assignment_status"] not in STATUSES:
            raise ValueError(f"unknown gene assignment_status: {gene}")
        if not re.fullmatch(r"0|[1-9]\d*", row["accepted_ko_count"]):
            raise ValueError(f"invalid accepted_ko_count: {gene}")
        if row["terminal_stop_stripped"] not in {"0", "1"}:
            raise ValueError(f"invalid terminal_stop_stripped: {gene}")
        genes[gene] = row
    hits = table(root / "gene_kos.tsv", HIT_FIELDS)
    pairs, accepted, hit_statuses = set(), defaultdict(set), defaultdict(set)
    for row in hits:
        gene, ko = row["gene_id"], row["ko"]
        if row["species"] != species or gene not in genes:
            raise ValueError(f"hit gene does not belong to species {species}: {gene}")
        if not re.fullmatch(r"K\d{5}", ko):
            raise ValueError(f"invalid KO identifier: {ko}")
        if (gene, ko) in pairs:
            raise ValueError(f"duplicate gene/KO hit: {gene}, {ko}")
        pairs.add((gene, ko))
        number(row["score"], "score", nonnegative=False)
        number(row["evalue"], "evalue")
        status = row["assignment_status"]
        if status not in {"accepted", "below_threshold", "threshold_missing"}:
            raise ValueError(f"invalid hit assignment_status: {gene}")
        if row["accepted"] != ("1" if status == "accepted" else "0"):
            raise ValueError(f"inconsistent accepted flag: {gene}")
        if status == "threshold_missing":
            if row["threshold"]:
                raise ValueError(f"threshold_missing hit has threshold: {gene}")
        else:
            number(row["threshold"], "threshold", nonnegative=False)
        # Kofam's marker precedes rounding in detail-tsv; do not recompute it.
        if status == "accepted":
            accepted[gene].add(ko)
        hit_statuses[gene].add(status)
    for gene, row in genes.items():
        kos = accepted[gene]
        expected = ("unique" if len(kos) == 1 else "ambiguous" if kos else
                    "threshold_missing" if "threshold_missing" in hit_statuses[gene] else
                    "below_threshold" if hit_statuses[gene] else "unannotated")
        selected = next(iter(kos)) if len(kos) == 1 else ""
        if (row["assignment_status"] != expected or int(row["accepted_ko_count"]) != len(kos)
                or row["selected_ko"] != selected):
            raise ValueError(f"gene summary disagrees with KO hits: {gene}")
    return genes, hits, file_record(provenance)


def abundance_values(path, run):
    values = {}
    for row in table(path, ["target_id", "tpm"]):
        gene = identifier(row["target_id"], "target_id")
        if gene in values:
            raise ValueError(f"duplicate target_id in {run}: {gene}")
        values[gene] = number(row["tpm"], f"TPM in {run}: {gene}")
    if not values:
        raise ValueError(f"{run}: abundance table has no targets")
    return values


def aggregate(samples, run, annotation_dir, output, qc, ambiguity="duplicate"):
    if ambiguity not in AMBIGUITY_POLICIES:
        raise ValueError("KO ambiguity must be duplicate, drop, or error; splitting is unsupported")
    manifest = selected_samples(samples)
    matches = [row for row in manifest if row["run"] == run]
    if len(matches) != 1:
        raise ValueError(f"expected exactly one manifest row for run {run}")
    sample = matches[0]
    genes, hits, source = load_annotation(annotation_dir, sample["species"])
    values = abundance_values(sample["abundance"], run)
    groups, contributions = defaultdict(list), defaultdict(list)
    status_counts = Counter(row["assignment_status"] for row in genes.values())
    target_counts, partitions = Counter(), defaultdict(list)
    retained_genes = set()
    # load_annotation verifies unique gene/KO pairs and the upstream acceptance marker.
    for hit in hits:
        gene, ko = hit["gene_id"], hit["ko"]
        if hit["accepted"] == "1" and (ambiguity == "duplicate" or genes[gene]["assignment_status"] == "unique"):
            groups[ko].append(gene)
            if gene in values:
                contributions[ko].append(values[gene])
                retained_genes.add(gene)
    for gene, value in values.items():
        row = genes.get(gene)
        status = row["assignment_status"] if row else "no_protein"
        target_counts[status] += 1
        partitions[status].append(value)
    if target_counts["ambiguous"] and ambiguity == "error":
        raise ValueError(f"{run}: {target_counts['ambiguous']} quantified genes have multiple accepted KOs")
    result = [{"species": sample["species"], "run": run, "ko": ko,
               "tpm_sum": summed(contributions[ko]) if contributions[ko] else "",
               "annotated_genes": len(members), "quantified_genes": len(contributions[ko])}
              for ko, members in sorted(groups.items())]
    total = summed(values.values())
    # Coverage counts each input gene once, even when its TPM supports several KOs.
    retained = summed(values[gene] for gene in retained_genes)
    report = {"created_at": now(), "species": sample["species"], "run": run, "ambiguity": ambiguity,
              "targets": len(values), "protein_genes": len(genes),
              "quantified_proteins": len(values) - target_counts["no_protein"],
              "retained_targets": len(retained_genes), "annotated_kos": len(result),
              "quantified_kos": sum(bool(contributions[ko]) for ko in groups),
              "total_tpm": total, "retained_tpm": retained,
              "retained_tpm_fraction": retained / total if total else None,
              "quantified_assignments": sum(len(items) for items in contributions.values()),
              "ko_tpm_sum": summed(row["tpm_sum"] for row in result if row["quantified_genes"]),
              "no_protein_targets": target_counts["no_protein"],
              "no_protein_tpm": summed(partitions["no_protein"]),
              "no_retained_kos": not any(contributions.values()),
              "abundance": file_record(sample["abundance"]), "annotation_provenance": source}
    for status in STATUSES:
        report.update({f"{status}_genes": status_counts[status], f"{status}_targets": target_counts[status],
                       f"{status}_tpm": summed(partitions[status])})
    write_tsv(output, KO_FIELDS, result)
    report["result"] = file_record(output)
    write_json(qc, report)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ["samples", "run", "annotation-dir", "output", "qc"]:
        parser.add_argument(f"--{flag}", required=True)
    parser.add_argument("--ambiguity", choices=AMBIGUITY_POLICIES, default="duplicate")
    aggregate(**vars(parser.parse_args()))
