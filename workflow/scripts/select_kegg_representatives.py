#!/usr/bin/env python3
"""Select one median-length protein per observed OrthoDB group."""
import argparse
from array import array
from collections import Counter, defaultdict
import hashlib
from pathlib import Path
import re
import tempfile

from aggregate_ko_tpm import selected_samples
from busco_phylogeny import fasta_records
from common import file_record, now, species_from_gene_id, write_json, write_tsv
from dataset_assets import record, stat_identity
from mapping_tables import load_tables, read_species
from run_kofam import _publish

FIELDS = ["orthogroup", "representative_gene_id", "representative_species", "length",
          "median_length", "unknown_residues", "member_count", "valid_member_count",
          "invalid_member_count", "selection_status", "sequence_sha256", "batch"]
STANDARD = set("ACDEFGHIKLMNPQRSTVWY")
SAFE_OG = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*\Z")


def normalized(sequence):
    value = sequence[:-1] if sequence.endswith("*") else sequence
    return value if value and re.fullmatch(r"[A-Z]+", value) else None


def members(snapshot, name, protein):
    """Validate ownership/coverage while retaining only one sample's mapping."""
    genes, _ = read_species(snapshot, name)
    for gene, sequence in fasta_records(protein):
        if gene not in genes or species_from_gene_id(gene) != name:
            raise ValueError(f"duplicate or unmapped FASTA identifier: {gene}")
        groups = genes.pop(gene)
        if len(groups) != len(set(groups)) or len(groups) > 1:
            raise ValueError(f"multiple ODB assignments for {gene}")
        if not groups:
            continue
        og = groups[0]
        if not SAFE_OG.fullmatch(og):
            raise ValueError(f"invalid orthogroup: {og}")
        yield og, gene, normalized(sequence)
    if genes:
        raise ValueError(f"mapping genes missing from FASTA: {next(iter(genes))}")


def select(samples, mapping, protein_dir, outdir, batch_size=5000):
    if type(batch_size) is not int or batch_size < 1:
        raise ValueError("batch_size must be positive")
    manifest = selected_samples(samples)
    species = {row["species"]: row["odb_species"] for row in manifest}
    snapshot = load_tables(mapping, verify_files=False)
    if set(snapshot["tables"]) != set(species):
        raise ValueError("mapping species differ from selected samples")
    inputs = [record(samples), record(mapping)]
    lengths, counts = defaultdict(lambda: array("I")), Counter()
    sources = []
    for name, odb_name in sorted(species.items()):
        path = Path(protein_dir) / f"{odb_name}_protein.fa"
        source = record(path)
        if source["sha256"] != snapshot["tables"][name]["protein_sha256"]:
            raise ValueError(f"protein differs from ODB mapping: {name}")
        sources.append((name, path, source))
        for og, gene, sequence in members(mapping, name, path):
            counts[og] += 1
            if sequence is not None:
                lengths[og].append(len(sequence))
    # Packed integer lengths need about 4 bytes per mapped, valid protein.
    # Only the largest single OG is expanded into Python integers for sorting.
    targets, medians = {}, {}
    for og, values in lengths.items():
        ordered = sorted(values)
        n = len(ordered)
        targets[og] = ordered[n // 2]  # upper middle wins an even-median tie
        medians[og] = (ordered[(n - 1) // 2] + ordered[n // 2]) / 2
    valid_counts = {og: len(values) for og, values in lengths.items()}
    del lengths
    chosen = {}
    for name, path, source in sources:
        if source["stat"] != stat_identity(path):
            raise ValueError(f"protein changed during representative selection: {name}")
        for og, gene, sequence in members(mapping, name, path):
            if sequence is None or len(sequence) != targets.get(og):
                continue
            rank = (sum(aa not in STANDARD for aa in sequence), gene)
            if og not in chosen or rank < chosen[og][0]:
                chosen[og] = (rank, gene, name, sequence)
    for entry in [*inputs, *(s for _, _, s in sources)]:
        if entry["stat"] != stat_identity(entry["path"]):
            raise ValueError("inputs changed during representative selection")
    out = Path(outdir).absolute()
    out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".ko-representatives-", dir=out.parent) as tmp:
        stage = Path(tmp) / "result"
        (stage / "batches").mkdir(parents=True)
        batches, rows = [], []
        groups = sorted(chosen)
        batch_by_og = {}
        for start in range(0, len(groups), batch_size):
            batch = f"batch_{start // batch_size:05d}"
            path = stage / "batches" / f"{batch}.faa"
            subset = groups[start:start + batch_size]
            with path.open("w") as handle:
                for og in subset:
                    handle.write(f">{og}\n{chosen[og][3]}\n")
                    batch_by_og[og] = batch
            batches.append({"batch": batch, "orthogroups": subset})
        for og in sorted(counts):
            value = chosen.get(og)
            rows.append(dict(zip(FIELDS, [
                og, value[1] if value else "", value[2] if value else "",
                len(value[3]) if value else "", medians.get(og, ""),
                value[0][0] if value else "", counts[og], valid_counts.get(og, 0),
                counts[og] - valid_counts.get(og, 0),
                "selected" if value else "no_valid_representative",
                hashlib.sha256(value[3].encode()).hexdigest() if value else "",
                batch_by_og.get(og, ""),
            ])))
        write_tsv(stage / "representatives.tsv", FIELDS, rows)
        write_json(stage / "provenance.json", {
            "schema_version": 1, "annotation_scope": "orthogroup",
            "method": "median_length_one_representative",
            "tie_break": "nearest_median,longer,fewer_nonstandard_residues,gene_id",
            "created_at": now(), "odb_version": snapshot["version"], "odb_node": snapshot["node"],
            "samples": inputs[0], "mapping": inputs[1], "proteins": [s for _, _, s in sources],
            "orthogroups": len(counts), "representatives": len(chosen),
            "unmapped_proteins": sum(e["qc"]["protein_genes"] for e in snapshot["tables"].values()) - sum(counts.values()),
            "batch_size": batch_size, "batches": batches,
            "results": [{**file_record(p), "path": str(p.relative_to(stage))}
                        for p in sorted(stage.rglob("*")) if p.is_file()],
        })
        _publish(stage, out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ["samples", "mapping", "protein-dir", "outdir"]:
        parser.add_argument(f"--{flag}", required=True)
    parser.add_argument("--batch-size", type=int, default=5000)
    select(**vars(parser.parse_args()))
