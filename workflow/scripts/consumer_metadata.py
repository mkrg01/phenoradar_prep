#!/usr/bin/env python3
"""Assemble the PhenoRadar input metadata from completed, explicit contrast branches."""
import argparse
import json
from pathlib import Path

from common import file_record, now, sha256, write_json, write_tsv
from phenoradar_metadata import fields, read_base, taxonomy_for, with_pairs


BRANCHES = {"all", "phenotyped", "representatives"}


def assemble(base, taxonomy, pair_tables, trait="C4"):
    if not isinstance(pair_tables, dict) or set(pair_tables) - BRANCHES:
        raise ValueError("contrast sources must name all, phenotyped, or representatives")
    rows = read_base(base, trait)
    ranks = taxonomy_for(taxonomy, {r["species"] for r in rows})
    for row in rows:
        for rank, value in ranks[row["species"]].items():
            if row.get(rank, "") and row[rank] != value:
                raise ValueError(f"base metadata {rank} differs from selected taxonomy: {row['species']}")
        row.update(ranks[row["species"]])
        row["contrast_pair_id"] = ""
    columns = fields(trait)
    # Start from canonical columns, so a former pair source never leaks through.
    rows = [{key: row[key] for key in columns} for row in rows]
    for branch, path in sorted(pair_tables.items()):
        joined = with_pairs(rows, path, trait)
        key = "contrast_pair_id" if len(pair_tables) == 1 else "contrast_pair_id_" + branch
        if key not in columns:
            columns.append(key)
        for row, paired in zip(rows, joined):
            row[key] = paired["contrast_pair_id"]
    return columns, rows


def create(base, taxonomy, output, pair_tables=None, trait="C4"):
    pair_tables = {} if pair_tables is None else pair_tables
    output = Path(output).resolve()
    root = output.parent.parent
    inputs = [Path(base).resolve(), Path(taxonomy).resolve(),
              *(Path(p).resolve() for p in pair_tables.values())]
    if output in inputs or output.with_suffix(".json") in inputs:
        raise ValueError("consumer metadata output must differ from its inputs")
    records = {p: file_record(p) for p in inputs}
    columns, rows = assemble(base, taxonomy, pair_tables, trait)
    for path, record in records.items():
        if sha256(path) != record["sha256"]:
            raise ValueError("metadata input changed during assembly: " + str(path))
    def relative(path):
        record = dict(records[Path(path).resolve()])
        record["path"] = str(Path(path).resolve().relative_to(root))
        return record
    sources = {"base": relative(base), "taxonomy": relative(taxonomy),
               "contrasts": {key: relative(path) for key, path in sorted(pair_tables.items())}}
    write_tsv(output, columns, rows)
    write_json(output.with_suffix(".json"), {
        "schema_version": 1, "created_at": now(), "trait": trait, "samples": len(rows),
        "pair_policy": "single_source" if len(pair_tables) == 1 else "separate_columns" if pair_tables else "no_contrasts",
        "sources": sources, "output": dict(file_record(output), path=str(output.relative_to(root))),
    })
    return rows


def validate(output, trait="C4"):
    """Check the join against its recorded sources without altering any file."""
    output = Path(output).resolve()
    root = output.parent.parent
    info = json.loads(output.with_suffix(".json").read_text())
    if info.get("schema_version") != 1 or info.get("trait") != trait:
        raise ValueError("consumer metadata provenance differs from the requested trait")
    def source(record):
        path = Path(record["path"])
        if path.is_absolute() or ".." in path.parts:
            raise ValueError("consumer metadata source must be inside the analysis")
        path = root / path
        if sha256(path) != record["sha256"]:
            raise ValueError("consumer metadata checksum mismatch: " + str(path))
        return path
    sources = info["sources"]
    policy = "single_source" if len(sources["contrasts"]) == 1 else "separate_columns" if sources["contrasts"] else "no_contrasts"
    if info.get("pair_policy") != policy:
        raise ValueError("consumer metadata pair policy differs from its sources")
    columns, expected = assemble(source(sources["base"]), source(sources["taxonomy"]),
                                 {branch: source(record) for branch, record in sources["contrasts"].items()}, trait)
    actual = read_base(output, trait)
    if (source(info["output"]) != output or not actual or list(actual[0]) != columns
            or actual != expected or info["samples"] != len(expected)):
        raise ValueError("consumer metadata differs from its source joins")
    return actual


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["base", "taxonomy", "output"]:
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--trait", default="C4")
    parser.add_argument("--pair-tables", type=json.loads, default={})
    create(**vars(parser.parse_args()))
