#!/usr/bin/env python3
"""Prepare minimal species metadata for PhenoRadar from selected workflow inputs."""
import argparse
import csv
import re

from common import write_tsv
from species_traits import read_species_traits


def fields(trait="C4"):
    if (not isinstance(trait, str) or not trait.strip() or trait != trait.strip()
            or any(c in trait for c in "\t\r\n")
            or trait in {"species", "contrast_pair_id", "family"}):
        raise ValueError("trait must name a non-reserved output column")
    return ["species", trait, "contrast_pair_id", "family"]


def _table(path, required):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = reader.fieldnames or []
        if len(columns) != len(set(columns)) or not set(required) <= set(columns):
            raise ValueError(f"missing or duplicate TSV columns: {path}; required {required}")
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"malformed TSV row: {path}:{reader.line_num}")
            yield row


def _species(value):
    if not isinstance(value, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", value):
        raise ValueError(f"invalid species identifier: {value!r}")
    return value


def _trait(value, species):
    value = value.strip()
    value = "" if value in {"", "NA", "NaN", "nan"} else value
    if value not in {"", "0", "1"}:
        raise ValueError(f"PhenoRadar trait must be 0, 1 or missing: {species}={value}")
    return value


def _base_rows(rows, trait):
    required = fields(trait)
    result, seen = [], set()
    for row in rows:
        if not set(required) <= set(row) or any(row[key] is None for key in required):
            raise ValueError(f"base metadata requires {required}")
        name = _species(row["species"])
        if name in seen:
            raise ValueError("duplicate species in PhenoRadar metadata: " + name)
        seen.add(name)
        result.append({"species": name, trait: _trait(row[trait], name),
                       "contrast_pair_id": row["contrast_pair_id"].strip(),
                       "family": row["family"].strip()})
    if not result:
        raise ValueError("empty PhenoRadar metadata")
    return result


def read_base(path, trait="C4"):
    """Read validated minimal metadata without consulting original raw inputs."""
    return _base_rows(_table(path, fields(trait)), trait)


def with_pairs(rows, pairs, trait="C4"):
    """Left-join one contrast result, replacing pair IDs and preserving base order."""
    result = _base_rows(rows, trait)
    by_species = {row["species"]: row for row in result}
    for row in result:
        row["contrast_pair_id"] = ""
    seen = set()
    for pair in _table(pairs, ["species", trait, "contrast_pair_id"]):
        name = _species(pair["species"])
        if name in seen:
            raise ValueError("duplicate species in contrast metadata: " + name)
        seen.add(name)
        if name not in by_species:
            raise ValueError("contrast metadata contains unknown species: " + name)
        row = by_species[name]
        if _trait(pair[trait], name) != row[trait]:
            raise ValueError("contrast trait differs from base metadata: " + name)
        row["contrast_pair_id"] = pair["contrast_pair_id"].strip()
    return result


def prepare(samples, metadata, traits, output, trait="C4", pairs=None):
    """Create one species row, retaining missing traits and unpaired species.

    Multiple runs may describe a species when identity and taxonomy agree. Run
    aggregation is a separate expression-input decision, not metadata selection.
    """
    columns = fields(trait)
    selected, runs = {}, {}
    identities = ["scientific_name", "taxid", "odb_species", "cds"]
    for row in _table(samples, ["species", "run"]):
        name = _species(row["species"])
        run = row["run"]
        if not run or run != run.strip() or run in runs:
            raise ValueError("empty or duplicate run in selected samples: " + run)
        if name in selected and any(selected[name].get(key) != row.get(key) for key in identities):
            raise ValueError("conflicting sample identity for species: " + name)
        selected[name] = row
        runs[run] = name
    if not selected:
        raise ValueError("empty selected sample table")

    families, metadata_runs = {}, set()
    has_runs = False
    for row in _table(metadata, ["species", "family"]):
        name = _species(row["species"])
        if name not in selected:
            raise ValueError("selected metadata contains unknown species: " + name)
        family = row["family"].strip()
        if name in families and families[name] != family:
            raise ValueError("conflicting family for species: " + name)
        families[name] = family
        if "run" in row:
            has_runs = True
            run = row["run"]
            if run in metadata_runs or runs.get(run) != name:
                raise ValueError("selected metadata run/species differs from samples: " + run)
            metadata_runs.add(run)
    if set(families) != set(selected):
        raise ValueError("species missing from selected taxonomy metadata: "
                         + ", ".join(sorted(set(selected) - set(families))))
    if has_runs and metadata_runs != set(runs):
        raise ValueError("runs missing from selected taxonomy metadata: "
                         + ", ".join(sorted(set(runs) - metadata_runs)))

    annotations = read_species_traits(traits, trait) if traits else {}
    if not traits:
        print("PhenoRadar metadata: no phenotype source supplied; trait values are blank", flush=True)
    rows = [{"species": name, trait: _trait(annotations.get(name, ""), name),
             "contrast_pair_id": "", "family": families[name]} for name in sorted(selected)]
    if pairs is not None:
        rows = with_pairs(rows, pairs, trait)
    write_tsv(output, columns, rows)
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["samples", "metadata", "output"]:
        parser.add_argument(f"--{name}", required=True)
    parser.add_argument("--traits")
    parser.add_argument("--trait", default="C4")
    parser.add_argument("--pairs")
    prepare(**vars(parser.parse_args()))
