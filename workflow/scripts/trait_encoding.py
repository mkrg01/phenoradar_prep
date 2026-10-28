"""Normalize selected binary traits once, retaining explicit source labels."""
import argparse
import json
from pathlib import Path

from common import file_record, now, read_tsv, write_json, write_tsv
from sample_identity import species_id
from species_traits import read_species_traits


MISSING = {"", "NA", "NaN", "nan"}


def validate_mapping(mapping):
    if mapping is None:
        return None
    if (not isinstance(mapping, dict) or len(mapping) != 2
            or any(not isinstance(label, str) or label != label.strip() or label in MISSING
                   for label in mapping)
            or any(type(value) is not int for value in mapping.values())
            or set(mapping.values()) != {0, 1}):
        raise ValueError("trait_encoding must map two distinct non-missing string labels to the integers 0 and 1")
    return mapping


def read_labels(path, trait):
    if path is None:
        return {}
    report = json.loads(Path(path).read_text())
    mapping = validate_mapping(report.get("mapping"))
    if report.get("trait") != trait or mapping is None:
        raise ValueError("trait encoding record differs from the requested trait")
    return {str(value): label for label, value in mapping.items()}


def normalize(samples, output, report, trait="C4", traits=None, mapping=None):
    """Unknown values stay missing; never guess which text label is positive."""
    mapping = validate_mapping(mapping)
    destinations = {Path(output).resolve(), Path(report).resolve()}
    if len(destinations) != 2 or any(Path(path).resolve() in destinations for path in [samples, traits] if path):
        raise ValueError("trait normalization outputs must be distinct from each other and their source inputs")
    taxa = sorted({species_id(row) for row in read_tsv(samples)})
    if not taxa:
        raise ValueError("trait normalization requires selected samples")
    annotations = read_species_traits(traits, trait) if traits else {}
    observed = sorted({annotations.get(name, "") for name in taxa} - {""})
    if mapping is None:
        if not set(observed) <= {"0", "1"}:
            raise ValueError(f"trait_encoding is required for text states of {trait}: {observed}; "
                             "map the negative/reference state to 0 and the positive/focal state to 1")
        mapping = {"0": 0, "1": 1}
    unknown = set(observed) - set(mapping)
    if unknown:
        raise ValueError(f"selected {trait} states are absent from trait_encoding: {sorted(unknown)}")
    rows = [{"species": name, trait: str(mapping[annotations[name]]) if annotations.get(name, "") else ""}
            for name in taxa]
    write_tsv(output, ["species", trait], rows)
    write_json(report, {"created_at": now(), "trait": trait, "mapping": mapping,
                       "observed_source_states": observed, "biological_species": len(taxa),
                       "missing_trait_species": sum(row[trait] == "" for row in rows),
                       "samples": file_record(samples), "source": file_record(traits) if traits else None,
                       "normalized": file_record(output)})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["samples", "output", "report", "trait"]:
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--traits")
    parser.add_argument("--mapping", type=json.loads)
    normalize(**vars(parser.parse_args()))
