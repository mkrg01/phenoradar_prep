"""Validate and collect binary traits for the selected species."""
import argparse
from pathlib import Path

from common import file_record, now, read_tsv, write_json, write_tsv
from sample_identity import species_id
from species_traits import read_species_traits


def normalize(samples, output, report, trait="C4", traits=None):
    destinations = {Path(output).resolve(), Path(report).resolve()}
    if len(destinations) != 2 or any(Path(path).resolve() in destinations for path in [samples, traits] if path):
        raise ValueError("trait normalization outputs must be distinct from each other and their source inputs")
    taxa = sorted({species_id(row) for row in read_tsv(samples)})
    if not taxa:
        raise ValueError("trait normalization requires selected samples")
    annotations = read_species_traits(traits, trait) if traits else {}
    observed = {annotations.get(name, "") for name in taxa} - {""}
    invalid = observed - {"0", "1"}
    if invalid:
        raise ValueError(f"{trait} must be 0, 1 or missing; invalid states for selected species: {sorted(invalid)}")
    rows = [{"species": name, trait: annotations.get(name, "")} for name in taxa]
    write_tsv(output, ["species", trait], rows)
    write_json(report, {"created_at": now(), "trait": trait,
                       "observed_states": sorted(observed), "biological_species": len(taxa),
                       "missing_trait_species": sum(row[trait] == "" for row in rows),
                       "samples": file_record(samples), "source": file_record(traits) if traits else None,
                       "normalized": file_record(output)})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["samples", "output", "report", "trait"]:
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--traits")
    normalize(**vars(parser.parse_args()))
