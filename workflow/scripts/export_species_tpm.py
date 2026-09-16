#!/usr/bin/env python3
"""Export species/orthogroup/tpm for PhenoRadar from a run-level TPM table."""
import argparse
from collections import Counter
import csv
import math
from pathlib import Path

from common import read_tsv, write_tsv


def export(samples, input, output):
    """Require one run per species and preserve the original numeric text.

    Input must have one contiguous block per run, as written by merge_tpm.
    Memory is bounded to one run's feature identifiers plus the manifest.
    """
    if Path(input).resolve() == Path(output).resolve():
        raise ValueError("input and output must be different files")
    runs = {}
    for row in read_tsv(samples):
        if not row.get("species") or not row.get("run") or row["run"] in runs:
            raise ValueError("empty species/run or duplicate run in sample manifest")
        runs[row["run"]] = row["species"]
    if not runs:
        raise ValueError("empty sample manifest")
    counts = Counter(runs.values())
    repeated = sorted(species for species, count in counts.items() if count > 1)
    if repeated:
        raise ValueError("multiple runs per species; select exactly one run per species before "
                         "PhenoRadar export (no automatic aggregation): " + ", ".join(repeated[:10]))

    def rows():
        seen, features = set(), set()
        current = None
        with open(input, newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = reader.fieldnames or []
            if (len(fields) != len(set(fields))
                    or not {"species", "run", "orthogroup", "tpm"} <= set(fields)):
                raise ValueError("input requires species, run, orthogroup, tpm columns")
            for row in reader:
                if None in row or any(value is None for value in row.values()):
                    raise ValueError(f"malformed TPM row at line {reader.line_num}")
                species, run, og = row["species"], row["run"], row["orthogroup"]
                if runs.get(run) != species:
                    raise ValueError(f"TPM run/species differs from manifest: {species}/{run}")
                if run != current:
                    if run in seen:
                        raise ValueError(f"duplicate/discontiguous run block: {run}")
                    seen.add(run)
                    current, features = run, set()
                if not og.strip() or og in features:
                    raise ValueError(f"empty or duplicate run/orthogroup coordinate: {run}/{og}")
                features.add(og)
                try:
                    number = float(row["tpm"])
                except ValueError as error:
                    raise ValueError(f"invalid TPM: {run}/{og}: {row['tpm']!r}") from error
                if not math.isfinite(number) or number < 0:
                    raise ValueError(f"invalid TPM: {run}/{og}: {row['tpm']!r}")
                yield dict(species=species, orthogroup=og, tpm=row["tpm"])
        if seen != set(runs):
            raise ValueError(f"TPM lacks selected runs: {sorted(set(runs) - seen)}")

    write_tsv(output, ["species", "orthogroup", "tpm"], rows())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["samples", "input", "output"]:
        parser.add_argument(f"--{name}", required=True)
    export(**vars(parser.parse_args()))
