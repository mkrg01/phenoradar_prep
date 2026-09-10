#!/usr/bin/env python3
"""Make deterministic ODB chunks from the selected sample manifest."""
import argparse
from pathlib import Path

from common import atomic_writer, read_tsv, write_json


def chunk_plan(rows, protein_dir, chunk_size):
    if chunk_size < 1:
        raise ValueError("chunk size must be positive")
    species = {}
    for row in rows:
        species[row["species"]] = row["odb_species"]
    if not species:
        raise ValueError("empty sample manifest")
    chunks = []
    ordered = sorted(species)
    for start in range(0, len(ordered), chunk_size):
        label = f"chunk_{start // chunk_size:03d}"
        members = ordered[start:start + chunk_size]
        proteins = [str((Path(protein_dir) / f"{species[s]}_protein.fa").resolve()) for s in members]
        chunks.append({"chunk": label, "species": members, "proteins": proteins})
    return chunks


def make(samples, protein_dir, outdir, chunk_size):
    chunks = chunk_plan(read_tsv(samples), protein_dir, chunk_size)
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    master = []
    for chunk in chunks:
        with atomic_writer(out / f'{chunk["chunk"]}.fs') as handle:
            handle.write("\n".join(chunk["proteins"]) + "\n")
        master.extend(chunk["proteins"])
    with atomic_writer(out / "master.fs") as handle:
        handle.write("\n".join(master) + "\n")
    write_json(out / "chunks.json", chunks)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ["samples", "protein-dir", "outdir"]:
        parser.add_argument(f"--{flag}", required=True)
    parser.add_argument("--chunk-size", type=int, default=50)
    args = vars(parser.parse_args())
    make(**args)
