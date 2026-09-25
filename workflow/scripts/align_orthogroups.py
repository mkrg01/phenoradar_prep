#!/usr/bin/env python3
"""Collect every mapped protein and save untrimmed, gene-labelled OG alignments."""
import argparse
from collections import OrderedDict, defaultdict
import json
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

from busco_phylogeny import fasta_records
from mapping_tables import load_tables, read_species
from common import atomic_writer, file_record, now, read_tsv, species_from_gene_id, write_json


SAFE_OG = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*\Z")
PROTEIN_ALPHABET = set("ABCDEFGHIJKLMNOPQRSTUVWXYZ*")


def validate_protein(gene, sequence):
    # Preserve ambiguity and stop symbols, including terminal stops. Reject
    # malformed/pre-gapped inputs rather than silently deleting their residues.
    if not sequence or set(sequence) - PROTEIN_ALPHABET:
        raise ValueError(f"empty or invalid ungapped protein: {gene}")


def collect(samples, mapping, protein_dir, outdir):
    """Read each species once, with at most 64 OG FASTA files open at a time."""
    species = {r["species"]: r["odb_species"] for r in read_tsv(samples)}
    outdir = Path(outdir)
    outdir.parent.mkdir(parents=True, exist_ok=True)
    sources = []
    data = load_tables(mapping, verify_files=False)
    if set(data['tables']) != set(species):
        raise ValueError('mapping species differ from selected species')
    groups = set()
    with tempfile.TemporaryDirectory(prefix=".og-inputs-", dir=outdir.parent) as temporary:
        staging = Path(temporary) / "inputs"
        staging.mkdir()
        handles = OrderedDict()
        try:
            for name, odb_name in sorted(species.items()):
                assignments, entry = read_species(mapping, name)
                assignments = {gene:ogs for gene,ogs in assignments.items() if ogs}
                for ogs in assignments.values():
                    if any(not SAFE_OG.fullmatch(og) for og in ogs):
                        raise ValueError('orthogroup IDs must be safe, nonempty filename components')
                    groups.update(ogs)
                source = Path(protein_dir) / f"{odb_name}_protein.fa"
                seen = set()
                for gene, sequence in fasta_records(source):
                    if gene in seen:
                        raise ValueError(f"duplicate protein ID: {gene}")
                    seen.add(gene)
                    gene_groups = assignments.pop(gene, [])
                    if not gene_groups:
                        continue
                    if species_from_gene_id(gene) != name:
                        raise ValueError(f"gene ID species differs from mapping/protein species: {gene}: {name}")
                    validate_protein(gene, sequence)
                    for og in gene_groups:
                        handle = handles.pop(og, None)
                        if handle is None:
                            if len(handles) >= 64:
                                handles.popitem(last=False)[1].close()
                            handle = open(staging / f"{og}.faa", "a", encoding="utf-8")
                        handles[og] = handle
                        handle.write(f">{gene}\n{sequence}\n")
                if assignments:
                    raise ValueError(f"mapped protein missing from {source}: {next(iter(assignments))}")
                sources.append(file_record(source))
        finally:
            for handle in handles.values():
                handle.close()
        write_json(staging / "provenance.json", {
            "created_at": now(), "orthogroups": sorted(groups),
            "samples": file_record(samples), "mapping": file_record(mapping), "proteins": sources,
        })
        # This is a workflow-owned checkpoint directory; replace it only
        # after collection succeeds, so removed OGs cannot survive a rerun.
        if outdir.exists():
            shutil.rmtree(outdir)
        staging.replace(outdir)

def align(fasta, output, provenance, threads=1, command="famsa"):
    inputs = list(fasta_records(fasta))
    original = dict(inputs)
    if not inputs or len(original) != len(inputs):
        raise ValueError("OG FASTA must contain unique, nonempty gene IDs")
    for gene, sequence in inputs:
        species_from_gene_id(gene)
        validate_protein(gene, sequence)
    if type(threads) is not int or threads < 1:
        raise ValueError("alignment threads must be positive")
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    argv, executable = None, None
    if len(inputs) == 1:
        records = original
    else:
        resolved = shutil.which(command)
        if resolved is None:
            raise FileNotFoundError(f"alignment executable not found: {command}")
        executable = file_record(resolved)
        with tempfile.TemporaryDirectory(prefix=".famsa-", dir=output.parent) as temporary:
            aligned = Path(temporary) / "alignment.faa"
            argv = [str(Path(resolved).resolve()), "-t", str(threads), "-v", str(Path(fasta).resolve()), str(aligned)]
            subprocess.run(argv, check=True)
            result = list(fasta_records(aligned))
            records = dict(result)
            if len(result) != len(inputs) or set(records) != set(original):
                raise ValueError("FAMSA changed the gene ID set")
            if len({len(seq) for seq in records.values()}) != 1:
                raise ValueError("FAMSA returned unequal sequence lengths")
            if any(records[gene].replace("-", "") != seq for gene, seq in inputs):
                raise ValueError("FAMSA changed input residues")
    with atomic_writer(output) as handle:
        for gene, _ in inputs:
            handle.write(f">{gene}\n{records[gene]}\n")
    write_json(provenance, {
        "created_at": now(), "input": file_record(fasta), "alignment": file_record(output),
        "sequences": len(inputs), "sites": len(next(iter(records.values()))),
        "method": "singleton" if len(inputs) == 1 else "famsa",
        "command": argv, "executable": executable, "threads": threads,
    })


def finish(inputs, outdir, reports):
    """Publish provenance after all current OG jobs have succeeded."""
    inputs, outdir, reports = Path(inputs), Path(outdir), Path(reports)
    plan = json.loads((inputs / "provenance.json").read_text())
    groups = plan["orthogroups"]
    alignments = []
    for og in groups:
        report_path = reports / f"{og}.json"
        report = json.loads(report_path.read_text())
        expected = outdir / f"{og}.faa"
        if not expected.is_file() or Path(report["alignment"]["path"]).name != expected.name:
            raise ValueError(f"missing or misplaced OG alignment: {og}")
        current = file_record(expected)
        if any(current[key] != report["alignment"][key] for key in ["bytes", "sha256"]):
            raise ValueError(f"OG alignment changed after its job completed: {og}")
        alignments.append({"orthogroup": og, "alignment": current,
                           "provenance": file_record(report_path)})
    outdir.mkdir(parents=True, exist_ok=True)
    # This output directory owns its *.faa files. Prune obsolete OGs only once
    # every current alignment exists; unrelated files are left alone.
    expected_names = {f"{og}.faa" for og in groups}
    for path in outdir.glob("*.faa"):
        if path.name not in expected_names:
            path.unlink()
    write_json(outdir / "provenance.json", {
        "created_at": now(), "collection": file_record(inputs / "provenance.json"),
        "alignments": alignments,
    })


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="action", required=True)
    for action, flags in {
        "collect": ["samples", "mapping", "protein-dir", "outdir"],
        "align": ["fasta", "output", "provenance"],
        "finish": ["inputs", "outdir", "reports"],
    }.items():
        sub = commands.add_parser(action)
        for flag in flags:
            sub.add_argument(f"--{flag}", required=True)
        if action == "align":
            sub.add_argument("--threads", type=int, default=1)
    args = vars(parser.parse_args())
    globals()[args.pop("action")](**args)
