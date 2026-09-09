#!/usr/bin/env python3
"""FAMSA, trimAl and alignment QC, VeryFastTree, and ASTRAL-IV/CASTLES-II."""
import argparse
from collections import Counter
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

from busco_phylogeny import AMINO, fasta_records
from common import atomic_writer, file_record, now, read_tsv, sha256, write_json, write_tsv


def executable(command):
    resolved = shutil.which(command)
    if not resolved:
        raise ValueError(f"executable unavailable: {command}; see docs/phylogeny.md")
    return str(Path(resolved).resolve())


def read_tree(path, expected=None):
    from ete4 import Tree
    text = Path(path).read_text().strip()
    if not text or text.count(";") != 1:
        raise ValueError(f"expected exactly one Newick tree: {path}")
    tree = Tree(text, parser=1)
    leaves = list(tree.leaf_names())
    if len(leaves) < 4 or len(leaves) != len(set(leaves)):
        raise ValueError(f"tree requires at least four unique species: {path}")
    if expected is not None and set(leaves) != set(expected):
        raise ValueError(f"tree leaf set mismatch: {path}")
    for node in tree.traverse():
        if node.is_root:
            continue
        if node.dist is None or not math.isfinite(node.dist) or node.dist < 0:
            raise ValueError(f"invalid/missing branch length: {path}")
    return tree


def alignment_qc(records, settings):
    """After trimAl, drop short taxa and all-missing columns, without gap cutoffs."""
    import numpy as np
    names = [r[0] for r in records]
    if not records:
        return [], {"status": "too_few_taxa", "taxa": 0, "sites": 0, "variable_sites": 0, "informative_sites": 0,
                    "retained_columns_0based": []}
    if len(names) != len(set(names)) or len({len(r[1]) for r in records}) != 1:
        raise ValueError("alignment must have unique names and equal sequence lengths")
    matrix = np.frombuffer("".join(seq for _, seq in records).encode("ascii"),
                           dtype=np.uint8).reshape(len(records), len(records[0][1]))
    original_taxa, original_sites = matrix.shape
    known = np.isin(matrix, [ord(c) for c in AMINO])
    keep = known.sum(axis=1) >= settings["min_protein_length"]
    taxa = np.flatnonzero(keep)
    columns = np.flatnonzero(known[keep].any(axis=0))
    matrix = matrix[keep][:, columns]
    observed_states = np.zeros(len(columns), dtype=np.uint8)
    repeated_states = np.zeros(len(columns), dtype=np.uint8)
    for aa in AMINO:
        counts = (matrix == ord(aa)).sum(axis=0)
        observed_states += counts > 0
        repeated_states += counts >= 2
    variable = int((observed_states >= 2).sum())
    informative = int((repeated_states >= 2).sum())
    n, length = matrix.shape
    status = ("too_few_taxa" if n < settings["min_taxa"] else
              "no_variable_sites" if variable == 0 else "retained")
    output = [(names[i], row.tobytes().decode("ascii")) for i, row in zip(taxa, matrix)]
    retained_taxa = set(taxa.tolist())
    return output, {"status": status, "taxa_before": original_taxa, "sites_before": original_sites,
                    "taxa": n, "sites": length, "variable_sites": variable, "informative_sites": informative,
                    "retained_columns_0based": columns.tolist(),
                    "removed_species": [name for i, name in enumerate(names) if i not in retained_taxa]}


def align(fasta, output, qc, command, threads, settings):
    command = executable(command)
    inputs = list(fasta_records(fasta))
    names = [name for name, _ in inputs]
    if len(names) != len(set(names)):
        raise ValueError("duplicated species in marker FASTA")
    argv = None
    if len(names) >= settings["min_taxa"]:
        Path(output).parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix=".famsa-", dir=Path(output).parent) as tmp:
            raw = Path(tmp) / "alignment.faa"
            argv = [command, "-t", str(threads), "-v", str(Path(fasta).resolve()), str(raw)]
            subprocess.run(argv, check=True)
            records = list(fasta_records(raw))
            if len(records) != len(names) or set(n for n, _ in records) != set(names):
                raise ValueError("FAMSA changed the species set")
            if len({len(seq) for _, seq in records}) != 1:
                raise ValueError("FAMSA returned unequal sequence lengths")
            original = dict(inputs)
            if any(seq.replace("-", "") != original[name] for name, seq in records):
                raise ValueError("FAMSA changed input residues")
            report = {"status": "retained", "taxa": len(records), "sites": len(records[0][1])}
    else:
        records = []
        report = {"status": "too_few_taxa", "taxa": len(names), "sites": 0, "informative_sites": 0}
    with atomic_writer(output) as handle:
        if report["status"] == "retained":
            for name, sequence in records:
                handle.write(f">{name}\n{sequence}\n")
    write_json(qc, {**report, "created_at": now(), "input": file_record(fasta), "command": argv,
                    "executable": file_record(command), "settings": settings})


def trim(alignment, raw_qc, output, qc, columns, command, mode, settings):
    if mode not in {"gappyout", "automated1"}:
        raise ValueError("trimAl mode must be gappyout or automated1")
    command = executable(command)
    source_report = json.loads(Path(raw_qc).read_text())
    records = list(fasta_records(alignment))
    names = [name for name, _ in records]
    if len(names) != len(set(names)) or len({len(seq) for _, seq in records}) > 1:
        raise ValueError("alignment must have unique names and equal sequence lengths")
    if any(set(seq) - AMINO - set("X-") for _, seq in records):
        raise ValueError("unexpected protein alphabet before trimAl")
    if source_report["status"] == "retained" and (
            len(names) != source_report["taxa"] or not records or len(records[0][1]) != source_report["sites"]):
        raise ValueError("raw alignment disagrees with FAMSA QC")
    argv, selected, final_columns, retained = None, [], [], []
    report = {"status": source_report["status"], "taxa": len(records), "sites": 0,
              "variable_sites": 0, "informative_sites": 0}
    if source_report["status"] == "retained":
        if not any(c in AMINO for _, seq in records for c in seq):
            report.update(status="all_missing")
        else:
            Path(output).parent.mkdir(parents=True, exist_ok=True)
            with tempfile.TemporaryDirectory(prefix=".trimal-", dir=Path(output).parent) as tmp:
                gap_view, raw = Path(tmp) / "missing_as_gaps.faa", Path(tmp) / "trimmed.faa"
                # Masked X must count as missing in gap statistics. The real
                # alignment is projected by column number afterwards to retain X.
                viewed = {name: seq.replace("X", "-") for name, seq in records}
                with open(gap_view, "w") as handle:
                    for name, seq in viewed.items():
                        handle.write(f">{name}\n{seq}\n")
                argv = [command, "-in", str(gap_view), "-out", str(raw), "-fasta", "-keepheader",
                        "-keepseqs", "-colnumbering", "-" + mode]
                result = subprocess.run(argv, check=True, stdout=subprocess.PIPE, text=True)
                maps = re.findall(r"^#ColumnsMap\t(.*)$", result.stdout, re.MULTILINE)
                if len(maps) != 1 or not re.fullmatch(r"\d+(?:,\s*\d+)*\s*", maps[0]):
                    raise ValueError("trimAl did not return a valid column map")
                selected = [int(value) for value in maps[0].split(",")]
                if selected != sorted(set(selected)) or selected[-1] >= len(records[0][1]):
                    raise ValueError("trimAl column map is duplicated, unordered or out of range")
                trimmed = list(fasta_records(raw))
                if len(trimmed) != len(names) or {n for n, _ in trimmed} != set(names):
                    raise ValueError("trimAl changed the species set despite -keepseqs")
                if any(seq != "".join(viewed[name][i] for i in selected) for name, seq in trimmed):
                    raise ValueError("trimAl residues disagree with its column map")
                projected = [(name, "".join(seq[i] for i in selected)) for name, seq in records]
                retained, report = alignment_qc(projected, settings)
                final_columns = [selected[i] for i in report.pop("retained_columns_0based")]
    with atomic_writer(output) as handle:
        if report["status"] == "retained":
            for name, seq in retained:
                handle.write(f">{name}\n{seq}\n")
    write_tsv(columns, ["trimmed_column_1based", "famsa_column_1based"],
              [{"trimmed_column_1based": i + 1, "famsa_column_1based": col + 1}
               for i, col in enumerate(final_columns)] if report["status"] == "retained" else [])
    write_json(qc, {**report, "created_at": now(), "input": file_record(alignment),
                    "raw_alignment_qc": file_record(raw_qc), "command": argv,
                    "executable": file_record(command), "mode": mode, "settings": settings,
                    "unknown_as_gap_for_selection": True, "original_residues_restored": True,
                    "raw_sites": source_report["sites"], "trimal_sites": len(selected),
                    "trimal_columns_0based": selected, "columns": file_record(columns)})


def gene_tree(alignment, alignment_qc, output, qc, command, threads, seed):
    report = json.loads(Path(alignment_qc).read_text())
    argv = None
    tool_record = None
    with atomic_writer(output) as handle:
        if report["status"] == "retained":
            command = executable(command)
            # Topology search uses LG+CAT; -gamma rescales lengths with Gamma20.
            # Double precision avoids FastTree single-precision short-branch issues.
            argv = [command, "-threads", str(threads), "-double-precision", "-lg", "-gamma",
                    "-seed", str(seed), str(Path(alignment).resolve())]
            subprocess.run(argv, stdout=handle, check=True)
            tool_record = file_record(command)
    if report["status"] == "retained":
        read_tree(output, [name for name, _ in fasta_records(alignment)])
    write_json(qc, {"created_at": now(), "status": report["status"], "sites": report["sites"],
                   "taxa": report["taxa"], "alignment": file_record(alignment),
                   "command": argv, "executable": tool_record,
                   "support_type": "SH-like local support, not bootstrap proportions"})


def merge(manifest, markers, tree_dir, output, coverage, qc):
    species = read_tsv(manifest)
    expected = {r["species"] for r in species}
    counts = Counter({s: 0 for s in expected})
    retained, excluded = [], []
    with atomic_writer(output) as handle:
        for row in read_tsv(markers):
            marker = row["marker"]
            path = Path(tree_dir) / f"{marker}.nwk"
            report = json.loads((Path(tree_dir) / f"{marker}.json").read_text())
            if report["status"] != "retained":
                excluded.append({"marker": marker, "reason": report["status"]})
                continue
            tree = read_tree(path)
            names = set(tree.leaf_names())
            if names - expected:
                raise ValueError(f"unexpected species in {marker}")
            occupancy = len(names) / len(expected)
            counts.update(names)
            retained.append({"marker": marker, "sites": report["sites"], "tree": file_record(path),
                             "occupancy": occupancy})
            handle.write(path.read_text().strip() + "\n")
    write_tsv(coverage, ["species", "gene_trees", "represented"],
              [{"species": s, "gene_trees": n, "represented": n > 0} for s, n in sorted(counts.items())])
    missing = sorted(s for s, n in counts.items() if n == 0)
    count_distribution = Counter(counts.values())
    report = {"created_at": now(), "retained": retained, "excluded": excluded,
              "total_gene_sites": sum(r["sites"] for r in retained),
              "mean_gene_length": sum(r["sites"] for r in retained) / len(retained) if retained else 0,
              "species": len(expected),
              "gene_tree_count_distribution": [{"gene_trees": n, "species": count}
                                               for n, count in sorted(count_distribution.items())],
              "occupancy_definition": "retained gene-tree tips / all selected unique species",
              "status": "insufficient_coverage" if missing or not retained else "retained",
              "unrepresented_species": missing}
    write_json(qc, report)
    # Leave diagnostics as successful outputs; the inference stage refuses a
    # failed coverage check, so Snakemake does not delete the diagnostic tables.


def astral(trees, merge_qc, manifest, output, qc, command, outgroup, threads, seed):
    command = executable(command)
    # A filename alone cannot demonstrate a LARGE_DATA build.
    expected = {r["species"] for r in read_tsv(manifest)}
    if len(expected) > 5000 and Path(command).name != "astral4_int128":
        raise ValueError("more than 5000 species requires astral4_int128; run prepare_phylogeny_tools.py")
    if len(expected) > 5000:
        from prepare_phylogeny_tools import SOURCES
        build_path = Path(command).parent.parent / "aster.json"
        if not build_path.is_file():
            raise ValueError("int128 build provenance missing; run prepare_phylogeny_tools.py")
        build = json.loads(build_path.read_text())
        if (build.get("source") != SOURCES["aster"] or "LARGE_DATA" not in build.get("command", [])
                or build.get("executable", {}).get("sha256") != sha256(command)):
            raise ValueError("int128 build provenance/checksum mismatch; run prepare_phylogeny_tools.py")
    if outgroup not in expected:
        raise ValueError("ASTRAL outgroup is absent from the selected species")
    report = json.loads(Path(merge_qc).read_text())
    if report.get("status") != "retained":
        raise ValueError(f"insufficient gene-tree coverage: no retained trees or species absent from all retained trees; "
                         f"see {merge_qc} and species_coverage.tsv")
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".astral-", dir=Path(output).parent) as tmp:
        raw = Path(tmp) / "species_tree.nwk"
        argv = [command, "-i", str(Path(trees).resolve()), "-o", str(raw), "-t", str(threads),
                "--root", outgroup, "--seed", str(seed), "--length", "SULength",
                "--genelength", str(report["mean_gene_length"]), "-u", "1"]
        subprocess.run(argv, check=True)
        tree = read_tree(raw, expected)
        if len(tree.children) != 2 or not any(c.is_leaf and c.name == outgroup for c in tree.children):
            raise ValueError("ASTRAL did not return the requested rooted species tree")
        with atomic_writer(output) as handle:
            handle.write(raw.read_text().strip() + "\n")
    write_json(qc, {"created_at": now(), "command": argv, "executable": file_record(command),
                   "input": file_record(trees), "species": len(expected), "outgroup": outgroup,
                   "branch_length_unit": "substitutions_per_site", "branch_length_method": "CASTLES-II",
                   "support_type": "ASTRAL local posterior probability",
                   "mean_gene_length": report["mean_gene_length"],
                   "total_gene_sites": sum(r["sites"] for r in report["retained"]), "concatenation_used": False})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    fields = {
        "align": ["fasta", "output", "qc", "command", "settings"],
        "trim": ["alignment", "raw-qc", "output", "qc", "columns", "command", "mode", "settings"],
        "gene_tree": ["alignment", "alignment-qc", "output", "qc", "command"],
        "merge": ["manifest", "markers", "tree-dir", "output", "coverage", "qc"],
        "astral": ["trees", "merge-qc", "manifest", "output", "qc", "command", "outgroup"],
    }
    for action, names in fields.items():
        p = sub.add_parser(action)
        for name in names:
            p.add_argument("--" + name, required=True)
        if action in {"align", "gene_tree", "astral"}:
            p.add_argument("--threads", type=int, required=True)
        if action in {"gene_tree", "astral"}:
            p.add_argument("--seed", type=int, default=12345)
    args = vars(parser.parse_args())
    action = args.pop("action")
    if "settings" in args:
        args["settings"] = json.loads(args["settings"])
    globals()[action](**args)
