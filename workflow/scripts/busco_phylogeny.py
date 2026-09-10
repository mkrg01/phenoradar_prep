#!/usr/bin/env python3
"""Reuse BUSCO full tables and original CDS/proteins, with cdskit preparation."""
import argparse
from collections import Counter, OrderedDict, defaultdict
from functools import lru_cache
import gzip
import hashlib
import json
from pathlib import Path
import re

from common import atomic_writer, file_record, now, read_tsv, write_json, write_tsv


SAFE = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*\Z")
COORDINATES = re.compile(r":\d+-\d+(?:\|[+-])?\Z")
AMINO = set("ACDEFGHIKLMNPQRSTVWY")
CDSKIT_COMMIT = "218f6ed61abcac7f11dd81b17c087cb16119c39e"
CDSKIT_MODULES = {
    "pad": "66352003e366f54caa90ac2facb41366c889a87ea54c950fa309a77b8a4e7f30",
    "mask": "6af8ea346a96ebde8e47cb5d584a98f70aec0a2db022abd54217d6ca8c2908c6",
    "translate": "ea0164d9aa629a9afd0dfa36baffc4d6ce95df7a956b1fcd4ff702b327498e88",
}
CDSKIT_CONDA_MASK_SHA256 = "9f137703c4f2cefd0a7c49571a6cafb099c0dc5adc6238d4258cb0d529c31d3d"


@lru_cache(maxsize=1)
def cdskit_backend():
    """Use the pinned command implementations without starting a process per CDS."""
    import importlib
    import cdskit
    modules = {name: importlib.import_module("cdskit." + name) for name in CDSKIT_MODULES}
    records = {name: file_record(module.__file__) for name, module in modules.items()}
    expected = dict(CDSKIT_MODULES)
    if cdskit.__version__ == "0.27.0":
        expected["mask"] = CDSKIT_CONDA_MASK_SHA256
    if cdskit.__version__ not in {"0.27.0", "0.29.2"} or any(records[n]["sha256"] != digest for n, digest in expected.items()):
        raise ValueError("cdskit API differs from the tested source; install workflow/envs/phylogeny.yaml")
    return modules, {"version": cdskit.__version__,
                     "source": "bioconda::cdskit=0.27.0" if cdskit.__version__ == "0.27.0" else CDSKIT_COMMIT,
                     "modules": records,
                     "interface": "Python functions used by cdskit pad, mask and translate"}


def prepare_cds(sequence, marker, codon_table):
    modules, _ = cdskit_backend()
    normalized = sequence.replace("U", "T").replace("X", "N").replace(".", "-")
    if not normalized or set(normalized) - set("ACGTRYSWKMBDHVN-"):
        raise ValueError(f"empty/invalid CDS alphabet: {marker}")
    # genegalleon's padding heuristic can change the reading frame to reduce
    # internal stops. Record that decision explicitly; it is not an ORF proof.
    result = modules["pad"].process_record_padding(marker, normalized, codon_table, "N")
    padded = result["new_seq"]
    head, tail = 0, 0
    if result["log"]:
        match = re.search(r"head_padding=(\d+), tail_padding=(\d+)", result["log"])
        if not match:
            raise ValueError("unrecognized cdskit padding report")
        head, tail = map(int, match.groups())
    if padded != "N" * head + normalized + "N" * tail or len(padded) % 3:
        raise ValueError("cdskit padding changed original nucleotides or did not preserve complete codons")
    translate = modules["translate"].translate_sequence_string
    before_mask = translate(padded, codon_table, False)
    masked = modules["mask"].mask_sequence_string(padded, codon_table, "NNN", True, True)
    if len(masked) != len(padded):
        raise ValueError("cdskit mask changed sequence length")
    changed = [i // 3 for i in range(0, len(padded), 3) if padded[i:i+3] != masked[i:i+3]]
    if any(masked[3*i:3*i+3] != "NNN" for i in changed):
        raise ValueError("cdskit mask made a change other than codon masking")
    protein = translate(masked, codon_table, False)
    if len(protein) != len(masked) // 3 or "*" in protein:
        raise ValueError("cdskit masked-CDS translation is incomplete or contains stop codons")
    return protein, {
        "input_length_nt": len(sequence), "head_padding_nt": head, "tail_padding_nt": tail,
        "reading_frame_changed": bool(head % 3), "padding_report": result["log"].strip(),
        "internal_stops_after_padding": before_mask[:-1].count("*"),
        "terminal_stop_masked": before_mask.endswith("*"),
        "mask_changed_codons": len(changed), "mask_changed_codon_indices_0based": changed,
        "padded_cds_sha256": hashlib.sha256(padded.encode()).hexdigest(),
        "masked_cds_sha256": hashlib.sha256(masked.encode()).hexdigest(),
    }


def fasta_records(path):
    """Stream a FASTA; memory is bounded by one sequence."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        name, parts = None, []
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts).upper()
                fields = line[1:].split()
                if not fields:
                    raise ValueError(f"empty FASTA ID: {path}")
                name, parts = fields[0], []
            elif name is None:
                raise ValueError(f"sequence before FASTA header: {path}")
            else:
                parts.append(line)
        if name is not None:
            yield name, "".join(parts).upper()


def busco_table(path, expected_lineage):
    """Accept transcriptome/protein full tables; duplicated hits are never rescued."""
    hits = defaultdict(list)
    lineage = None
    header = False
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith("#"):
                match = re.search(r"lineage dataset is:\s*(\S+)", line)
                if match:
                    lineage = match[1]
                fields = line.lstrip("# ").strip().split("\t")
                if fields[:3] == ["Busco id", "Status", "Sequence"]:
                    if fields[3:5] != ["Score", "Length"]:
                        raise ValueError(f"genome BUSCO tables are unsupported; supply transcriptome/protein tables: {path}")
                    header = True
                continue
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2 or not SAFE.fullmatch(fields[0]):
                raise ValueError(f"invalid BUSCO row: {path}:{number}")
            marker, status = fields[:2]
            if status not in {"Complete", "Duplicated", "Fragmented", "Missing"}:
                raise ValueError(f"unknown BUSCO status {status}: {path}:{number}")
            if status != "Missing":
                if len(fields) < 5 or not fields[2]:
                    raise ValueError(f"incomplete BUSCO hit: {path}:{number}")
                score, length = float(fields[3]), int(fields[4])
                if not (0 <= score < float("inf")) or length <= 0:
                    raise ValueError(f"invalid BUSCO score/length: {path}:{number}")
                hits[marker].append((status, fields[2], score, length))
            else:
                hits[marker].append((status, "", 0, 0))
    if lineage != expected_lineage or not header or not hits:
        raise ValueError(f"BUSCO lineage/header mismatch: {path}; expected {expected_lineage}, found {lineage}")
    complete = {m: h[0] for m, h in hits.items() if len(h) == 1 and h[0][0] == "Complete"}
    # A single original sequence assigned to different markers is not independent evidence.
    originals = Counter(COORDINATES.sub("", h[1]) for h in complete.values())
    complete = {m: h for m, h in complete.items() if originals[COORDINATES.sub("", h[1])] == 1}
    return complete, set(hits)


def unique_species(samples):
    result = {}
    for row in read_tsv(samples):
        species = row["species"]
        if not SAFE.fullmatch(species):
            raise ValueError(f"unsafe species label: {species}")
        if species in result and row["cds"] != result[species]["cds"]:
            raise ValueError(f"conflicting CDS inputs for {species}")
        result[species] = row
    if len(result) < 4:
        raise ValueError("phylogeny requires at least four selected species")
    return dict(sorted(result.items()))


def plan(samples, outdir, settings):
    species = unique_species(samples)
    if settings["outgroup"] not in species:
        raise ValueError("phylogeny.outgroup must name a selected species; it is required for CASTLES-II branch lengths")
    counts, lengths = Counter(), Counter()
    manifest, universe = [], None
    for name, row in species.items():
        table = Path(settings["busco_full_dir"]) / (name + settings["busco_full_suffix"])
        complete, ids = busco_table(table, settings["lineage"])
        if universe is not None and ids != universe:
            raise ValueError(f"BUSCO marker set differs between species: {table}")
        universe = ids
        sequence = (Path(settings["sequence_dir"]) / (name + settings["sequence_suffix"])) if settings["sequence_dir"] else Path(row["cds"])
        if not sequence.is_file() or sequence.stat().st_size == 0:
            raise ValueError(f"missing sequence input: {sequence}")
        for marker, hit in complete.items():
            counts[marker] += 1
            lengths[marker] += hit[3]
        manifest.append({"species": name, "busco_table": str(table.resolve()),
                         "sequences": str(sequence.resolve()),
                         "single_copy_hits": len(complete), "busco_sha256": file_record(table)["sha256"]})
    stats = []
    for marker in sorted(universe):
        occupancy = counts[marker] / len(species)
        mean_length = lengths[marker] / counts[marker] if counts[marker] else 0
        eligible = counts[marker] >= settings["min_taxa"]
        stats.append({"marker": marker, "species": counts[marker], "occupancy": occupancy,
                      "mean_busco_length": mean_length,
                      "eligible": eligible, "selection_rank": "", "selected": False})
    # Coverage is the only biological ranking criterion. Break ties by ID for
    # reproducibility, independently of taxonomy, match length and input order.
    ranked = sorted((r for r in stats if r["eligible"]),
                    key=lambda r: (-r["occupancy"], r["marker"]))
    for rank, row in enumerate(ranked, 1):
        row["selection_rank"] = rank
    chosen = ranked[:settings["max_markers"]]
    if not chosen:
        raise ValueError("no BUSCO markers pass selection; inspect per-marker species counts, min_taxa and input lineage")
    for row in chosen:
        row["selected"] = True
    out = Path(outdir)
    write_tsv(out / "species.tsv", list(manifest[0]), manifest)
    write_tsv(out / "markers.tsv", list(stats[0]), chosen)
    write_tsv(out / "marker_stats.tsv", list(stats[0]), stats)
    write_json(out / "provenance.json", {"created_at": now(), "settings": settings,
               "samples": file_record(samples), "species": len(species),
               "markers": len(chosen), "eligible_markers": len(ranked),
               "occupancy_definition": "species with an admissible single-copy Complete hit / all selected unique species",
               "ranking": ["occupancy descending", "BUSCO marker ID ascending (lexicographic)"],
               "post_qc_backfill": False})


def extract(species, table, sequences, markers, output, qc, settings):
    complete, _ = busco_table(table, settings["lineage"])
    wanted = {r["marker"] for r in read_tsv(markers)}
    selected = {m: h for m, h in complete.items() if m in wanted}
    # These are already oriented in-frame CDS, not genomic intervals. MetaEuk's
    # coordinates omit exon/strand details in full_table.tsv: do NOT slice by them.
    aliases = defaultdict(set)
    for marker, hit in selected.items():
        aliases[hit[1]].add(marker)
        aliases[COORDINATES.sub("", hit[1])].add(marker)
    found, proteins, rejected = set(), {}, Counter()
    records = []
    backend = cdskit_backend()[1] if settings["sequence_mode"] == "cds" else None
    for gene, sequence in fasta_records(sequences):
        if gene not in aliases:
            continue
        for marker in aliases[gene]:
            if marker in found:
                raise ValueError(f"ambiguous or duplicated sequence ID for {species}/{marker}: {gene}")
            found.add(marker)
            if settings["sequence_mode"] == "cds":
                protein, preparation = prepare_cds(sequence, marker, settings["translation_table"])
            else:
                preparation, protein = {"source": "supplied protein"}, sequence.removesuffix("*")
            # Apply one missing-data alphabet to proteins from either source.
            known = sum(c in AMINO for c in protein)
            unknown_fraction = (len(protein) - known) / len(protein) if protein else 1.0
            reason = ""
            if "*" in protein:  # no CDS context is available to repair supplied proteins
                reason = "internal_stop"
            elif known < settings["min_protein_length"]:
                reason = "short_protein"
            elif unknown_fraction > settings["max_unknown_fraction"]:
                reason = "ambiguous_protein"
            if reason:
                rejected[reason] += 1
            else:
                proteins[marker] = "".join(c if c in AMINO else "X" for c in protein)
            records.append({"marker": marker, "gene": gene, "busco_sequence": selected[marker][1],
                            "protein_length": len(protein), "known_residues": known,
                            "unknown_fraction": unknown_fraction, "preparation": preparation,
                            "status": reason or "retained"})
    missing = sorted(set(selected) - found)
    if missing:
        raise ValueError(f"BUSCO IDs absent from original sequences for {species}: {missing[:10]}")
    with atomic_writer(output) as handle:
        for marker, protein in sorted(proteins.items()):
            handle.write(f">{marker}\n{protein}\n")
    write_json(qc, {"created_at": now(), "species": species, "table": file_record(table),
                   "sequences": file_record(sequences), "settings": settings,
                   "cdskit": backend,
                   "retained": len(proteins), "rejected": dict(rejected), "records": records})


def collect(manifest, markers, species_dir, outdir):
    """Transpose once, rather than scanning 5,000 proteomes once per marker."""
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    allowed = {r["marker"] for r in read_tsv(markers)}
    for marker in allowed:
        (out / f"{marker}.faa").write_text("")
    handles = OrderedDict()
    try:
        for row in read_tsv(manifest):
            seen = set()
            for marker, sequence in fasta_records(Path(species_dir) / f'{row["species"]}.faa'):
                if marker not in allowed or marker in seen:
                    raise ValueError(f"unexpected/duplicated marker for {row['species']}: {marker}")
                seen.add(marker)
                if marker not in handles:
                    if len(handles) >= 128:
                        handles.popitem(last=False)[1].close()
                    handles[marker] = open(out / f"{marker}.faa", "a")
                handles.move_to_end(marker)
                handles[marker].write(f'>{row["species"]}\n{sequence}\n')
    finally:
        for handle in handles.values():
            handle.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    p = sub.add_parser("plan")
    for name in ["samples", "outdir", "settings"]:
        p.add_argument("--" + name, required=True)
    p = sub.add_parser("extract")
    for name in ["species", "table", "sequences", "markers", "output", "qc", "settings"]:
        p.add_argument("--" + name, required=True)
    p = sub.add_parser("collect")
    for name in ["manifest", "markers", "species-dir", "outdir"]:
        p.add_argument("--" + name, required=True)
    args = vars(parser.parse_args())
    action = args.pop("action")
    if "settings" in args:
        args["settings"] = json.loads(args["settings"])
    globals()[action](**args)
