#!/usr/bin/env python3
"""Annotate one species with KofamScan while retaining uncertain and absent calls."""
import argparse
from collections import Counter
import csv
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

from common import file_record, now, sha256, write_json, write_tsv
from prepare_kegg_reference import read_ko_list
from verify_kegg_reference import verify


HIT_FIELDS = ["species", "gene_id", "ko", "score", "threshold", "evalue",
              "assignment_status", "accepted"]
GENE_FIELDS = ["species", "gene_id", "assignment_status", "accepted_ko_count",
               "selected_ko", "terminal_stop_stripped"]
RESULT_NAMES = ["detail.tsv", "gene_kos.tsv", "genes.tsv", "stdout.log", "stderr.log",
                "execution_config.json"]


def normalize_protein(source, destination):
    """Preserve FASTA headers/IDs and remove exactly one terminal translation stop."""
    records, seen, header, sequence = [], set(), None, []

    def finish(handle):
        if header is None:
            return
        gene = header[1:].split()[0]
        if gene in seen:
            raise ValueError(f"duplicate protein FASTA identifier: {gene}")
        seen.add(gene)
        value = "".join(sequence)
        stripped = value.endswith("*")
        if stripped:
            value = value[:-1]
        if "*" in value:
            raise ValueError(f"internal stop in protein FASTA: {gene}")
        if not value or not re.fullmatch(r"[A-Za-z]+", value):
            raise ValueError(f"empty or invalid amino acid sequence: {gene}")
        handle.write(header + "\n")
        for start in range(0, len(value), 80):
            handle.write(value[start:start + 80] + "\n")
        records.append({"gene_id": gene, "terminal_stop_stripped": int(stripped)})

    with open(source, encoding="utf-8") as inp, open(destination, "w", encoding="utf-8") as out:
        for line in inp:
            line = line.rstrip("\r\n")
            if not line.strip():
                continue
            if line.startswith(">"):
                if not line[1:] or line[1].isspace() or "\x00" in line:
                    raise ValueError("empty or invalid protein FASTA header")
                finish(out)
                header, sequence = line, []
            elif header is None:
                raise ValueError("protein FASTA sequence precedes its header")
            else:
                sequence.append("".join(line.split()))
        finish(out)
    if not records:
        raise ValueError("protein FASTA contains no sequences")
    return records


def _number(value, field, line, nonnegative=False):
    try:
        number = float(value)
    except ValueError as error:
        raise ValueError(f"invalid {field} at KofamScan line {line}: {value!r}") from error
    if not math.isfinite(number) or (nonnegative and number < 0):
        raise ValueError(f"invalid {field} at KofamScan line {line}: {value!r}")
    return number


def parse_detail(path, species, proteins, ko_entries, profile_kos=None):
    """Read upstream seven-column detail-tsv, retaining its pre-rounding decision.

    Upstream prints scores to one decimal and thresholds to two decimals AFTER
    calculating the '*' decision. Recomputing score >= threshold from that text
    would change near-threshold calls. Check consistency within rounding bounds.
    """
    genes = {record["gene_id"]: record for record in proteins}
    pairs, no_hits, header_seen = {}, set(), False
    with open(path, encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t", strict=True)
        for row in reader:
            line = reader.line_num
            if not row:
                continue
            if row[0] == "#":
                if row == ["#", "gene name", "KO", "thrshld", "score", "E-value", "KO definition"]:
                    if header_seen:
                        raise ValueError("duplicate KofamScan detail-tsv header")
                    header_seen = True
                elif not header_seen or not all(re.fullmatch(r"-+", value) for value in row[1:]):
                    raise ValueError(f"unexpected KofamScan header at line {line}")
                continue
            if not header_seen or len(row) != 7:
                raise ValueError(f"expected seven-column KofamScan detail-tsv at line {line}")
            marker, gene, ko, threshold, score, evalue, _ = row
            if marker not in {"", "*"}:
                raise ValueError(f"invalid KofamScan acceptance marker at line {line}")
            if gene not in genes:
                raise ValueError(f"KofamScan gene does not belong to input species {species}: {gene}")
            if not ko:
                if marker or any(row[2:]):
                    raise ValueError(f"malformed unannotated row at line {line}")
                no_hits.add(gene)
                continue
            if ko not in ko_entries or (profile_kos is not None and ko not in profile_kos):
                raise ValueError(f"KofamScan KO is absent from reference profiles: {ko}")
            numeric_score = _number(score, "score", line)
            _number(evalue, "E-value", line, nonnegative=True)
            expected = ko_entries[ko]["threshold"]
            missing = threshold in {"", "-"}
            if missing != (expected == "-"):
                raise ValueError(f"KofamScan threshold differs from reference for {ko}")
            if missing:
                if marker:
                    raise ValueError(f"accepted KO has no threshold: {ko}")
                threshold = ""
                status = "threshold_missing"
            else:
                numeric_threshold = _number(threshold, "threshold", line)
                if abs(numeric_threshold - float(expected)) > 0.00500001:
                    raise ValueError(f"KofamScan threshold differs from reference for {ko}")
                delta = numeric_score - float(expected)
                if (marker and delta < -0.05000001) or (not marker and delta > 0.05000001):
                    raise ValueError(f"KofamScan acceptance marker disagrees with threshold for {gene}/{ko}")
                status = "accepted" if marker else "below_threshold"
            result = dict(zip(HIT_FIELDS, [species, gene, ko, score, threshold, evalue,
                                          status, int(bool(marker))]))
            key = (gene, ko)
            if key in pairs and result != pairs[key]:
                raise ValueError(f"conflicting duplicate KofamScan hit: {gene}/{ko}")
            pairs[key] = result
    if not header_seen:
        raise ValueError("missing KofamScan detail-tsv header")
    hit_genes = {gene for gene, _ in pairs}
    if no_hits & hit_genes:
        raise ValueError("KofamScan gene has both unannotated and hit rows")
    by_gene = {gene: [] for gene in genes}
    hits = [pairs[key] for key in sorted(pairs)]
    for hit in hits:
        by_gene[hit["gene_id"]].append(hit)
    rows = []
    for protein in proteins:
        gene = protein["gene_id"]
        candidates = by_gene[gene]
        accepted = sorted({hit["ko"] for hit in candidates if hit["accepted"]})
        if len(accepted) == 1:
            status = "unique"
        elif accepted:
            status = "ambiguous"
        elif any(hit["assignment_status"] == "threshold_missing" for hit in candidates):
            status = "threshold_missing"
        elif candidates:
            status = "below_threshold"
        else:
            status = "unannotated"
        rows.append({"species": species, **protein, "assignment_status": status,
                     "accepted_ko_count": len(accepted),
                     "selected_ko": accepted[0] if status == "unique" else ""})
    return hits, rows


def _software(command):
    executable = shutil.which(str(command))
    if executable is None:
        raise FileNotFoundError(f"KofamScan executable not found: {command}")
    executable = Path(executable).resolve(strict=True)
    records = {"command": file_record(executable)}
    # KofamScan is a Ruby entry point: its adjacent library is part of the tool.
    library = executable.parent / "lib" / "kofam_scan"
    if library.is_dir():
        files = sorted(library.rglob("*.rb"))
        top = library.with_suffix(".rb")
        if top.is_file():
            files.append(top)
        records["kofam_library"] = [file_record(path) for path in files]
    for dependency in ["ruby", "hmmsearch", "parallel"]:
        path = shutil.which(dependency)
        records[dependency] = file_record(path) if path else None
    return executable, records


def _completed(output, fingerprint):
    try:
        record = json.loads((output / "provenance.json").read_text(encoding="utf-8"))
        if record.get("schema_version") != 1 or record.get("fingerprint") != fingerprint:
            return None
        results = record.get("results", [])
        if (len(results) != len(RESULT_NAMES)
                or {Path(item["path"]).name for item in results} != set(RESULT_NAMES)):
            return None
        for item in results:
            path = output / Path(item["path"]).name
            if path.is_symlink() or Path(item["path"]) != path or file_record(path) != item:
                return None
        return record
    except (OSError, ValueError, KeyError, TypeError):
        return None


def _publish(source, destination):
    """Replace a generated directory only after all files and provenance exist."""
    staging = Path(tempfile.mkdtemp(prefix=f".{destination.name}.publish-", dir=destination.parent))
    try:
        shutil.copytree(source, staging / "new")
        if destination.exists():
            destination.rename(staging / "old")
        try:
            (staging / "new").rename(destination)
        except BaseException:
            if (staging / "old").exists():
                (staging / "old").rename(destination)
            raise
    finally:
        shutil.rmtree(staging)


def run(protein, species, reference, output_dir, work_dir, command="exec_annotation", threads=1):
    if not species or any(character.isspace() for character in species) or "\x00" in species:
        raise ValueError("species must be a nonempty identifier without whitespace")
    if threads < 1 or threads > int(os.environ.get("SLURM_CPUS_PER_TASK", threads)):
        raise ValueError("invalid KofamScan thread count or exceeds Slurm CPU allocation")
    protein = Path(protein).resolve(strict=True)
    reference = Path(reference).resolve(strict=True)
    ref = verify(reference, full=False)
    executable, software = _software(command)
    scripts = Path(__file__).resolve().parent
    identity = {"species": species, "protein": file_record(protein),
                "reference": file_record(reference), "reference_id": ref["reference_id"],
                "reference_inventory": file_record(ref["inventory_path"]),
                "software": software,
                "options": {"threads": threads, "format": "detail-tsv", "threshold_scale": 1,
                            "report_unannotated": False, "terminal_stop": "strip_one_reject_internal"},
                "implementation": [file_record(scripts / name) for name in
                                   ["run_kofam.py", "common.py", "prepare_kegg_reference.py",
                                    "verify_kegg_reference.py"]]}
    fingerprint = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    out, base = Path(output_dir).absolute(), Path(work_dir).absolute()
    if out.is_symlink() or base.is_symlink() or (out.exists() and not out.is_dir()):
        raise ValueError("KofamScan output/work must be directories without a symlink leaf")
    out.parent.mkdir(parents=True, exist_ok=True)
    base.mkdir(parents=True, exist_ok=True)
    out, base = out.parent.resolve() / out.name, base.resolve()
    if out == base or out.is_relative_to(base) or base.is_relative_to(out):
        raise ValueError("KofamScan work and output directories must not overlap")
    with open(base / ".lock", "a") as lock, open(out.parent / f".{out.name}.lock", "a") as out_lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        fcntl.flock(out_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        status_path = base / "status.json"
        cached = _completed(out, fingerprint)
        if cached is not None:
            write_json(status_path, {"state": "success", "completed_at": now(),
                                     "fingerprint": fingerprint, "reused": True})
            return cached
        work = Path(tempfile.mkdtemp(prefix=f"{fingerprint}.", dir=base))
        result = work / "result"
        result.mkdir()
        write_json(work / "identity.json", identity)
        write_json(status_path, {"state": "running", "started_at": now(),
                                 "fingerprint": fingerprint, "work": str(work)})
        try:
            normalized = work / "protein.faa"
            proteins = normalize_protein(protein, normalized)
            ko_entries = read_ko_list(ref["ko_list"])
            profile_kos = {path.stem for path in Path(ref["profiles_dir"]).glob("*.hmm")}
            config = {"profile": ref["profiles_dir"], "ko_list": ref["ko_list"], "cpu": threads}
            for dependency in ["hmmsearch", "parallel"]:
                if software[dependency]:
                    config[dependency] = software[dependency]["path"]
            # JSON is valid YAML; an explicit config prevents ambient config.yml changes.
            write_json(result / "execution_config.json", config)
            argv = [str(executable), "-c", str(result / "execution_config.json"),
                    "-p", ref["profiles_dir"], "-k", ref["ko_list"], "--cpu", str(threads),
                    "--tmp-dir", str(work / "hmmsearch"), "-f", "detail-tsv", "-T", "1",
                    "--no-report-unannotated", "-o", str(result / "detail.tsv"), str(normalized)]
            with open(result / "stdout.log", "w") as stdout, open(result / "stderr.log", "w") as stderr:
                subprocess.run(argv, cwd=work, stdout=stdout, stderr=stderr, check=True)
            hits, genes = parse_detail(result / "detail.tsv", species, proteins, ko_entries, profile_kos)
            write_tsv(result / "gene_kos.tsv", HIT_FIELDS, hits)
            write_tsv(result / "genes.tsv", GENE_FIELDS, genes)
            if (file_record(protein) != identity["protein"]
                    or file_record(reference) != identity["reference"]
                    or file_record(ref["inventory_path"]) != identity["reference_inventory"]):
                raise ValueError("KofamScan inputs changed during annotation")
            record = {"schema_version": 1, "completed_at": now(), "fingerprint": fingerprint,
                      "identity": identity, "reference_release": ref["release"], "command": argv,
                      "normalized_protein_sha256": sha256(normalized),
                      "terminal_stop_stripped_count": sum(p["terminal_stop_stripped"] for p in proteins),
                      "gene_status_counts": dict(Counter(g["assignment_status"] for g in genes)),
                      "quantitative_policy": "unique accepted KO only; ambiguous assignments excluded",
                      "acceptance_policy": "upstream pre-rounding marker, verified against KO threshold",
                      "results": [{**file_record(result / name), "path": str(out / name)}
                                  for name in RESULT_NAMES]}
            write_json(result / "provenance.json", record)
            _publish(result, out)
            write_json(status_path, {"state": "success", "completed_at": now(),
                                     "fingerprint": fingerprint, "reused": False})
        except BaseException as error:
            write_json(status_path, {"state": "failed", "time": now(), "fingerprint": fingerprint,
                                     "work": str(work), "error": str(error)})
            raise
        else:
            shutil.rmtree(work)
            return record


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ["protein", "species", "reference", "output-dir", "work-dir"]:
        parser.add_argument(f"--{flag}", required=True)
    parser.add_argument("--command", default="exec_annotation")
    parser.add_argument("--threads", type=int, default=1)
    run(**vars(parser.parse_args()))
