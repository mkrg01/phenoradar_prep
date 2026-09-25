"""Species/reference-bound receipts for imported and GeneGalleon RNA-seq products."""
import csv
import fcntl
import gzip
import hashlib
import json
import math
import os
import re
import shutil
from contextlib import contextmanager
from pathlib import Path

from common import file_record, read_tsv, sha256, write_json
from translate_cds import fasta_ids

COUNTS = ["busco_cds_single", "busco_cds_duplicated", "busco_cds_fragmented", "busco_cds_missing", "busco_cds_total"]
SAFE = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*\Z")


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def identities(metadata):
    with open(metadata, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames
        if not fields or len(set(fields)) != len(fields) or not {"scientific_name", "run", "taxid"} <= set(fields):
            raise ValueError("metadata requires unique columns including scientific_name, run, taxid")
        rows = list(reader)
    if not rows:
        raise ValueError("metadata must contain at least one species")
    species, runs, odb_names = set(), set(), set()
    result = []
    for row in rows:
        if None in row or any(v is None for v in row.values()):
            raise ValueError("malformed metadata row")
        name = row["scientific_name"]
        label = name.replace(" ", "_")
        odb = label.replace("-", "_")
        if not SAFE.fullmatch(label) or not SAFE.fullmatch(row["run"]):
            raise ValueError(f"unsafe species/run: {name!r}, {row['run']!r}")
        if label in species or row["run"] in runs or odb in odb_names:
            raise ValueError("dataset metadata requires one row/run per species and unique normalized identities")
        if not row["taxid"].isdigit() or int(row["taxid"]) < 1:
            raise ValueError(f"positive taxid required: {name}")
        if row.get("reference_id") and not re.fullmatch(r"[0-9a-f]{64}", row["reference_id"]):
            raise ValueError("reference_id must be a registered CDS SHA256")
        species.add(label); runs.add(row["run"]); odb_names.add(odb)
        result.append({"species": label, "odb_species": odb, "row": row})
    return fields, sorted(result, key=lambda r: r["species"])


def normalize_private_paths(items, metadata):
    base = Path(metadata).resolve().parent
    for item in items:
        row = item["row"]
        if row.get("private_file", "").lower() == "yes":
            for key in ("read1_path", "read2_path"):
                if row.get(key):
                    path = Path(row[key])
                    row[key] = str((base / path).resolve())


def stat_identity(path):
    s = Path(path).stat()
    return [s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns]


def record(path):
    before = stat_identity(path)
    entry = file_record(path)
    if before != stat_identity(path):
        raise ValueError(f"file changed while registering: {path}")
    return dict(entry, stat=before)


def verify(entry):
    path = Path(entry["path"])
    if not path.is_file():
        raise ValueError(f"registered file missing: {path}")
    if entry.get("stat") != stat_identity(path) and sha256(path) != entry["sha256"]:
        raise ValueError(f"registered file changed: {path}")
    return path


@contextmanager
def locked(path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        yield


def link_file(source, target):
    source, target = Path(source), Path(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        if sha256(source) != sha256(target):
            raise ValueError(f"refusing to replace different staged input: {target}")
        return
    try:
        os.link(source, target)
    except OSError:
        shutil.copy2(source, target)


def counts(values):
    result = {}
    for key in COUNTS:
        value = str(values[key])
        if not re.fullmatch(r"[0-9]+", value):
            raise ValueError(f"invalid BUSCO count: {key}={value}")
        result[key] = int(value)
    if result[COUNTS[-1]] <= 0 or sum(result[c] for c in COUNTS[:-1]) != result[COUNTS[-1]]:
        raise ValueError("BUSCO counts must sum to a positive total")
    return result


def busco_full(path, lineage, cds):
    """Read full tables, including duplicated marker rows, without rerunning BUSCO."""
    ids = set(fasta_ids(cds, compressed=str(cds).endswith(".gz")))
    markers, observed_lineage = {}, None
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("# The lineage dataset is:"):
                observed_lineage = line.split(":", 1)[1].strip().split()[0]
            if not line.strip() or line.startswith("#"):
                continue
            row = line.rstrip("\n").split("\t")
            if len(row) < 2 or row[1] not in {"Complete", "Duplicated", "Fragmented", "Missing"}:
                raise ValueError(f"invalid BUSCO full table: {path}")
            marker, status = row[:2]
            if marker in markers and (markers[marker] != status or status != "Duplicated"):
                raise ValueError(f"conflicting BUSCO marker: {marker}")
            markers[marker] = status
            if status != "Missing":
                if len(row) < 3 or row[2].rsplit(":", 1)[0] not in ids:
                    raise ValueError(f"BUSCO sequence does not belong to CDS: {row}")
    if observed_lineage != lineage:
        raise ValueError(f"BUSCO lineage mismatch: expected {lineage}, found {observed_lineage}: {path}")
    result = {key: sum(v == status for v in markers.values()) for key, status in zip(COUNTS[:-1], ["Complete", "Duplicated", "Fragmented", "Missing"])}
    result[COUNTS[-1]] = len(markers)
    return counts(result)


def validate_abundance(path, cds):
    ids = set(fasta_ids(cds, compressed=str(cds).endswith(".gz")))
    observed, total = set(), 0.0
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"target_id", "tpm"} <= set(reader.fieldnames or []):
            raise ValueError(f"abundance requires target_id and tpm: {path}")
        for row in reader:
            target, value = row["target_id"], float(row["tpm"])
            if target in observed or target not in ids or not math.isfinite(value) or value < 0:
                raise ValueError(f"invalid abundance target/value: {path}: {target}")
            observed.add(target); total += value
    if observed != ids or total <= 0:
        raise ValueError(f"abundance target set differs from CDS or has no positive TPM: {path}")


def reference_path(store, species, ref):
    return Path(store) / species / ref / "reference.json"


def register_reference(store, item, cds, provenance=None):
    ids = fasta_ids(cds, compressed=str(cds).endswith(".gz"))
    if any(not x.startswith(item["species"] + "_") for x in ids):
        raise ValueError(f"CDS IDs must have species prefix: {item['species']}")
    entry = record(cds)
    path = reference_path(store, item["species"], entry["sha256"])
    with locked(path.parent / ".lock"):
        if path.exists():
            ref = json.loads(path.read_text())
            verify(ref["cds"])
            if ref["taxid"] != item["row"]["taxid"]:
                raise ValueError("registered species taxid differs from metadata")
        else:
            ref = {"schema_version": 1, "species": item["species"], "taxid": item["row"]["taxid"],
                   "reference_id": entry["sha256"], "cds": entry, "provenance": provenance or {"source": "legacy"}}
            write_json(path, ref)
    return ref


def register_busco(store, ref, summary=None, full=None, short=None, lineage="embryophyta_odb12", provenance=None):
    cds = verify(ref["cds"])
    values = busco_full(full, lineage, cds) if full else counts(summary)
    if summary and counts(summary) != values:
        raise ValueError(f"BUSCO summary/full disagree: {ref['species']}")
    record_path = reference_path(store, ref["species"], ref["reference_id"]).parent / "busco.json"
    receipt = {"schema_version": 1, "reference_id": ref["reference_id"], "lineage": lineage,
               "counts": values, "full": record(full) if full else None,
               "short": record(short) if short else None, "provenance": provenance or {"source": "legacy"}}
    with locked(record_path.parent / ".lock"):
        if record_path.exists():
            old = json.loads(record_path.read_text())
            if old["counts"] != values or old["lineage"] != lineage:
                raise ValueError("BUSCO result conflicts with registered reference")
            for key in ("full", "short"):
                if old.get(key):
                    verify(old[key])
                    if receipt[key] and old[key]["sha256"] != receipt[key]["sha256"]:
                        raise ValueError("BUSCO file differs from registered result")
                    receipt[key] = old[key]
        write_json(record_path, receipt)
    return receipt


def register_quant(store, ref, item, abundance, provenance=None):
    validate_abundance(abundance, verify(ref["cds"]))
    path = reference_path(store, ref["species"], ref["reference_id"]).parent / "quant" / (item["row"]["run"] + ".json")
    receipt = {"schema_version": 1, "reference_id": ref["reference_id"], "run": item["row"]["run"],
               "abundance": record(abundance), "sample": {k: item["row"].get(k, "") for k in
               ("private_file", "lib_layout", "read1_path", "read2_path")},
               "provenance": provenance or {"source": "legacy"}}
    with locked(path.parent / ".lock"):
        if path.exists():
            old = json.loads(path.read_text())
            verify(old["abundance"])
            if old["abundance"]["sha256"] != receipt["abundance"]["sha256"]:
                raise ValueError("quantification differs from registered run/reference")
            return old
        write_json(path, receipt)
    return receipt


def resolve(store, item, lineage, need_full=False):
    root = Path(store) / item["species"]
    requested = item["row"].get("reference_id", "")
    candidates = [root / requested / "reference.json"] if requested else sorted(root.glob("*/reference.json"))
    if requested and not candidates[0].exists():
        raise ValueError(f"unknown reference_id for {item['species']}: {requested}")
    if len(candidates) > 1:
        raise ValueError(f"multiple references for {item['species']}; specify reference_id in metadata")
    if not candidates:
        return {"reference": None, "busco": None, "quant": None}
    ref = json.loads(candidates[0].read_text())
    if ref["taxid"] != item["row"]["taxid"] or ref["species"] != item["species"]:
        raise ValueError(f"registered identity differs: {item['species']}")
    verify(ref["cds"])
    bus = candidates[0].parent / "busco.json"
    bus = json.loads(bus.read_text()) if bus.exists() else None
    if bus:
        if bus["reference_id"] != ref["reference_id"] or bus["lineage"] != lineage:
            raise ValueError(f"registered BUSCO reference/lineage differs: {item['species']}")
        counts(bus["counts"])
        for key in ("full", "short"):
            if bus.get(key): verify(bus[key])
    quant = candidates[0].parent / "quant" / (item["row"]["run"] + ".json")
    quant = json.loads(quant.read_text()) if quant.exists() else None
    if quant:
        if quant["reference_id"] != ref["reference_id"] or quant["run"] != item["row"]["run"]:
            raise ValueError("quantification receipt identity differs")
        verify(quant["abundance"])
        for key, value in quant.get("sample", {}).items():
            supplied = item["row"].get(key, "")
            if supplied and value and supplied != value:
                raise ValueError(f"run metadata changed for cached quantification: {item['row']['run']}: {key}")
        # Retired FASTQs need not remain present. If retained, detect replacing
        # their bytes under the same run/path before reusing native quantification.
        for entry in quant.get("provenance", {}).get("raw_inputs", {}).values():
            if Path(entry["path"]).exists(): verify(entry)
    return {"reference": ref, "busco": bus if bus and (not need_full or bus.get("full")) else None,
            "quant": quant, "assessment": bus}


def import_existing(store, input_dir, metadata, lineage="embryophyta_odb12", excluded_runs=()):
    _, items = identities(metadata)
    normalize_private_paths(items, metadata)
    root = Path(input_dir)
    summaries = {}
    summary_path = root / "busco/summary.tsv"
    if summary_path.exists():
        rows = read_tsv(summary_path)
        summaries = {r["Species"]: r for r in rows}
        if len(rows) != len(summaries):
            raise ValueError("duplicate BUSCO species in legacy summary")
    result = []
    for item in items:
        species, run = item["species"], item["row"]["run"]
        if run in excluded_runs:
            result.append({"species": species, "run": run, "status": "excluded"}); continue
        cds = root / "cds" / f"{species}_longestCDS.fa.gz"
        if not cds.is_file():
            result.append({"species": species, "status": "no_cds"}); continue
        ref = register_reference(store, item, cds)
        fulls = [p for p in (root / "busco/full").glob(f"{species}*full*") if p.is_file()]
        # Match exact supported basenames rather than a species-prefix wildcard.
        fulls = [p for p in fulls if p.name in {f"{species}.busco.full.tsv", f"{species}_busco.full.tsv", f"{species}.busco.full.tsv.gz", f"{species}_busco.full.tsv.gz"}]
        if len(fulls) > 1:
            raise ValueError(f"multiple legacy BUSCO full tables: {species}")
        summary = summaries.get(item["row"]["scientific_name"])
        if summary or fulls:
            register_busco(store, ref, summary, fulls[0] if fulls else None, lineage=lineage)
        abundance = root / "quant" / species / run / f"{run}_abundance.tsv"
        if abundance.is_file():
            register_quant(store, ref, item, abundance)
        protein = root / "proteins" / f"{item['odb_species']}_protein.fa"
        if protein.is_file() and protein.with_suffix('.json').is_file():
            from protein_cache import register_translation
            provenance = json.loads(protein.with_suffix('.json').read_text())
            register_translation(cds, protein, protein.with_suffix('.json'),
                                 Path(store) / '.proteins', provenance['translation_table'])
        result.append({"species": species, "status": "registered", "reference_id": ref["reference_id"]})
    return result
