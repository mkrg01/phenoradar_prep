"""Internal stage receipts and verification for reusable pipeline products."""
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
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextvars import ContextVar
from pathlib import Path

from common import file_record, read_tsv, sha256, write_json
from translate_cds import fasta_ids
from sample_identity import annotate
from verification_cache import active_cache, compact_rows, expand_rows

COUNTS = ["busco_cds_single", "busco_cds_duplicated", "busco_cds_fragmented", "busco_cds_missing", "busco_cds_total"]
_CACHE_ONLY = ContextVar('cache_only_verification', default=False)


class ProductIntegrityError(ValueError):
    """A registered artifact is missing or no longer matches its recorded content."""


class UncachedVerification(ValueError):
    """Preparation needs a file read and must move to an allocated Slurm job."""


@contextmanager
def cached_verification():
    token = _CACHE_ONLY.set(True)
    try:
        yield
    finally:
        _CACHE_ONLY.reset(token)


def verification_workers(requested=4):
    allocated = os.environ.get('SLURM_CPUS_PER_TASK')
    return max(1, min(requested, int(allocated))) if allocated else max(1, requested)


SAFE = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*\Z")
TAXONOMY_ID_FIELDS = ("taxid", "taxid_domain", "taxid_kingdom", "taxid_phylum", "taxid_class",
                      "taxid_order", "taxid_family", "taxid_genus", "taxid_species")


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def normalize_taxonomy_ids(row):
    """Canonicalize AMALGKIT taxids without inferring taxonomy; taxid is required."""
    row = dict(row)
    for field in TAXONOMY_ID_FIELDS:
        if field not in row:
            continue
        original = row[field]
        value = original.strip()
        if not value and field != "taxid":
            row[field] = ""
            continue
        match = re.fullmatch(r"([0-9]+)(?:\.0+)?", value)
        integer = int(match.group(1)) if match else 0
        if integer < 1:
            expected = "positive taxid" if field == "taxid" else f"positive integer {field}"
            raise ValueError(f"{expected} required: run={row['run']!r}, value={original!r}")
        row[field] = str(integer)
    return row


def identities(metadata, source=None):
    path = Path(metadata).resolve()
    before = stat_identity(path)
    if source is not None and source.get('stat') != before:
        raise ValueError(f'file changed while reading sample metadata: {path}')
    cache = active_cache()
    checksum = source['sha256'] if source is not None else None
    indexed = cache.document(path, 'sample-metadata-v1', before, checksum) if cache is not None else None
    if indexed is not None:
        try:
            items = indexed['items']
            rows = expand_rows(indexed['columns'], [i['row'] for i in items])
            result = [dict(item, row=row) for item, row in zip(items, rows)]
            if before != stat_identity(path):
                raise ValueError(f'file changed while reading sample metadata: {path}')
            return indexed['fields'], result
        except (KeyError, TypeError, IndexError):
            pass
    fields, result = _read_identities(path)
    if cache is not None:
        checksum = checksum or sha256(path)
        columns, rows = compact_rows(i['row'] for i in result)
        indexed = {'fields': fields, 'columns': columns,
                   'items': [dict(item, row=row) for item, row in zip(result, rows)]}
        if before != stat_identity(path):
            raise ValueError(f'file changed while reading sample metadata: {path}')
        cache.remember_document(path, 'sample-metadata-v1', before, checksum, indexed)
    return fields, result


def _read_identities(metadata):
    with open(metadata, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames
        if not fields or len(set(fields)) != len(fields) or not {"scientific_name", "run", "taxid"} <= set(fields):
            raise ValueError("metadata requires unique columns including scientific_name, run, taxid")
        rows = list(reader)
    if not rows:
        raise ValueError("metadata must contain at least one species")
    species, runs, odb_names, taxids = set(), set(), set(), {}
    result = []
    for row in rows:
        if None in row or any(v is None for v in row.values()):
            raise ValueError("malformed metadata row")
        name = row["scientific_name"]
        row = annotate(normalize_taxonomy_ids(row))
        label = row["analysis_sample_id"]
        odb = label.replace("-", "_")
        if not SAFE.fullmatch(label) or not SAFE.fullmatch(row["run"]):
            raise ValueError(f"unsafe species/run: {name!r}, {row['run']!r}")
        if label in species or row["run"] in runs or odb in odb_names:
            raise ValueError("dataset metadata requires unique runs and normalized sample identities")
        taxon = row["species_id"]
        if taxon in taxids and taxids[taxon] != row["taxid"]:
            raise ValueError(f"conflicting taxids for the same species: {name}")
        taxids[taxon] = row["taxid"]
        if row.get("reference_id") and not re.fullmatch(r"[0-9a-f]{64}", row["reference_id"]):
            raise ValueError("reference_id must be a registered CDS SHA256")
        species.add(label); runs.add(row["run"]); odb_names.add(odb)
        result.append({"species": label, "odb_species": odb, "row": row})
    fields = list(fields) + [k for k in ("species_id", "analysis_sample_id") if k not in fields]
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
    path = Path(path).resolve()
    before = stat_identity(path)
    cache = active_cache()
    checksum = cache.checksum(before) if cache is not None else None
    entry = {'path': str(path), 'bytes': before[2], 'sha256': checksum} if checksum else file_record(path)
    if before != stat_identity(path):
        raise ValueError(f"file changed while registering: {path}")
    if cache is not None:
        cache.remember(path, entry['sha256'], before)
    return dict(entry, stat=before)


def verify(entry):
    path = Path(entry["path"])
    if not path.is_file():
        raise ProductIntegrityError(f"registered file missing: {path}")
    before = stat_identity(path)
    if entry.get("stat") == before:
        cache = active_cache()
        if cache is not None:
            cache.remember(path, entry['sha256'], before)
        return path
    cache = active_cache()
    if cache is not None and cache.matches(path, entry["sha256"], before):
        return path
    if _CACHE_ONLY.get():
        raise UncachedVerification(f'uncached file verification: {path}')
    if sha256(path) != entry["sha256"]:
        raise ProductIntegrityError(f"registered file changed: {path}")
    if cache is not None:
        try:
            if before == stat_identity(path):
                cache.remember(path, entry["sha256"], before)
        except OSError:
            pass  # Preserve verify's result, but never cache a file that changed after hashing.
    return path


def verification_cached(entry):
    """Check receipts/cache using stat alone, without reading file contents."""
    path = Path(entry['path'])
    if not path.is_file():
        return False
    try:
        identity = stat_identity(path)
    except OSError:
        return False
    cache = active_cache()
    return (entry.get('stat') == identity or
            (cache is not None and cache.matches(path, entry['sha256'], identity)))


def verify_many(entries, workers=4, progress=None):
    """Hash cache misses concurrently; all SQLite access stays on the caller thread.

    Return an error per entry, so a damaged sample does not discard healthy ones.
    Aliases with the same checksum and full file identity share a single read.
    """
    entries = list(entries)
    errors, pending = [None] * len(entries), {}
    cache = active_cache()
    checked = 0
    for index, entry in enumerate(entries):
        path = Path(entry['path'])
        try:
            if not path.is_file():
                raise ProductIntegrityError(f'registered file missing: {path}')
            identity = stat_identity(path)
            if (entry.get('stat') == identity or
                    (cache is not None and cache.matches(path, entry['sha256'], identity))):
                if cache is not None:
                    cache.remember(path, entry['sha256'], identity)
                checked += 1
                continue
            key = (entry['sha256'], tuple(identity))
            pending.setdefault(key, []).append((index, path))
        except OSError as error:
            errors[index] = ProductIntegrityError(f'registered file unavailable: {path}: {error}')
            checked += 1
        except ProductIntegrityError as error:
            errors[index] = error
            checked += 1
    if progress:
        progress(checked, len(entries))
    if pending:
        if _CACHE_ONLY.get():
            raise UncachedVerification('uncached product verification')
        with ThreadPoolExecutor(max_workers=verification_workers(workers)) as executor:
            futures = {executor.submit(sha256, aliases[0][1]): (key, aliases)
                       for key, aliases in pending.items()}
            for future in as_completed(futures):
                (expected, identity), aliases = futures[future]
                try:
                    actual = future.result()
                    error = None
                except OSError as failure:
                    actual, error = None, str(failure)
                for index, path in aliases:
                    if actual != expected:
                        errors[index] = ProductIntegrityError(
                            f'registered file unavailable: {path}: {error}' if error else
                            f'registered file changed: {path}')
                    elif cache is not None:
                        try:
                            if tuple(stat_identity(path)) == identity:
                                cache.remember(path, expected, identity)
                        except OSError:
                            pass
                checked += len(aliases)
                if progress:
                    progress(checked, len(entries))
    return errors


def product_file(entry, verify_files=True):
    """Allow previews to check presence without claiming content verification."""
    if verify_files:
        return verify(entry)
    path = Path(entry['path'])
    if not path.is_file():
        raise ProductIntegrityError(f'registered file missing: {path}')
    return path


def product_json(entry, verify_files=True):
    path = product_file(entry, verify_files)
    try:
        return json.loads(path.read_text())
    except (UnicodeError, json.JSONDecodeError) as error:
        raise ProductIntegrityError(f'registered JSON file damaged: {path}') from error


@contextmanager
def locked(path):
    """Hold an exclusive lock; its handle can be inherited by a child process."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        yield handle


def link_file(source, target):
    source, target = Path(source), Path(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        if os.path.samefile(source, target):
            return
        if sha256(source) != sha256(target):
            raise ValueError(f"refusing to replace different staged input: {target}")
        return
    try:
        os.link(source, target)
    except OSError:
        shutil.copy2(source, target)


def stage_file(entry, target, final=None):
    """Import verified bytes and record the identity after our hard-link operation."""
    source = Path(entry['path'])
    before = stat_identity(source)
    verify(entry)
    if before != stat_identity(source):
        raise ProductIntegrityError(f'registered file changed while staging: {source}')
    target = Path(target)
    link_file(source, target)
    after, identity = stat_identity(source), stat_identity(target)
    if before[:4] != after[:4]:
        raise ProductIntegrityError(f'registered file changed while staging: {source}')
    if identity[:2] == after[:2]:
        # os.link changes ctime without writing bytes. Register that controlled
        # change, then keep requiring the complete identity on subsequent use.
        if identity != after:
            raise ProductIntegrityError(f'registered file changed while staging: {source}')
    else:
        # Cross-filesystem copies must be checked; only a hard link shares proof.
        if _CACHE_ONLY.get():
            raise UncachedVerification(f'copied file needs verification: {target}')
        if before != after or sha256(target) != entry['sha256']:
            raise ProductIntegrityError(f'registered file changed while staging: {source}')
        if identity != stat_identity(target):
            raise ProductIntegrityError(f'copied file changed while staging: {target}')
    cache = active_cache()
    if cache is not None:
        cache.remember(source, entry['sha256'], after)
        cache.remember(final or target, entry['sha256'], identity)
    return dict(entry, path=str(final or target), stat=identity)


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
            if ref.get("schema_version", 1) >= 2 and any(ref[k] != item["row"][k] for k in ("species_id", "run")):
                raise ValueError("registered sample identity differs from metadata")
        else:
            ref = {"schema_version": 2, "species": item["species"], "taxid": item["row"]["taxid"],
                   "species_id": item["row"]["species_id"], "run": item["row"]["run"],
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


def resolve(store, item, lineage, need_full=False, verify_files=True):
    root = Path(store) / item["species"]
    requested = item["row"].get("reference_id", "")
    candidates = [root / requested / "reference.json"] if requested else sorted(root.glob("*/reference.json"))
    if requested and not candidates[0].exists():
        raise ValueError(f"unknown reference_id for {item['species']}: {requested}")
    if len(candidates) > 1:
        raise ValueError(f"multiple references for {item['species']}; specify reference_id in metadata")
    if not candidates:
        return {"reference": None, "busco": None, "quant": None}
    try:
        ref = json.loads(candidates[0].read_text())
    except (UnicodeError, json.JSONDecodeError) as error:
        raise ProductIntegrityError(f'registered JSON file damaged: {candidates[0]}') from error
    if ref["taxid"] != item["row"]["taxid"] or ref["species"] != item["species"]:
        raise ValueError(f"registered identity differs: {item['species']}")
    if ref.get("schema_version", 1) >= 2 and any(ref[k] != item["row"][k] for k in ("species_id", "run")):
        raise ValueError("registered sample identity differs from metadata")
    product_file(ref["cds"], verify_files)
    bus = candidates[0].parent / "busco.json"
    try:
        bus = json.loads(bus.read_text()) if bus.exists() else None
    except (UnicodeError, json.JSONDecodeError) as error:
        raise ProductIntegrityError(f'registered JSON file damaged: {bus}') from error
    if bus:
        if bus["reference_id"] != ref["reference_id"] or bus["lineage"] != lineage:
            raise ValueError(f"registered BUSCO reference/lineage differs: {item['species']}")
        counts(bus["counts"])
        for key in ("full", "short"):
            if bus.get(key): product_file(bus[key], verify_files)
    quant = candidates[0].parent / "quant" / (item["row"]["run"] + ".json")
    try:
        quant = json.loads(quant.read_text()) if quant.exists() else None
    except (UnicodeError, json.JSONDecodeError) as error:
        raise ProductIntegrityError(f'registered JSON file damaged: {quant}') from error
    if quant:
        if quant["reference_id"] != ref["reference_id"] or quant["run"] != item["row"]["run"]:
            raise ValueError("quantification receipt identity differs")
        product_file(quant["abundance"], verify_files)
        for key, value in quant.get("sample", {}).items():
            supplied = item["row"].get(key, "")
            if supplied and value and supplied != value:
                raise ValueError(f"run metadata changed for cached quantification: {item['row']['run']}: {key}")
    return {"reference": ref, "busco": bus if bus and (not need_full or bus.get("full")) else None,
            "quant": quant, "assessment": bus}
