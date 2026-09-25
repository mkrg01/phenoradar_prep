#!/usr/bin/env python3
"""Plan missing ODB species and publish immutable, reusable mapping snapshots."""
import argparse
import fcntl
import json
import os
import shutil
import tempfile
from pathlib import Path

from common import file_record, read_tsv, sha256, write_json
from make_manifests import make


def load_snapshot(root, version, node):
    root = Path(root)
    record = json.loads((root / "snapshot.json").read_text())
    if record.get("schema_version") != 1 or record.get("version") != version or record.get("node") != node:
        raise ValueError(f"existing ODB snapshot schema/version/node differs: {root}")
    proteins = {p["species"]: p for p in record["proteins"]}
    if len(proteins) != len(record["proteins"]):
        raise ValueError(f"existing ODB snapshot has duplicate or missing selected species: {root}")
    if record["annotations"]["path"] != "annotations.tsv":
        raise ValueError(f"snapshot must contain annotations.tsv: {root}")
    return record, proteins


def plan(samples, protein_dir, outdir, cache_dir, existing=None, reference=None,
         version="v12", node=3193, chunk_size=20):
    rows = {r["species"]: r for r in read_tsv(samples)}
    inputs = {s: file_record(Path(protein_dir) / f'{r["odb_species"]}_protein.fa') for s, r in rows.items()}
    roots = [Path(existing)] if existing else []
    roots += sorted((Path(cache_dir) / f"{version}_{node}").glob("*/snapshot.json"))
    roots = [p.parent if p.name == "snapshot.json" else p for p in roots]
    roots = [p for p in roots if not p.name.startswith(".")]
    assigned, sources = {}, []
    reference_hash = sha256(reference) if reference and Path(reference).is_file() else None
    for root in roots:
        record, proteins = load_snapshot(root, version, node)
        members = sorted(set(rows) & set(proteins))
        if not members:
            continue
        if record.get("reference_sha256") and reference_hash and record["reference_sha256"] != reference_hash:
            raise ValueError(f"ODB reference changed since cached mapping: {root}; use a new reference/cache namespace")
        for name in members:
            original = proteins[name]
            if original["odb_species"] != rows[name]["odb_species"] or original["sha256"] != inputs[name]["sha256"]:
                raise ValueError(f"protein differs from existing ODB input: {name}; select a matching snapshot/cache")
        members = [s for s in members if s not in assigned]
        if not members:
            continue
        if sha256(root / "annotations.tsv") != record["annotations"]["sha256"]:
            raise ValueError(f"ODB result changed after completion: {root}")
        index = len(sources)
        assigned.update({s: index for s in members})
        sources.append({"kind": "existing", "root": str(root.resolve()), "species": members,
                        "snapshot_sha256": sha256(root / "snapshot.json")})
    missing = sorted(set(rows) - set(assigned))
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    if missing:
        from common import write_tsv
        path = out / "missing.tsv"
        write_tsv(path, list(next(iter(rows.values()))), [rows[s] for s in missing])
        make(path, protein_dir, out, chunk_size)
        chunks = json.loads((out / "chunks.json").read_text())
    else:
        chunks = []
        write_json(out / "chunks.json", [])
    for chunk in chunks:
        sources.append({"kind": "mapped", **chunk})
    result = {"schema_version": 1, "version": version, "node": node,
              "proteins": inputs, "sources": sources, "chunks": chunks,
              "reused_species": sorted(assigned), "mapped_species": missing}
    write_json(out / "plan.json", result)
    return result


def publish(samples, chunk_dir, protein_dir, cache_dir, label, version="v12", node=3193):
    root = Path(chunk_dir)
    provenance = json.loads((root / "provenance.json").read_text())
    identity = provenance["identity"]
    if identity["version"] != version or identity["node"] != node:
        raise ValueError("completed mapping has a different node/version")
    by_odb = {r["odb_species"]: r["species"] for r in read_tsv(samples)}
    proteins = []
    for entry in identity["proteins"]:
        odb_name = Path(entry["path"]).name.removesuffix("_protein.fa")
        path = Path(protein_dir) / f"{odb_name}_protein.fa"
        if odb_name not in by_odb or sha256(path) != entry["sha256"]:
            raise ValueError("mapping inputs changed before publishing cache")
        proteins.append({"species": by_odb[odb_name], "odb_species": odb_name, **entry})
    for entry in provenance["results"]:
        if sha256(root / Path(entry["path"]).name) != entry["sha256"]:
            raise ValueError("mapping results changed before publishing cache")
    base = Path(cache_dir) / f"{version}_{node}"
    base.mkdir(parents=True, exist_ok=True)
    destination = base / provenance["fingerprint"]
    with open(base / ".publish.lock", "a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if destination.exists():
            saved, _ = load_snapshot(destination, version, node)
            if sha256(destination / "annotations.tsv") != saved["annotations"]["sha256"]:
                raise ValueError(f"corrupt existing cache: {destination}")
            return destination
        staging = Path(tempfile.mkdtemp(prefix=".publishing-", dir=base))
        try:
            shutil.copy2(root / f"{label}.og.annotations", staging / "annotations.tsv")
            shutil.copy2(root / "provenance.json", staging / "mapping_provenance.json")
            record = {"schema_version": 1, "version": version, "node": node,
                      "proteins": proteins, "reference_sha256": identity["reference"]["sha256"],
                      "annotations": dict(file_record(staging / "annotations.tsv"), path="annotations.tsv"),
                      "mapping_provenance": dict(file_record(staging / "mapping_provenance.json"), path="mapping_provenance.json")}
            write_json(staging / "snapshot.json", record)
            os.rename(staging, destination)
        finally:
            if staging.exists():
                shutil.rmtree(staging)
    return destination


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    p = sub.add_parser("plan")
    for name in ("samples", "protein-dir", "outdir", "cache-dir"):
        p.add_argument(f"--{name}", required=True)
    p.add_argument("--existing")
    p.add_argument("--reference")
    p.add_argument("--chunk-size", type=int, default=20)
    q = sub.add_parser("publish")
    for name in ("samples", "chunk-dir", "protein-dir", "cache-dir", "label"):
        q.add_argument(f"--{name}", required=True)
    for command in (p, q):
        command.add_argument("--version", default="v12")
        command.add_argument("--node", type=int, default=3193)
    args = vars(parser.parse_args())
    action = args.pop("action")
    (plan if action == "plan" else publish)(**args)
