#!/usr/bin/env python3
"""Verify a frozen KEGG/KOfam snapshot, without changing or downloading data."""
import argparse
import json
from pathlib import Path
import re

from common import now, sha256, write_json


def _contained_path(root, relative):
    if not isinstance(relative, str) or not relative or Path(relative).is_absolute():
        raise ValueError(f"invalid relative reference path: {relative!r}")
    path = root / relative
    if ".." in Path(relative).parts or path.resolve() != path.absolute():
        raise ValueError(f"reference path traverses outside snapshot or uses a symlink: {relative}")
    if not path.is_relative_to(root):
        raise ValueError(f"reference path is outside snapshot: {relative}")
    return path


def _check_entry(root, entry, full):
    path = _contained_path(root, entry.get("relative_path"))
    size, digest = entry.get("bytes"), entry.get("sha256")
    if (not isinstance(size, int) or isinstance(size, bool) or size < 0
            or not isinstance(digest, str) or not re.fullmatch(r"[0-9a-f]{64}", digest)):
        raise ValueError(f"invalid reference inventory record: {path}")
    if not path.is_file() or path.stat().st_size != size:
        raise ValueError(f"reference file missing or size mismatch: {path}")
    if full and sha256(path) != digest:
        raise ValueError(f"reference checksum mismatch: {path}")
    return path


def verify(reference, full=True):
    """Return metadata with resolved paths; full=False checks inventory and sizes only."""
    reference = Path(reference).resolve(strict=True)
    root = reference.parent
    record = json.loads(reference.read_text(encoding="utf-8"))
    if record.get("format_version") != 1 or record.get("kind") != "phenoradar-kegg-kofam":
        raise ValueError("unsupported KEGG reference format")
    if not isinstance(record.get("release"), str) or not record["release"].strip():
        raise ValueError("missing KEGG reference release label")
    inventory = record.get("inventory")
    if not isinstance(inventory, dict):
        raise ValueError("missing reference inventory")
    inventory_path = _check_entry(root, inventory, full=True)
    entries = json.loads(inventory_path.read_text(encoding="utf-8"))
    if not isinstance(entries, list) or not entries:
        raise ValueError("empty or malformed reference inventory")
    listed = set()
    for entry in entries:
        if not isinstance(entry, dict):
            raise ValueError("malformed reference inventory record")
        path = _check_entry(root, entry, full=full)
        if path in listed or path in {reference, inventory_path}:
            raise ValueError(f"duplicate or recursive reference inventory entry: {path}")
        listed.add(path)

    result = dict(record)
    keys = ["ko_list", "ko_modules", "ko_pathways", "raw_module_links", "raw_pathway_links"]
    for key in keys:
        path = _contained_path(root, record.get(key))
        if path not in listed:
            raise ValueError(f"required reference artifact is absent from inventory: {key}")
        result[key] = str(path)
    profiles_dir = _contained_path(root, record.get("profiles_dir"))
    profiles = {p for p in listed if p.parent == profiles_dir and p.suffix == ".hmm"}
    if not profiles or not profiles_dir.is_dir():
        raise ValueError("reference inventory contains no KOfam profiles")
    # An unrecorded HMM would silently change KofamScan's search database.
    actual = {p for p in root.rglob("*") if p.is_file() or p.is_symlink()}
    if actual != listed | {reference, inventory_path}:
        raise ValueError("snapshot contains unrecorded or missing reference files")
    if record.get("counts", {}).get("profiles") != len(profiles):
        raise ValueError("profile count does not match reference inventory")
    result.update(profiles_dir=str(profiles_dir), inventory_path=str(inventory_path),
                  reference_json=str(reference), reference_id=sha256(reference),
                  verified_files=len(entries), verification_mode="sha256" if full else "sizes")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True, help="Path to reference.json")
    parser.add_argument("--quick", action="store_true", help="Check inventory hash and file sizes only")
    parser.add_argument("--output", help="Write a verification JSON report outside the snapshot")
    args = parser.parse_args()
    metadata = verify(args.reference, full=not args.quick)
    if args.output:
        if Path(args.output).resolve().is_relative_to(Path(args.reference).resolve().parent):
            raise ValueError("verification output must be outside the frozen snapshot")
        write_json(args.output, {key: metadata[key] for key in [
            "reference_json", "reference_id", "release", "verified_files", "verification_mode"
        ]} | {"verified_at": now()})
    print(f"Verified KEGG reference {metadata['reference_id']} ({'sizes' if args.quick else 'SHA256'})")
