#!/usr/bin/env python3
"""Freeze local KOfam profiles and KEGG membership tables in a new snapshot."""
import argparse
import csv
import fcntl
import math
import shutil
import tempfile
import time
import urllib.error
import urllib.request
from pathlib import Path
import re

from common import now, sha256, write_json, write_tsv


FORMAT_VERSION = 1
REFERENCE_KIND = "phenoradar-kegg-kofam"
LINK_URLS = {
    "module": "https://rest.kegg.jp/link/module/ko",
    "pathway": "https://rest.kegg.jp/link/pathway/ko",
}
KO_PATTERN = re.compile(r"K[0-9]{5}")


def read_ko_list(path):
    """Validate the positional, 12-column table expected by KofamScan."""
    entries = {}
    with open(path, encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
        if (len(header) != 12 or header[:3] != ["knum", "threshold", "score_type"]
                or header[-1] != "definition"):
            raise ValueError(f"invalid KOfam ko_list header (expected 12 columns): {path}")
        for lineno, row in enumerate(reader, 2):
            if len(row) != 12 or not KO_PATTERN.fullmatch(row[0]):
                raise ValueError(f"invalid ko_list row at {path}:{lineno}")
            ko, threshold, score_type = row[:3]
            if ko in entries:
                raise ValueError(f"duplicate KO in ko_list: {ko}")
            if threshold != "-":
                try:
                    value = float(threshold)
                except ValueError as error:
                    raise ValueError(f"invalid threshold for {ko}") from error
                if not math.isfinite(value) or score_type not in {"full", "domain"}:
                    raise ValueError(f"invalid threshold/score_type for {ko}")
            elif score_type not in {"full", "domain", "-"}:
                raise ValueError(f"invalid score_type for {ko}")
            entries[ko] = {"threshold": threshold, "score_type": score_type}
    if not entries:
        raise ValueError(f"empty KOfam ko_list: {path}")
    return entries


def validate_profile(path):
    """Reject empty, truncated, non-protein, or incorrectly named HMM files."""
    ko = path.stem
    if not KO_PATTERN.fullmatch(ko):
        raise ValueError(f"profile filename must be KNNNNN.hmm: {path}")
    name, alphabet, last, models = None, None, "", 0
    with open(path, encoding="utf-8") as handle:
        first = next(handle, "")
        if not first.startswith("HMMER3/"):
            raise ValueError(f"expected an HMMER3 profile: {path}")
        for line in handle:
            if line.strip():
                last = line.strip()
            fields = line.split()
            if fields and fields[0] == "NAME":
                models += 1
                name = fields[1] if len(fields) == 2 else None
            elif fields and fields[0] == "ALPH":
                alphabet = fields[1] if len(fields) == 2 else None
    if name != ko or alphabet != "amino" or last != "//" or models != 1:
        raise ValueError(f"incomplete or mismatched protein HMM profile: {path}")
    return ko


def normalize_links(path, kind):
    """Read raw KEGG links in either direction; retain all memberships."""
    if kind not in LINK_URLS:
        raise ValueError(f"unsupported KEGG membership kind: {kind}")
    pairs = set()
    with open(path, encoding="utf-8", newline="") as handle:
        for lineno, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.strip().split("\t")
            if len(fields) != 2:
                raise ValueError(f"invalid KEGG {kind} links at {path}:{lineno}")
            a, b = fields
            if re.fullmatch(r"(?:ko:)?K[0-9]{5}", b):
                a, b = b, a
            if not re.fullmatch(r"(?:ko:)?K[0-9]{5}", a):
                raise ValueError(f"invalid KO in {path}:{lineno}")
            ko = a.removeprefix("ko:")
            if kind == "module":
                match = re.fullmatch(r"(?:(?:md|module):)?(M[0-9]{5})", b)
                member = match.group(1) if match else None
            else:
                match = re.fullmatch(r"(?:(?:path|pathway):)?(?:map|ko)([0-9]{5})", b)
                member = "map" + match.group(1) if match else None
            if member is None:
                raise ValueError(f"invalid reference {kind} ID in {path}:{lineno}")
            pairs.add((ko, member))
    if not pairs:
        raise ValueError(f"empty KEGG {kind} links: {path}")
    return [{"ko": ko, kind: member} for ko, member in sorted(pairs)]


def download_links(url, destination, attempts=4):
    """Fetch only small KEGG REST mappings, at most two requests per second."""
    if attempts < 1:
        raise ValueError("download attempts must be positive")
    for attempt in range(attempts):
        # A delay before every request also covers consecutive mapping downloads.
        time.sleep(0.5)
        request = urllib.request.Request(url, headers={"User-Agent": "phenoradar_prep/KEGG-reference"})
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                if response.status != 200:
                    raise ValueError(f"unexpected HTTP status {response.status} for {url}")
                payload = response.read(64 * 1024 * 1024 + 1)
                if len(payload) > 64 * 1024 * 1024:
                    raise ValueError(f"KEGG mapping exceeds 64 MiB: {url}")
                destination.write_bytes(payload)
                return {"url": url, "retrieved_at": now(),
                        "last_modified": response.headers.get("Last-Modified"),
                        "etag": response.headers.get("ETag")}
        except (urllib.error.URLError, TimeoutError, ConnectionError) as error:
            if isinstance(error, urllib.error.HTTPError) and error.code not in {429, 500, 502, 503, 504}:
                raise
            if attempt + 1 == attempts:
                raise
            delay = 2 ** attempt
            if isinstance(error, urllib.error.HTTPError):
                retry_after = error.headers.get("Retry-After", "") if error.headers else ""
                if retry_after.isdigit():
                    delay = max(delay, min(int(retry_after), 60))
            time.sleep(delay)


def _inventory_entry(root, path):
    return {"relative_path": path.relative_to(root).as_posix(),
            "bytes": path.stat().st_size, "sha256": sha256(path)}


def _copied_source_record(source, frozen):
    # Hash the bytes actually frozen, not a source another process could modify.
    return {"path": str(source), "bytes": frozen.stat().st_size, "sha256": sha256(frozen)}


def prepare(profiles_dir, ko_list, reference_dir, module_links=None, pathway_links=None,
            release="unspecified"):
    """Build a relocatable snapshot and publish only after successful validation."""
    root = Path(reference_dir).absolute()
    if root.exists() or root.is_symlink():
        raise FileExistsError(f"reference path already exists; use a new directory: {root}")
    if not isinstance(release, str) or not release.strip():
        raise ValueError("release label must not be empty")
    source_dir = Path(profiles_dir).resolve(strict=True)
    source_ko_list = Path(ko_list).resolve(strict=True)
    if not source_dir.is_dir() or not source_ko_list.is_file():
        raise ValueError("profiles-dir must be a directory and ko-list must be a file")
    ko_entries = read_ko_list(source_ko_list)
    profiles = sorted(source_dir.rglob("*.hmm"))
    if not profiles:
        raise ValueError(f"no extracted *.hmm profiles found: {source_dir}")
    profile_ids = set()
    for path in profiles:
        ko = validate_profile(path)
        if ko not in ko_entries:
            raise ValueError(f"profile {ko} is absent from ko_list")
        if ko in profile_ids:
            raise ValueError(f"duplicate profile for {ko}")
        profile_ids.add(ko)

    root.parent.mkdir(parents=True, exist_ok=True)
    root = root.parent.resolve() / root.name
    # Keep the lock file: unlinking an acquired lock permits competing lock inodes.
    with open(root.parent / f".{root.name}.prepare.lock", "a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if root.exists() or root.is_symlink():
            raise FileExistsError(f"reference path already exists; use a new directory: {root}")
        staging = Path(tempfile.mkdtemp(prefix=f".{root.name}.preparing-", dir=root.parent))
        try:
            (staging / "profiles").mkdir()
            (staging / "raw").mkdir()
            for source in profiles:
                shutil.copyfile(source, staging / "profiles" / source.name)
            shutil.copyfile(source_ko_list, staging / "ko_list")
            # Validate the copied input, including consistency if a source changed during setup.
            frozen_entries = read_ko_list(staging / "ko_list")
            for path in sorted((staging / "profiles").glob("*.hmm")):
                if validate_profile(path) not in frozen_entries:
                    raise ValueError(f"copied profile absent from ko_list: {path.name}")

            sources = {"profiles_dir": str(source_dir),
                       "ko_list": _copied_source_record(source_ko_list, staging / "ko_list")}
            counts = {"profiles": len(profiles), "ko_list_entries": len(frozen_entries)}
            for kind, local in [("module", module_links), ("pathway", pathway_links)]:
                raw = staging / "raw" / f"ko_{kind}_links.tsv"
                if local is None:
                    sources[f"{kind}_links"] = download_links(LINK_URLS[kind], raw)
                else:
                    source = Path(local).resolve(strict=True)
                    shutil.copyfile(source, raw)
                    sources[f"{kind}_links"] = {
                        "local_file": _copied_source_record(source, raw), "copied_at": now()}
                rows = normalize_links(raw, kind)
                write_tsv(staging / f"ko_{kind}s.tsv", ["ko", kind], rows)
                counts[f"{kind}_memberships"] = len(rows)

            paths = sorted(p for p in staging.rglob("*") if p.is_file())
            write_json(staging / "files.json", [_inventory_entry(staging, path) for path in paths])
            write_json(staging / "reference.json", {
                "format_version": FORMAT_VERSION, "kind": REFERENCE_KIND,
                "created_at": now(), "release": release.strip(),
                "profiles_dir": "profiles", "ko_list": "ko_list",
                "ko_modules": "ko_modules.tsv", "ko_pathways": "ko_pathways.tsv",
                "raw_module_links": "raw/ko_module_links.tsv",
                "raw_pathway_links": "raw/ko_pathway_links.tsv",
                "inventory": _inventory_entry(staging, staging / "files.json"),
                "sources": sources, "counts": counts,
                "module_semantics": "KO membership only; no reaction completeness or flux inference",
                "pathway_identifiers": "mapNNNNN (koNNNNN aliases deduplicated)",
            })
            from verify_kegg_reference import verify
            verify(staging / "reference.json", full=True)
            if root.exists() or root.is_symlink():
                raise FileExistsError(f"reference path appeared during setup: {root}")
            staging.rename(root)
        finally:
            if staging.exists():
                shutil.rmtree(staging)
    return root / "reference.json"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profiles-dir", required=True, help="Extracted local KOfam *.hmm profiles")
    parser.add_argument("--ko-list", required=True, help="Local decompressed KOfam ko_list")
    parser.add_argument("--reference-dir", "--output", dest="reference_dir", required=True)
    parser.add_argument("--module-links", help="Local raw KEGG KO/MODULE links; otherwise fetch via REST")
    parser.add_argument("--pathway-links", help="Local raw KEGG KO/PATHWAY links; otherwise fetch via REST")
    parser.add_argument("--release", default="unspecified", help="User-supplied KOfam release/source label")
    print(prepare(**vars(parser.parse_args())))
