#!/usr/bin/env python3
"""Download and freeze a missing KOfam/KEGG reference; reuse completed downloads."""
import argparse
import fcntl
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import tarfile
import tempfile
import time
from urllib.request import Request, urlopen

from common import atomic_writer, now, sha256, write_json
from prepare_kegg_reference import LINK_URLS, prepare
from verify_kegg_reference import verify

DOWNLOADS = {
    "profiles": ("https://www.genome.jp/ftp/db/kofam/profiles.tar.gz", "profiles.tar.gz"),
    "ko_list": ("https://www.genome.jp/ftp/db/kofam/ko_list.gz", "ko_list.gz"),
    "module_links": (LINK_URLS["module"], "ko_module_links.tsv"),
    "pathway_links": (LINK_URLS["pathway"], "ko_pathway_links.tsv"),
}


def download_cached(url, path):
    path = Path(path)
    receipt = path.with_name(path.name + ".json")
    if receipt.is_file() and path.is_file():
        record = json.loads(receipt.read_text())
        if (record.get("url") != url or record.get("bytes") != path.stat().st_size
                or record.get("sha256") != sha256(path)):
            raise ValueError(f"cached download checksum/source mismatch: {path}")
        print(f"Reusing download: {path}", flush=True)
        return record
    print(f"Downloading {url}", flush=True)
    time.sleep(0.5)
    request = Request(url, headers={"User-Agent": "phenoradar_prep/KEGG-reference"})
    with urlopen(request, timeout=60) as response, atomic_writer(path, "wb") as handle:
        if response.status != 200:
            raise ValueError(f"unexpected HTTP status {response.status} for {url}")
        digest, size = hashlib.sha256(), 0
        for block in iter(lambda: response.read(1024 * 1024), b""):
            handle.write(block)
            digest.update(block)
            size += len(block)
        length = response.headers.get("Content-Length")
        if size == 0 or (length is not None and size != int(length)):
            raise ValueError(f"empty or incomplete download: {url}")
        record = {"url": url, "retrieved_at": now(), "bytes": size, "sha256": digest.hexdigest(),
                  "last_modified": response.headers.get("Last-Modified"), "etag": response.headers.get("ETag")}
    write_json(receipt, record)
    return record


def extract_profiles(archive, destination):
    # Stream extraction and reject links/devices; validate the entire gzip stream.
    with gzip.open(archive, "rb") as compressed:
        with tarfile.open(fileobj=compressed, mode="r|") as handle:
            for member in handle:
                if not (member.isfile() or member.isdir()):
                    raise ValueError(f"unsupported KOfam archive member: {member.name}")
                handle.extract(member, destination, filter="data")
        while compressed.read(1024 * 1024):
            pass


def bootstrap(reference_dir):
    root = Path(reference_dir).absolute()
    root.parent.mkdir(parents=True, exist_ok=True)
    with (root.parent / f".{root.name}.bootstrap.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if root.exists() or root.is_symlink():
            verify(root / "reference.json")
            print(f"Keeping existing reference snapshot: {root}", flush=True)
            return root / "reference.json"
        cache = root.parent / "downloads"
        records = {key: download_cached(url, cache / filename)
                   for key, (url, filename) in DOWNLOADS.items()}
        with tempfile.TemporaryDirectory(prefix=f".{root.name}.extracting-", dir=root.parent) as temporary:
            staging = Path(temporary)
            extract_profiles(cache / "profiles.tar.gz", staging / "profiles")
            with gzip.open(cache / "ko_list.gz", "rb") as source, (staging / "ko_list").open("wb") as output:
                shutil.copyfileobj(source, output)
            return prepare(staging / "profiles", staging / "ko_list", root,
                           cache / "ko_module_links.tsv", cache / "ko_pathway_links.tsv",
                           release="downloaded-" + records["profiles"]["retrieved_at"][:10],
                           source_downloads=records)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", default="resources/kegg/snapshot_v1")
    print(bootstrap(**vars(parser.parse_args())))
