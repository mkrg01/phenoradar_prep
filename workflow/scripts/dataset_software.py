"""Acquire pinned GeneGalleon source and SIF files before dataset jobs are submitted."""
import copy
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import urllib.request
from pathlib import Path

from common import now, sha256, write_json
from dataset_assets import digest, locked, record, verify

SOURCE_BASE = "https://codeload.github.com/kfuku52/genegalleon/tar.gz/"
ENTRYPOINT = "workflow/gg_transcriptome_generation_entrypoint.sh"
KEYS = {"repository", "image", "settings", "version", "revision", "image_uri", "image_sha256", "cache_dir"}


def validate(config):
    unknown = set(config) - KEYS
    if unknown:
        raise ValueError(f"unknown genegalleon settings: {sorted(unknown)}")
    if config.get("version") is not None and not re.fullmatch(r"[0-9]+\.[0-9]+\.[0-9]+", str(config["version"])):
        raise ValueError("genegalleon.version must be an explicit major.minor.patch version")
    if config.get("revision") is not None and not re.fullmatch(r"[0-9a-f]{40}", str(config["revision"])):
        raise ValueError("genegalleon.revision must be a full 40-character commit SHA")
    uri = config.get("image_uri")
    if uri is not None:
        if not isinstance(uri, str) or not (re.fullmatch(r"docker://[A-Za-z0-9./:_-]+@sha256:[0-9a-f]{64}", uri)
                                            or uri.startswith("https://")):
            raise ValueError("genegalleon.image_uri must use an OCI SHA256 digest or an HTTPS SIF URL")
        if uri.startswith("https://") and not config.get("image_sha256"):
            raise ValueError("HTTPS SIF downloads require genegalleon.image_sha256")
    if config.get("image_sha256") is not None and not re.fullmatch(r"[0-9a-f]{64}", str(config["image_sha256"])):
        raise ValueError("genegalleon.image_sha256 must be a SHA256 checksum")


def progress(message):
    print(f"dataset software: {message}", file=sys.stderr, flush=True)


def download(url, target):
    request = urllib.request.Request(url, headers={"User-Agent": "phenoradar-prep-dataset"})
    with urllib.request.urlopen(request, timeout=30) as response, Path(target).open("wb") as output:
        if not response.geturl().startswith("https://"):
            raise ValueError("software downloads must remain on HTTPS")
        shutil.copyfileobj(response, output, length=8 * 1024 * 1024)


def source_records(repository):
    repository = Path(repository)
    if not (repository / ENTRYPOINT).is_file():
        raise ValueError(f"GeneGalleon transcriptome entrypoint not found: {repository}")
    runtime_dirs = {"__pycache__", ".pytest_cache", ".mypy_cache", ".snakemake"}
    paths = sorted(p for p in (repository / "workflow").rglob("*")
                   if p.is_file() and not runtime_dirs.intersection(p.parts) and p.suffix not in {".pyc", ".pyo"})
    if (repository / "VERSION").is_file(): paths.append(repository / "VERSION")
    return [record(p) for p in paths]


def read_receipt(destination, identity):
    receipt = json.loads((destination / "receipt.json").read_text())
    if receipt.get("schema_version") != 1 or receipt.get("identity") != identity:
        raise ValueError(f"software cache identity differs: {destination}")
    for entry in receipt["files"]: verify(entry)
    return receipt


def relocate(records, old, new):
    return [dict(entry, path=str(new / Path(entry["path"]).relative_to(old))) for entry in records]


def fetch_source(config):
    revision, version = config.get("revision"), config.get("version")
    if not revision or not version:
        raise ValueError("automatic GeneGalleon source retrieval requires version and revision")
    base = Path(config["cache_dir"]) / "source"
    destination = base / revision
    identity = {"revision": revision, "version": version, "url": SOURCE_BASE + revision}
    with locked(base / (revision + ".lock")):
        if destination.exists():
            receipt = read_receipt(destination, identity)
            return destination / "repository", receipt
        staging = Path(tempfile.mkdtemp(prefix=".fetch-", dir=base))
        try:
            progress(f"fetching GeneGalleon {version} source at {revision}")
            archive = staging / "source.tar.gz"
            download(identity["url"], archive)
            extracted = staging / "unpacked"
            with tarfile.open(archive, "r:gz") as bundle:
                members = bundle.getmembers()
                prefix = "genegalleon-" + revision
                if not members or any(Path(m.name).parts[0] != prefix for m in members):
                    raise ValueError("unexpected GeneGalleon source archive layout")
                # Explicit data filter rejects paths/links escaping the extraction directory.
                try:
                    bundle.extractall(extracted, filter="data")
                except tarfile.TarError as error:
                    raise ValueError(f"unsafe GeneGalleon source archive: {error}") from error
            repository = staging / "repository"
            (extracted / prefix).rename(repository)
            actual = (repository / "VERSION").read_text().strip()
            if actual != version:
                raise ValueError(f"GeneGalleon version mismatch: requested {version}, source contains {actual}")
            files = source_records(repository)
            receipt = {"schema_version": 1, "kind": "downloaded_source", "created_at": now(),
                       "identity": identity, "archive_sha256": sha256(archive),
                       "files": relocate(files, staging, destination)}
            archive.unlink()
            shutil.rmtree(extracted)
            write_json(staging / "receipt.json", receipt)
            os.rename(staging, destination)
            return destination / "repository", receipt
        finally:
            if staging.exists(): shutil.rmtree(staging)


def architecture():
    value = platform.machine().lower()
    try:
        return {"x86_64": "amd64", "amd64": "amd64", "aarch64": "arm64", "arm64": "arm64"}[value]
    except KeyError as error:
        raise ValueError(f"unsupported GeneGalleon container architecture: {value}") from error


def inspect_image(runtime, path, config):
    data = json.loads(subprocess.check_output([runtime, "inspect", "--json", str(path)], text=True))
    labels = data.get("data", {}).get("attributes", {}).get("labels", {})
    # OCI labels are retained by standard Apptainer pulls. Direct SIF downloads
    # are additionally bound to the configured file checksum.
    check_image_labels(labels, config)
    return labels


def fetch_image(config):
    uri = config.get("image_uri")
    if not uri:
        raise ValueError("automatic GeneGalleon SIF retrieval requires image_uri")
    identity = {"uri": uri, "architecture": architecture(), "expected_sha256": config.get("image_sha256")}
    base = Path(config["cache_dir"]) / "images"
    key = digest(identity)
    destination = base / key
    with locked(base / (key + ".lock")):
        if destination.exists():
            receipt = read_receipt(destination, identity)
            check_image_labels(receipt["labels"], config)
            return destination / "genegalleon.sif", receipt
        runtime = shutil.which("apptainer") or shutil.which("singularity")
        if not runtime:
            raise ValueError("install Apptainer or Singularity to fetch/validate the GeneGalleon SIF")
        staging = Path(tempfile.mkdtemp(prefix=".fetch-", dir=base))
        try:
            image = staging / "genegalleon.sif"
            progress(f"fetching {uri} for {identity['architecture']} (first OCI pull converts to SIF)")
            if uri.startswith("https://"):
                download(uri, image)
            else:
                runtime_cache = Path(config["cache_dir"]) / "runtime-cache"
                runtime_cache.mkdir(parents=True, exist_ok=True)
                scratch = staging / "tmp"
                scratch.mkdir()
                env = dict(os.environ)
                for prefix in ("APPTAINER", "SINGULARITY"):
                    env[prefix + "_CACHEDIR"] = str(runtime_cache)
                    env[prefix + "_TMPDIR"] = str(scratch)
                env["TMPDIR"] = str(scratch)
                subprocess.run([runtime, "pull", "--arch", identity["architecture"], str(image), uri], env=env, check=True)
                shutil.rmtree(scratch)
            entry = record(image)
            if not entry["bytes"]: raise ValueError("downloaded GeneGalleon SIF is empty")
            if config.get("image_sha256") and entry["sha256"] != config["image_sha256"]:
                raise ValueError("GeneGalleon SIF SHA256 mismatch")
            labels = inspect_image(runtime, image, config)
            receipt = {"schema_version": 1, "kind": "downloaded_image", "created_at": now(),
                       "identity": identity, "runtime": subprocess.check_output([runtime, "--version"], text=True).strip(),
                       "labels": labels, "files": relocate([entry], staging, destination)}
            write_json(staging / "receipt.json", receipt)
            os.rename(staging, destination)
            return destination / "genegalleon.sif", receipt
        finally:
            if staging.exists(): shutil.rmtree(staging)


def check_image_labels(labels, config):
    for key, setting in (("org.opencontainers.image.revision", "revision"),
                         ("org.opencontainers.image.version", "version")):
        if labels.get(key) and config.get(setting) and labels[key] != config[setting]:
            raise ValueError(f"GeneGalleon image/source mismatch: {key}={labels[key]}, expected {config[setting]}")


def resolve(config):
    """Return concrete paths, immutable file receipts, and the resolved software lock."""
    validate(config)
    cfg = copy.deepcopy(config)
    if cfg.get("repository"):
        repository = Path(cfg["repository"])
        source = {"kind": "local_source", "files": source_records(repository)}
        if (repository / "VERSION").is_file():
            source["version"] = (repository / "VERSION").read_text().strip()
    else:
        repository, source = fetch_source(cfg)
    local_image = cfg.get("image")
    if not local_image and cfg.get("repository") and (repository / "genegalleon.sif").is_file():
        local_image = repository / "genegalleon.sif"
    if local_image:
        image = Path(local_image)
        entry = record(image)
        if cfg.get("image_sha256") and entry["sha256"] != cfg["image_sha256"]:
            raise ValueError("local GeneGalleon SIF SHA256 mismatch")
        container = {"kind": "local_image", "files": [entry]}
    else:
        if source["kind"] == "local_source" and source.get("version") != cfg.get("version"):
            raise ValueError("local GeneGalleon version differs from pinned image; provide its matching local image or update pins")
        image, container = fetch_image(cfg)
    cfg["repository"], cfg["image"] = str(repository), str(image)
    lock = {"schema_version": 1, "source": source, "container": container}
    return cfg, [*source["files"], *container["files"]], lock
