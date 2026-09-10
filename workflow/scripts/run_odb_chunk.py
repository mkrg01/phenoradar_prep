#!/usr/bin/env python3
"""Run one ODB chunk with content-based isolation of resumable work."""
import argparse
import fcntl
import hashlib
import json
import os
import re
import shutil
import subprocess
import tempfile
import uuid
from pathlib import Path

from common import atomic_writer, file_record, now, sha256, write_json
from odb_environment import check_storage, command_path, odb_environment, software_records


def publish_tree(source, destination):
    """Copy fully before replacing a prior generated tree; never mix old files."""
    destination = Path(destination)
    staging = Path(tempfile.mkdtemp(prefix=f".{destination.name}.", dir=destination.parent))
    try:
        shutil.copytree(source, staging / "new", symlinks=False)
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


def run(manifest, reference, output_dir, work_dir, label, command="ODB-mapper_v12", prefix="",
        version="v12", node=3193, jobs=16, batch_size=64, min_free_gb=750,
        allow_nonlocal=False, keep_work=False):
    if not re.fullmatch(r"chunk_[0-9]+", label) or jobs < 1 or batch_size < jobs:
        raise ValueError("invalid chunk label or concurrency settings")
    if jobs > int(os.environ.get("SLURM_CPUS_PER_TASK", jobs)):
        raise ValueError("ODB workers exceed Slurm CPU allocation")
    ref = json.loads(Path(reference).read_text())
    if ref["version"] != version or ref["node"] != node:
        raise ValueError("reference node/version differs from requested mapping")
    data = Path(ref["data_dir"]).resolve()
    inventory_path = Path(ref["inventory"]["path"])
    if not data.is_dir() or sha256(inventory_path) != ref["inventory"]["sha256"]:
        raise ValueError("reference data/inventory missing or modified")
    # Verify the recorded files still exist and have their recorded sizes.
    # Full SHA256 verification is available separately; hashing TB per chunk is avoided.
    for entry in json.loads(inventory_path.read_text()):
        path = data / entry["relative_path"]
        if not path.is_file() or path.stat().st_size != entry["bytes"]:
            raise ValueError(f"reference file missing or size changed: {path}")
    proteins = Path(manifest).read_text().splitlines()
    if not proteins or any(not p for p in proteins) or len(set(proteins)) != len(proteins):
        raise ValueError("manifest must contain unique, nonempty FASTA paths")
    env = odb_environment(prefix, version)
    command = command_path(command, env)
    # Upstream ODB shell code does not support whitespace/shell metacharacters in paths.
    for path in [manifest, reference, output_dir, work_dir, command, str(data), *proteins]:
        if not re.fullmatch(r"[A-Za-z0-9_./+-]+", str(Path(path).absolute())):
            raise ValueError(f"ODB path contains unsupported characters: {path}")
    scripts = Path(__file__).resolve().parent
    identity = {"manifest": file_record(manifest), "proteins": [file_record(p) for p in proteins],
                "reference": file_record(reference), "node": node, "version": version,
                "jobs": jobs, "batch_size": batch_size, "software": software_records(command, prefix),
                "implementation": [file_record(scripts / name) for name in
                                   ["run_odb_chunk.py", "odb_map.sh", "odb_environment.py", "common.py"]]}
    fingerprint = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    base = Path(work_dir).resolve() / label
    check_storage(base, min_free_gb, allow_nonlocal)
    with open(base / ".lock", "w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        work = base / fingerprint
        # Only failed/incomplete work is reusable. A forced rerun after success is fresh.
        if (work / "completed.json").exists():
            work.rename(base / f"{fingerprint}.completed.{uuid.uuid4().hex[:8]}")
        work.mkdir(parents=True, exist_ok=True)
        write_json(work / "identity.json", identity)
        status = base / "status.json"
        write_json(status, {"state": "running", "started_at": now(), "fingerprint": fingerprint, "work": str(work)})
        try:
            subprocess.run(["bash", str(scripts / "odb_map.sh"), command, str(Path(manifest).resolve()),
                            label, str(node), version, str(data), str(work), str(jobs), str(batch_size)],
                           env=env, check=True)
            project = work / "odbmapper" / version / "pipeline"
            required = [project / "Results" / f"{label}.{suffix}" for suffix in
                        ["og.annotations", "og.hits", "summary.txt"]]
            for path in required:
                if not path.is_file() or not path.stat().st_size:
                    raise ValueError(f"ODB result missing or empty: {path}")
            out = Path(output_dir)
            out.mkdir(parents=True, exist_ok=True)
            # Publish validated files individually and write provenance last.
            for path in required + [project / "orthologer_conf.sh", work / "report.txt", work / "odbmapper_config.txt"]:
                with open(path, "rb") as src, atomic_writer(out / path.name, "wb") as dst:
                    shutil.copyfileobj(src, dst)
            # Retain the other native outputs and logs before removing successful work.
            publish_tree(project / "Results", out / "native_results")
            if (project / "RunLogs").exists():
                publish_tree(project / "RunLogs", out / "RunLogs")
            record = {"completed_at": now(), "fingerprint": fingerprint, "identity": identity,
                      "results": [file_record(out / p.name) for p in required]}
            write_json(out / "provenance.json", record)
            write_json(work / "completed.json", record)
            write_json(status, {"state": "success", "completed_at": now(), "fingerprint": fingerprint})
            if not keep_work:
                # The only removed directory is our just-validated fingerprint work.
                assert work.parent == base and work.name == fingerprint
                shutil.rmtree(work)
        except BaseException:
            write_json(status, {"state": "failed", "time": now(), "fingerprint": fingerprint, "work": str(work)})
            raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ["manifest", "reference", "output-dir", "work-dir", "label"]:
        parser.add_argument(f"--{flag}", required=True)
    parser.add_argument("--command", default="ODB-mapper_v12")
    parser.add_argument("--prefix", default="")
    parser.add_argument("--version", default="v12")
    parser.add_argument("--node", type=int, default=3193)
    parser.add_argument("--jobs", type=int, default=16)
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--min-free-gb", type=float, default=750)
    parser.add_argument("--allow-nonlocal", action="store_true")
    parser.add_argument("--keep-work", action="store_true")
    run(**vars(parser.parse_args()))
