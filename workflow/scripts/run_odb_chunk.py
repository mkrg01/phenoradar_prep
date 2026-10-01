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
from pathlib import Path

from common import atomic_writer, file_record, now, sha256, write_json
from odb_environment import command_path, odb_environment, software_records
from cleanup_work import cleanup


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


def completed_output(output, fingerprint, label):
    """Reuse verified published files even after internal work was removed."""
    output = Path(output)
    try:
        record = json.loads((output / 'provenance.json').read_text())
        expected = {f'{label}.{suffix}' for suffix in ('og.annotations', 'og.hits', 'summary.txt')}
        results = record['results']
        if (record['fingerprint'] != fingerprint or len(results) != len(expected)
                or {Path(entry['path']).name for entry in results} != expected):
            return None
        for entry in results + record.get('artifacts', []):
            path = Path(entry['path'])
            if not path.is_relative_to(output) or path.is_symlink() or file_record(path) != entry:
                return None
        return record
    except (OSError, ValueError, KeyError, TypeError):
        return None


def run(manifest, reference, output_dir, work_dir, label, command="ODB-mapper_v12", prefix="",
        version="v12", node=3193, jobs=16, *, batch_size, keep_intermediates=False, defer_cleanup=False):
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
                                   ["run_odb_chunk.py", "odb_map.sh", "odb_environment.py", "common.py", "cleanup_work.py"]]}
    fingerprint = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    out, work_root = Path(output_dir).absolute(), Path(work_dir).absolute()
    if out.is_symlink() or work_root.is_symlink():
        raise ValueError('ODB work and output must not be symlinks')
    out, work_root = out.resolve(), work_root.resolve()
    if out == work_root or out.is_relative_to(work_root) or work_root.is_relative_to(out):
        raise ValueError('ODB work and output directories must not overlap')
    base = work_root / label
    if base.is_symlink() or (base / fingerprint).is_symlink():
        raise ValueError('ODB scratch must not be a symlink')
    work = base / fingerprint
    if any(Path(p).resolve().is_relative_to(work) for p in [manifest, reference, *proteins]):
        raise ValueError('ODB scratch must not contain inputs')
    base.mkdir(parents=True, exist_ok=True)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(base / ".lock", "a") as lock, open(out.parent / f'.{out.name}.lock', 'a') as out_lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        fcntl.flock(out_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        status = base / "status.json"
        cached = completed_output(out, fingerprint, label)
        if cached is not None:
            write_json(status, {"state": "success", "completed_at": now(),
                                "fingerprint": fingerprint, "work": str(work), "reused": True})
            if not keep_intermediates and not defer_cleanup:
                cleanup(base, [fingerprint], base / 'cleanup.json')
            return cached
        # Failed work remains available for ODB's internal resume.
        work.mkdir(parents=True, exist_ok=True)
        (work / "completed.json").unlink(missing_ok=True)
        write_json(work / "identity.json", identity)
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
            out.mkdir(parents=True, exist_ok=True)
            # Publish validated files individually and write provenance last.
            published = required + [project / "orthologer_conf.sh", work / "report.txt", work / "odbmapper_config.txt"]
            for path in published:
                source_record = file_record(path)
                with open(path, "rb") as src, atomic_writer(out / path.name, "wb") as dst:
                    shutil.copyfileobj(src, dst)
                saved = file_record(out / path.name)
                if any(saved[key] != source_record[key] for key in ('bytes', 'sha256')):
                    raise ValueError(f'ODB publication differs from source: {path}')
            # Native exports are diagnostic; standard results suffice downstream.
            if keep_intermediates:
                publish_tree(project / "Results", out / "native_results")
            if (project / "RunLogs").exists():
                publish_tree(project / "RunLogs", out / "RunLogs")
            record = {"completed_at": now(), "fingerprint": fingerprint, "identity": identity,
                      "results": [file_record(out / p.name) for p in required],
                      "artifacts": [file_record(out / p.name) for p in published[len(required):]] +
                                   [file_record(p) for p in sorted((out / 'RunLogs').rglob('*')) if p.is_file()]}
            write_json(out / "provenance.json", record)
            if completed_output(out, fingerprint, label) is None:
                raise ValueError("ODB published results failed verification")
            write_json(work / "completed.json", record)
            write_json(status, {"state": "success", "completed_at": now(),
                                "fingerprint": fingerprint, "work": str(work)})
        except BaseException:
            write_json(status, {"state": "failed", "time": now(), "fingerprint": fingerprint, "work": str(work)})
            raise
        if not keep_intermediates and not defer_cleanup:
            cleanup(base, [fingerprint], base / "cleanup.json")
        return record


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ["manifest", "reference", "output-dir", "work-dir", "label"]:
        parser.add_argument(f"--{flag}", required=True)
    parser.add_argument("--command", default="ODB-mapper_v12")
    parser.add_argument("--prefix", default="")
    parser.add_argument("--version", default="v12")
    parser.add_argument("--node", type=int, default=3193)
    parser.add_argument("--jobs", type=int, default=16)
    parser.add_argument("--batch-size", type=int, required=True)
    parser.add_argument("--keep-intermediates", action="store_true")
    parser.add_argument("--defer-cleanup", action="store_true", help="Retain scratch until the workflow publishes its reuse cache")
    run(**vars(parser.parse_args()))
