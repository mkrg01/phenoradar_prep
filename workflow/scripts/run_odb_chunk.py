#!/usr/bin/env python3
"""Run one ODB chunk with content-based isolation of resumable work."""
import argparse
import fcntl
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import tempfile
import tarfile
from pathlib import Path

from common import atomic_writer, file_record, now, sha256, write_json
from odb_environment import command_path, odb_environment, software_records
from cleanup_work import cleanup
from product_storage import decoded_sha256, stage_product, product_record
from odb_settings import pin_reference


def stage_reference_data(source, destination):
    """Keep mapper-writable todo lists private, without copying large databases.

    Orthologer check_todo sorts lists in place (cp -f follows symlinks).
    Real directories also keep newly generated lists out of the shared snapshot.
    """
    source, destination = Path(source).resolve(), Path(destination)
    if destination.is_symlink():
        raise ValueError('ODB reference workspace must not be a symlink')
    destination.mkdir(parents=True, exist_ok=True)
    for entry in sorted(source.rglob('*')):
        target = destination / entry.relative_to(source)
        if entry.is_dir():
            if target.is_symlink():
                raise ValueError(f'ODB reference workspace directory is a symlink: {target}')
            target.mkdir(exist_ok=True)
        elif entry.is_file():
            if entry.suffix == '.todo':
                if target.is_symlink() or (target.exists() and target.samefile(entry)):
                    raise ValueError(f'ODB todo workspace shares the frozen reference: {target}')
                if not target.exists():
                    with entry.open('rb') as src, atomic_writer(target, 'wb') as dst:
                        shutil.copyfileobj(src, dst)
            elif target.is_symlink():
                if target.resolve() != entry.resolve():
                    raise ValueError(f'ODB reference workspace link changed: {target}')
            elif target.exists():
                raise ValueError(f'ODB reference workspace has unexpected file: {target}')
            else:
                target.symlink_to(entry.resolve())
    return destination


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
                or {Path(entry['path']).name.removesuffix('.gz') for entry in results} != expected):
            return None
        for entry in results + record.get('artifacts', []):
            path = Path(entry['path'])
            if not path.is_relative_to(output) or path.is_symlink():
                return None
            actual = file_record(path)
            if any(actual[key] != entry[key] for key in ('path', 'bytes', 'sha256')):
                return None
            if path.suffix == '.gz' and entry in results and (
                    entry.get('content_sha256') != decoded_sha256(path)):
                return None
        return record
    except (OSError, ValueError, KeyError, TypeError):
        return None


def compress_output(output, fingerprint, label):
    """Migrate a verified publication; commit its receipt before retiring plain files.

    The caller holds the job/output locks. Scientific content identities survive
    encoding changes, and a crash before the receipt commit leaves plain reuse valid.
    """
    output = Path(output)
    saved = completed_output(output, fingerprint, label)
    if saved is None:
        raise ValueError('cannot compress unverified ODB results')
    result, retired = [], []
    for entry in saved['results']:
        path = Path(entry['path'])
        if path.name.endswith(('.og.annotations', '.og.hits')):
            result.append(stage_product(entry, Path(str(path) + '.gz')))
            retired.append(path)
        else:
            result.append(entry)
    if retired:
        saved = dict(saved, results=result)
        write_json(output / 'provenance.json', saved)
        if completed_output(output, fingerprint, label) is None:
            raise ValueError('compressed ODB publication failed verification')
        for path in retired:
            path.unlink()
    return saved


def run(manifest, reference, output_dir, work_dir, label, command="ODB-mapper", prefix="",
        version="v12", node=3193, jobs=16, *, batch_size, keep_intermediates=False, defer_cleanup=False,
        reference_pin=None, gzip_results=False):
    if not re.fullmatch(r"chunk_[0-9]+", label) or jobs < 1 or batch_size < jobs:
        raise ValueError("invalid chunk label or concurrency settings")
    if jobs > int(os.environ.get("SLURM_CPUS_PER_TASK", jobs)):
        raise ValueError("ODB workers exceed Slurm CPU allocation")
    ref = json.loads(Path(reference).read_text())
    if ref["version"] != version or ref["node"] != node:
        raise ValueError("reference node/version differs from requested mapping")
    reference_entry = file_record(reference)
    if reference_pin:
        pin_reference(reference_pin, version, node, reference_entry['sha256'])
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
    protein_entries = [dict(file_record(p), content_sha256=decoded_sha256(p)) if p.endswith('.gz')
                       else file_record(p) for p in proteins]
    identity = {"manifest": file_record(manifest), "proteins": protein_entries,
                "reference": reference_entry, "node": node, "version": version,
                "jobs": jobs, "batch_size": batch_size, "software": software_records(command, prefix),
                "implementation": [file_record(scripts / name) for name in
                                   ["run_odb_chunk.py", "odb_map.sh", "odb_environment.py", "odb_settings.py", "common.py", "cleanup_work.py", "product_storage.py"]]}
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
    persistent_work = work
    if os.environ.get('PHENORADAR_TASK_SCRATCH'):
        local = Path(os.environ['PHENORADAR_TASK_SCRATCH'])
        if not local.is_absolute() or local.is_symlink() or not local.is_dir():
            raise ValueError('ODB task scratch must be an existing absolute directory')
        work = local / 'odb' / label / fingerprint
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
            if gzip_results:
                cached = compress_output(out, fingerprint, label)
            write_json(status, {"state": "success", "completed_at": now(),
                                "fingerprint": fingerprint, "work": str(persistent_work), "reused": True})
            if not keep_intermediates and not defer_cleanup:
                cleanup(base, [fingerprint], base / 'cleanup.json')
            return cached
        # Failed work remains available for ODB's internal resume.
        work.mkdir(parents=True, exist_ok=True)
        (work / "completed.json").unlink(missing_ok=True)
        write_json(work / "identity.json", identity)
        write_json(status, {"state": "running", "started_at": now(), "fingerprint": fingerprint,
                            "work": str(persistent_work), "local_work": str(work)})
        try:
            tool_manifest = Path(manifest).resolve()
            if any(p.endswith('.gz') for p in proteins):
                tool_manifest = work / 'proteins.fs'
                tool_inputs = []
                for entry in protein_entries:
                    path = Path(entry['path'])
                    if path.suffix == '.gz':
                        destination = work / 'inputs' / path.name.removesuffix('.gz')
                        stage_product(entry, destination)
                        path = destination
                    tool_inputs.append(str(path))
                with atomic_writer(tool_manifest) as handle:
                    handle.write('\n'.join(tool_inputs) + '\n')
            local_data = stage_reference_data(data, work / "reference-data")
            subprocess.run(["bash", str(scripts / "odb_map.sh"), command, str(tool_manifest),
                            label, str(node), version, str(local_data), str(work), str(jobs), str(batch_size)],
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
            result_records = []
            for path in published:
                compressed = gzip_results and path.name.endswith(('.og.annotations', '.og.hits'))
                target = out / (path.name + ('.gz' if compressed else ''))
                # Never hardlink node scratch into a durable publication.
                if compressed:
                    source_record = file_record(path)
                    with path.open('rb') as source, atomic_writer(target, 'wb') as raw:
                        with gzip.GzipFile(filename='', fileobj=raw, mode='wb', mtime=0, compresslevel=6) as zipped:
                            shutil.copyfileobj(source, zipped)
                    saved = product_record(target)
                    if saved['content_sha256'] != source_record['sha256']:
                        raise ValueError(f'compressed ODB publication differs from source: {path}')
                    result_records.append(saved)
                else:
                    source_record = file_record(path)
                    with open(path, 'rb') as src, atomic_writer(target, 'wb') as dst:
                        shutil.copyfileobj(src, dst)
                    saved = file_record(target)
                    if any(saved[key] != source_record[key] for key in ('bytes', 'sha256')):
                        raise ValueError(f'ODB publication differs from source: {path}')
                    if path in required:
                        result_records.append(saved)
            # Native exports are diagnostic; standard results suffice downstream.
            if keep_intermediates:
                publish_tree(project / "Results", out / "native_results")
            log_artifacts = []
            if (project / "RunLogs").exists():
                if os.environ.get('PHENORADAR_COMPACT_ODB_LOGS') == '1':
                    archive = out / 'RunLogs.tar.gz'
                    with atomic_writer(archive, 'wb') as handle:
                        with tarfile.open(fileobj=handle, mode='w:gz') as packed:
                            packed.add(project / 'RunLogs', arcname='RunLogs')
                    log_artifacts = [file_record(archive)]
                else:
                    publish_tree(project / "RunLogs", out / "RunLogs")
                    log_artifacts = [file_record(p) for p in sorted((out / 'RunLogs').rglob('*')) if p.is_file()]
            record = {"completed_at": now(), "fingerprint": fingerprint, "identity": identity,
                      "results": result_records,
                      "artifacts": [file_record(out / p.name) for p in published[len(required):]] +
                                   log_artifacts}
            write_json(out / "provenance.json", record)
            if completed_output(out, fingerprint, label) is None:
                raise ValueError("ODB published results failed verification")
            write_json(work / "completed.json", record)
            write_json(status, {"state": "success", "completed_at": now(),
                                "fingerprint": fingerprint, "work": str(persistent_work), "local_work": str(work)})
        except BaseException:
            write_json(status, {"state": "failed", "time": now(), "fingerprint": fingerprint,
                                "work": str(persistent_work), "local_work": str(work)})
            raise
        if work != persistent_work and not keep_intermediates:
            shutil.rmtree(work)
        if not keep_intermediates and not defer_cleanup:
            cleanup(base, [fingerprint], base / "cleanup.json")
        return record


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ["manifest", "reference", "output-dir", "work-dir", "label"]:
        parser.add_argument(f"--{flag}", required=True)
    parser.add_argument("--command", default="ODB-mapper")
    parser.add_argument("--reference-pin")
    parser.add_argument("--gzip-results", action="store_true", help="Publish lossless compressed annotations and hits")
    parser.add_argument("--prefix", default="")
    parser.add_argument("--version", default="v12")
    parser.add_argument("--node", type=int, default=3193)
    parser.add_argument("--jobs", type=int, default=16)
    parser.add_argument("--batch-size", type=int, required=True)
    parser.add_argument("--keep-intermediates", action="store_true")
    parser.add_argument("--defer-cleanup", action="store_true", help="Retain scratch until the workflow publishes its reuse cache")
    run(**vars(parser.parse_args()))
