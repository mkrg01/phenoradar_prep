#!/usr/bin/env python3
"""Initialize a new OrthoDB snapshot; completed snapshots are never updated."""
import argparse
import fcntl
import json
import re
import subprocess
from pathlib import Path

from common import atomic_writer, file_record, now, sha256, write_json
from odb_environment import command_path, odb_environment, software_records
from odb_settings import pin_reference, validate as validate_odb
from check_odb_runtime import check as check_runtime


def _initialize(reference_dir, command="ODB-mapper", prefix="", version="v12", node=3193):
    root = Path(reference_dir).resolve()
    validate_odb({'version': version})
    if type(node) is not int or node < 1:
        raise ValueError("OrthoDB requires a positive node ID")
    if not re.fullmatch(r"[A-Za-z0-9_./+-]+", str(root)):
        raise ValueError(f"ODB path contains unsupported characters: {root}")
    root.mkdir(parents=True, exist_ok=True)
    with open(root / ".prepare.lock", "w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        marker = root / "reference.json"
        if marker.exists():
            record = json.loads(marker.read_text())
            if record["version"] != version or record["node"] != node:
                raise ValueError("reference settings changed; use a new reference directory")
            if not Path(record["data_dir"]).is_dir():
                raise ValueError("reference marker exists but its data directory is missing")
            print(f"Keeping existing reference snapshot: {marker}")
            return
        env = odb_environment(prefix, version)
        command = command_path(command, env)
        # Check new downloads without adding run-specific dependencies to a frozen snapshot.
        check_runtime(command=command, prefix=prefix, version=version, work_dir=root)
        settings = {"node": node, "version": version, "command": command}
        setup = root / "setup_identity.json"
        if setup.exists() and json.loads(setup.read_text()) != settings:
            raise ValueError("incomplete reference has different settings; use a new reference directory")
        write_json(setup, settings)
        env["ODBMAPPER_WORK"] = str(root / "odbmapper")
        subprocess.run([command, 'SETUP'], env=env, check=True)
        api = subprocess.check_output([command, "CONFIG", "apiver"], env=env, text=True).strip()
        if not re.fullmatch(re.escape(version) + r'(?:\.\d+)*', api):
            raise ValueError(f'ODB-mapper selected {api!r}, requested {version}; install a mapper supporting this release')
        subprocess.run([command, 'DOWNLOAD', str(node)], env=env, check=True)
        data = root / "odbmapper" / version / "data"
        files = sorted(p for p in data.rglob("*") if p.is_file())
        if sum(p.stat().st_size for p in files) < 1024**2:
            raise ValueError(f"downloaded reference looks incomplete: {data}")
        for name, args in [("dbinfo.txt", ["DBINFO"]), ("config.txt", ["CONFIG"])]:
            with atomic_writer(root / name) as handle:
                subprocess.run([command, *args], env=env, stdout=handle, check=True)
        print(f"Recording SHA256 checksums for {len(files)} reference files", flush=True)
        inventory = [{"relative_path": str(p.relative_to(data)), **file_record(p)} for p in files]
        write_json(root / "files.json", inventory)
        write_json(marker, {"created_at": now(), "version": version, "api_version": api, "node": node,
                            "data_dir": str(data), "inventory": file_record(root / "files.json"),
                            "software": software_records(command, prefix)})


def prepare(reference_dir, command="ODB-mapper", prefix="", version="v12", node=3193,
            receipt=None, reference_pin=None):
    """Validate the immutable snapshot before pinning it for this build.

    Snakemake owns only the per-run receipt, so reruns cannot remove imported
    reference files before the existing-snapshot check executes.
    """
    root = Path(reference_dir).resolve()
    validate_odb({'version': version})
    if type(node) is not int or node < 1:
        raise ValueError("OrthoDB requires a positive node ID")
    marker = root / 'reference.json'
    if reference_pin:
        pinned = pin_reference(reference_pin, version, node)
        if pinned['reference_sha256']:
            if not marker.is_file():
                raise ValueError('pinned reference is missing; restore the original reference')
            pin_reference(reference_pin, version, node, sha256(marker))
    _initialize(root, command, prefix, version, node)
    reference = file_record(marker)
    data = json.loads(marker.read_text())
    inventory = Path(data['inventory']['path'])
    if not inventory.is_file() or sha256(inventory) != data['inventory']['sha256']:
        raise ValueError('reference inventory missing or modified')
    entries = json.loads(inventory.read_text())
    if not entries:
        raise ValueError('reference inventory is empty')
    directory = Path(data['data_dir']).resolve()
    for entry in entries:
        path = (directory / entry['relative_path']).resolve()
        if not path.is_relative_to(directory):
            raise ValueError(f'reference file escapes data directory: {path}')
        if not path.is_file() or path.stat().st_size != entry['bytes']:
            raise ValueError(f'reference file missing or size changed: {path}')
    if sha256(marker) != reference['sha256']:
        raise ValueError('reference marker changed during validation')
    if reference_pin:
        pin_reference(reference_pin, version, node, reference['sha256'])
    result = {'schema_version': 1, 'checked_at': now(), 'version': version, 'node': node,
              'reference': reference, 'inventory': file_record(inventory), 'files_checked': len(entries)}
    if receipt:
        write_json(receipt, result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", required=True)
    parser.add_argument("--command", default="ODB-mapper")
    parser.add_argument("--prefix", default="")
    parser.add_argument("--version", default="v12")
    parser.add_argument("--node", type=int, default=3193)
    parser.add_argument("--receipt")
    parser.add_argument("--reference-pin")
    prepare(**vars(parser.parse_args()))
