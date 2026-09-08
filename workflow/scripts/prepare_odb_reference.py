#!/usr/bin/env python3
"""Initialize a new OrthoDB snapshot; completed snapshots are never updated."""
import argparse
import fcntl
import json
import re
import subprocess
from pathlib import Path

from common import atomic_writer, file_record, now, write_json
from odb_environment import check_storage, command_path, odb_environment, software_records


def prepare(reference_dir, command="ODB-mapper_v12", prefix="", version="v12", node=3193,
            min_free_gb=200, allow_nonlocal=False):
    root = Path(reference_dir).resolve()
    if version != "v12" or node < 1:
        raise ValueError("this workflow supports ODB v12 with a positive node ID")
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
        check_storage(root, min_free_gb, allow_nonlocal)
        env = odb_environment(prefix, version)
        command = command_path(command, env)
        settings = {"node": node, "version": version, "command": command}
        setup = root / "setup_identity.json"
        if setup.exists() and json.loads(setup.read_text()) != settings:
            raise ValueError("incomplete reference has different settings; use a new reference directory")
        write_json(setup, settings)
        env["ODBMAPPER_WORK"] = str(root / "odbmapper")
        for args in [["SETUP"], ["DOWNLOAD", str(node)]]:
            subprocess.run([command, *args], env=env, check=True)
        data = root / "odbmapper" / version / "data"
        files = sorted(p for p in data.rglob("*") if p.is_file())
        if sum(p.stat().st_size for p in files) < 1024**2:
            raise ValueError(f"downloaded reference looks incomplete: {data}")
        for name, args in [("dbinfo.txt", ["DBINFO"]), ("config.txt", ["CONFIG"])]:
            with atomic_writer(root / name) as handle:
                subprocess.run([command, *args], env=env, stdout=handle, check=True)
        api = subprocess.check_output([command, "CONFIG", "apiver"], env=env, text=True).strip()
        print(f"Recording SHA256 checksums for {len(files)} reference files", flush=True)
        inventory = [{"relative_path": str(p.relative_to(data)), **file_record(p)} for p in files]
        write_json(root / "files.json", inventory)
        write_json(marker, {"created_at": now(), "version": version, "api_version": api, "node": node,
                            "data_dir": str(data), "inventory": file_record(root / "files.json"),
                            "software": software_records(command, prefix)})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", required=True)
    parser.add_argument("--command", default="ODB-mapper_v12")
    parser.add_argument("--prefix", default="")
    parser.add_argument("--version", default="v12")
    parser.add_argument("--node", type=int, default=3193)
    parser.add_argument("--min-free-gb", type=float, default=200)
    parser.add_argument("--allow-nonlocal", action="store_true")
    prepare(**vars(parser.parse_args()))
