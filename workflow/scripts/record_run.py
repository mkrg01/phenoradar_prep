#!/usr/bin/env python3
"""Record the resolved configuration, code checksums, and Python environment."""
import argparse
import json
import platform
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

from common import file_record, now, write_json


def record(config_json, selection, workflow_dir, output):
    root = Path(workflow_dir)
    packages = {}
    for name in ["pandas", "ete4", "matplotlib"]:
        try:
            packages[name] = version(name)
        except PackageNotFoundError:
            packages[name] = None
    code = sorted(p for p in root.rglob("*") if p.is_file() and
                  (p.suffix in {".py", ".sh", ".smk", ".yaml"} or p.name == "Snakefile"))
    write_json(output, {"created_at": now(), "config": json.loads(config_json),
                        "python": platform.python_version(), "packages": packages,
                        "selection": file_record(selection), "workflow_files": [file_record(p) for p in code]})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ["config-json", "selection", "workflow-dir", "output"]:
        parser.add_argument(f"--{flag}", required=True)
    record(**vars(parser.parse_args()))
