#!/usr/bin/env python3
"""Build the pinned, unmodified official ASTRAL-IV int128 source."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import tarfile
import tempfile
from urllib.request import urlopen

from common import atomic_writer, file_record, now, sha256, write_json


SOURCES = {
    "aster": {"repository": "chaoszhang/ASTER", "commit": "ddc3dc6b28b20cecfec18cede8c6247decc82b77",
              "sha256": "ca4cb0f7160f8b1aaf73a4c09a4c6c3aba571b1ba42c7b64a108188c3a3c2beb", "binary": "astral4_int128"},
}


def prepare(destination, archives=None):
    root = Path(destination).resolve()
    (root / "bin").mkdir(parents=True, exist_ok=True)
    compiler = shutil.which("g++")
    if not compiler:
        raise ValueError("g++ is required to build the phylogeny tools")
    for name, spec in SOURCES.items():
        target = root / "bin" / spec["binary"]
        manifest = root / f"{name}.json"
        if manifest.is_file() and target.is_file():
            old = json.loads(manifest.read_text())
            if old.get("source") == spec and old["executable"]["sha256"] == sha256(target):
                print(f"Reusing {target}")
                continue
            raise ValueError(f"existing tool differs from pinned build: {target}; use a new destination")
        with tempfile.TemporaryDirectory(prefix=f".{name}-", dir=root) as temporary:
            tmp = Path(temporary)
            archive = Path(archives) / f"{name}.tar.gz" if archives else tmp / f"{name}.tar.gz"
            url = f'https://codeload.github.com/{spec["repository"]}/tar.gz/{spec["commit"]}'
            if not archives:
                with urlopen(url, timeout=60) as response, open(archive, "wb") as output:
                    shutil.copyfileobj(response, output)
            if sha256(archive) != spec["sha256"]:
                raise ValueError(f"source archive checksum mismatch: {archive}")
            with tarfile.open(archive) as handle:
                handle.extractall(tmp / "source", filter="data")
            source = next((tmp / "source").iterdir())
            built = tmp / spec["binary"]
            command = [compiler, "-std=gnu++17", "-O3", "-pthread", "-D", "LARGE_DATA",
                       str(source / "src/astral.cpp"), "-o", str(built)]
            subprocess.run(command, check=True)
            if not built.is_file():
                raise ValueError(f"build did not produce {built}")
            with atomic_writer(target, "wb") as output, open(built, "rb") as handle:
                shutil.copyfileobj(handle, output)
            os.chmod(target, 0o755)
            write_json(manifest, {"created_at": now(), "source": spec, "url": url,
                       "command": command, "compiler": subprocess.check_output([compiler, "--version"], text=True),
                       "executable": file_record(target),
                       "patch": None})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--destination", default="resources/phylogeny_tools")
    parser.add_argument("--archives", help="offline directory with verified aster.tar.gz")
    prepare(**vars(parser.parse_args()))
