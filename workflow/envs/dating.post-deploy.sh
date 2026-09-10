#!/usr/bin/env bash
set -euo pipefail

# Self-contained: Snakemake copies and hashes this script with the environment.
# For offline deployment, set LSD2_SOURCE_ARCHIVE to the verified tar.gz file.
"${CONDA_PREFIX:?}/bin/python" - <<'PY'
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tarfile
import tempfile
from urllib.request import urlopen

commit = "c61110f3a4fa05325b45c97b2134792ff9d55d4c"
checksum = "9bbeaa0f8f35783c1d8dec74df6c93a804dbca808fa04484f9123de4e7258b53"
url = f"https://codeload.github.com/tothuhien/lsd2/tar.gz/{commit}"
prefix = Path(os.environ["CONDA_PREFIX"])
compiler = os.environ.get("CXX")
if not compiler or not shutil.which(compiler):
    raise RuntimeError("Activate dating.yaml first: its Conda C++ compiler is required")
with tempfile.TemporaryDirectory(prefix=".lsd2-build-", dir=prefix) as temporary:
    work = Path(temporary)
    archive = Path(os.environ["LSD2_SOURCE_ARCHIVE"]) if os.environ.get("LSD2_SOURCE_ARCHIVE") else work / "source.tar.gz"
    if not os.environ.get("LSD2_SOURCE_ARCHIVE"):
        with urlopen(url, timeout=60) as response, archive.open("wb") as handle:
            shutil.copyfileobj(response, handle)
    if hashlib.sha256(archive.read_bytes()).hexdigest() != checksum:
        raise RuntimeError("LSD2 source archive checksum mismatch")
    with tarfile.open(archive) as source:
        source.extractall(work / "source", filter="data")
    source = work / "source" / f"lsd2-{commit}" / "src"
    command = ["make", "-j2", f"CXX={compiler}"]
    subprocess.run(command, cwd=source, check=True)
    binary = source / "lsd2"
    record = {"version": "2.4.4", "commit": commit, "url": url,
              "archive_sha256": checksum, "source_patch_applied": False,
              "command": command,
              "compiler": subprocess.check_output([compiler, "--version"], text=True),
              "executable_sha256": hashlib.sha256(binary.read_bytes()).hexdigest()}
    destination = prefix / "bin" / "lsd2"
    shutil.copy2(binary, destination.with_suffix(".tmp"))
    os.replace(destination.with_suffix(".tmp"), destination)
    manifest = prefix / "share" / "lsd2" / "build.json"
    manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.with_suffix(".tmp").write_text(json.dumps(record, indent=2) + "\n")
    os.replace(manifest.with_suffix(".tmp"), manifest)
PY
