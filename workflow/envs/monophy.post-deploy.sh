#!/usr/bin/env bash
set -euo pipefail

# Self-contained: Snakemake copies and hashes this file with the environment.
"${CONDA_PREFIX:?}/bin/python" - <<'PY'
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from urllib.request import urlopen

version = "1.3.2"
checksum = "e141147f0d11928a8d59bdbeb0a67ff0e80a845d0d1ca30846b0b1773ede2b8c"
prefix = Path(os.environ["CONDA_PREFIX"])
urls = [f"https://cran.r-project.org/src/contrib/MonoPhy_{version}.tar.gz",
        f"https://cran.r-project.org/src/contrib/Archive/MonoPhy/MonoPhy_{version}.tar.gz"]
with tempfile.TemporaryDirectory(prefix=".monophy-install-", dir=prefix) as temporary:
    archive = Path(os.environ.get("MONOPHY_SOURCE_ARCHIVE", Path(temporary) / "MonoPhy.tar.gz"))
    origin = str(archive)
    if not os.environ.get("MONOPHY_SOURCE_ARCHIVE"):
        for url in urls:
            try:
                with urlopen(url, timeout=60) as source, archive.open("wb") as target:
                    shutil.copyfileobj(source, target)
                origin = url
                break
            except OSError:
                if url == urls[-1]:
                    raise
    if hashlib.sha256(archive.read_bytes()).hexdigest() != checksum:
        raise RuntimeError("MonoPhy source archive checksum mismatch")
    subprocess.run([str(prefix / "bin/R"), "CMD", "INSTALL", "--no-multiarch",
                    "--library=" + str(prefix / "lib/R/library"), str(archive)], check=True)
    record = {"package": "MonoPhy", "version": version, "source": origin,
              "archive_sha256": checksum, "source_patch_applied": False}
    manifest = prefix / "share/monophy/install.json"
    manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_text(json.dumps(record, indent=2) + "\n")
PY
