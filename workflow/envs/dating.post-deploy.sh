#!/usr/bin/env bash
set -euo pipefail

# Self-contained: Snakemake copies and hashes this script with the environment.
# Offline installation accepts TREEPL_SOURCE_ARCHIVE, with the same checksum.
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

commit = "f41af04ae7cc830deadbe83a1217ed9feca60c86"
checksum = "9d77514f03fa1583c1b207183e60fe11e2812218cbf83a77e37c643424d423d4"
url = f"https://codeload.github.com/blackrim/treePL/tar.gz/{commit}"
prefix = Path(os.environ["CONDA_PREFIX"])
compiler, cc = os.environ.get("CXX"), os.environ.get("CC")
if not compiler or not cc or not all(shutil.which(c) for c in (compiler, cc)):
    raise RuntimeError("Activate dating.yaml first: its Conda C/C++ compilers are required")
dependencies = {
    "nlopt-2.4.2.tar.gz": "8099633de9d71cbc06cd435da993eb424bbcdbded8f803cdaa9fb8c6e09c8e89",
    "ADOL-C-2.6.3.tgz": "6ed74580695a0d2c960581e5430ebfcd380eb5da9337daf488bf2e89039e9c21",
}
commands = []
def run(argv, cwd):
    commands.append({"argv": argv, "cwd": str(cwd)})
    subprocess.run(argv, cwd=cwd, check=True)

def extract(archive, target, expected):
    if hashlib.sha256(archive.read_bytes()).hexdigest() != expected:
        raise RuntimeError(f"treePL source archive checksum mismatch: {archive.name}")
    with tarfile.open(archive) as source:
        source.extractall(target, filter="data")

with tempfile.TemporaryDirectory(prefix=".treepl-build-", dir=prefix) as temporary:
    work = Path(temporary)
    archive = Path(os.environ["TREEPL_SOURCE_ARCHIVE"]) if os.environ.get("TREEPL_SOURCE_ARCHIVE") else work / "source.tar.gz"
    if not os.environ.get("TREEPL_SOURCE_ARCHIVE"):
        with urlopen(url, timeout=60) as response, archive.open("wb") as handle:
            shutil.copyfileobj(response, handle)
    extract(archive, work / "source", checksum)
    source = work / "source" / f"treePL-{commit}"
    for name, expected in dependencies.items():
        extract(source / "deps" / name, work / "dependencies", expected)
    install = work / "install"
    shared = [f"--prefix={install}", f"--libdir={install / 'lib'}", "--enable-static", "--disable-shared",
              f"CC={cc}", f"CXX={compiler}", "CFLAGS=-O2", "CXXFLAGS=-O2 -std=gnu++11"]
    nlopt = work / "dependencies/nlopt-2.4.2"
    run(["./configure", *shared, "--without-guile", "--without-python", "--without-octave", "--without-matlab"], nlopt)
    run(["make", "-j2"], nlopt)
    run(["make", "install"], nlopt)
    adolc = work / "dependencies/ADOL-C-2.6.3"
    run(["./configure", *shared, "--without-boost", "--without-colpack", "--with-openmp-flag=-fopenmp"], adolc)
    run(["make", "-j2"], adolc)
    run(["make", "install"], adolc)
    src = source / "src"
    run([cc, "-O3", "-fopenmp", "-c", "tnc.c", "-o", "tnc.o"], src)
    # Explicit compiler flags and library paths; no source patches are applied.
    files = ["main.cpp", "myradops.cpp", "node.cpp", "optim_options.cpp", "optimize_nlopt.cpp",
             "optimize_tnc.cpp", "pl_calc_parallel.cpp", "siman_calc_par.cpp", "tree.cpp",
             "tree_reader.cpp", "tree_utils.cpp", "utils.cpp", "tnc.o"]
    binary = work / "treePL"
    run([compiler, "-std=gnu++11", "-O3", "-fopenmp", f"-I{install / 'include'}", *files,
         str(install / "lib/libadolc.a"), str(install / "lib/libnlopt.a"), "-lm", "-o", str(binary)], src)
    record = {"commit": commit, "url": url, "archive_sha256": checksum,
              "source_patch_applied": False, "dependency_archives": dependencies,
              "commands": commands, "compiler": subprocess.check_output([compiler, "--version"], text=True),
              "executable_sha256": hashlib.sha256(binary.read_bytes()).hexdigest()}
    destination = prefix / "bin/treePL"
    shutil.copy2(binary, destination.with_suffix(".tmp"))
    os.replace(destination.with_suffix(".tmp"), destination)
    manifest = prefix / "share/treepl/build.json"
    manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.with_suffix(".tmp").write_text(json.dumps(record, indent=2) + "\n")
    os.replace(manifest.with_suffix(".tmp"), manifest)
PY
