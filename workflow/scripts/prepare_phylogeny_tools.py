#!/usr/bin/env python3
"""Build pinned ASTRAL-IV int128, patched treePL and trimAl locally."""
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
    "treepl": {"repository": "blackrim/treePL", "commit": "f41af04ae7cc830deadbe83a1217ed9feca60c86",
               "sha256": "9d77514f03fa1583c1b207183e60fe11e2812218cbf83a77e37c643424d423d4", "binary": "treePL"},
    "trimal": {"repository": "inab/trimal", "commit": "d637091abe33595775f40480970d1a18d87a7bcb",
               "sha256": "1a5e0b24cfd1a5912725cbcabc86a72e768be2541a1ad2a412885ea5d29b55ba", "binary": "trimal"},
}
TREEPL_PATCH = Path(__file__).resolve().parents[1] / "patches/treepl.patch"


def build_treepl(source, tmp, compiler):
    """Use dependencies bundled in the verified archive; link them statically."""
    commands = []
    def run(argv, cwd):
        commands.append({"argv": argv, "cwd": str(cwd)})
        subprocess.run(argv, cwd=cwd, check=True)
    run(["patch", "-p1", "--batch", "-i", str(TREEPL_PATCH)], source)
    prefix = tmp / "dependencies"
    for archive, name, options in [
        ("nlopt-2.4.2.tar.gz", "nlopt-2.4.2", ["--without-guile", "--without-python", "--without-octave", "--without-matlab"]),
        ("ADOL-C-2.6.3.tgz", "ADOL-C-2.6.3", ["--without-boost", "--without-colpack", "--with-openmp-flag=-fopenmp"]),
    ]:
        with tarfile.open(source / "deps" / archive) as handle:
            handle.extractall(tmp / "dependency_sources", filter="data")
        cwd = tmp / "dependency_sources" / name
        run(["./configure", "--prefix=" + str(prefix), "--libdir=" + str(prefix / "lib"), "--enable-static", "--disable-shared",
             *options, "CXX=" + compiler, "CXXFLAGS=-O2 -std=gnu++11"], cwd)
        run(["make", "-j2"], cwd)
        run(["make", "install"], cwd)
    cwd, built = source / "src", tmp / "treePL"
    run(["gcc", "-O3", "-fopenmp", "-c", "tnc.c", "-o", "tnc.o"], cwd)
    run([compiler, "-std=gnu++11", "-O3", "-fopenmp", "-I" + str(prefix / "include"),
         *sorted(p.name for p in cwd.glob("*.cpp") if p.name != "radops.cpp"), "tnc.o",
         str(prefix / "lib/libadolc.a"), str(prefix / "lib/libnlopt.a"), "-lm", "-o", str(built)], cwd)
    return built, commands


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
            patch_matches = name != "treepl" or old.get("patch", {}).get("sha256") == sha256(TREEPL_PATCH)
            if old.get("source") == spec and old["executable"]["sha256"] == sha256(target) and patch_matches:
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
            if name == "aster":
                built = tmp / spec["binary"]
                command = [compiler, "-std=gnu++17", "-O3", "-pthread", "-D", "LARGE_DATA",
                           str(source / "src/astral.cpp"), "-o", str(built)]
                subprocess.run(command, check=True)
            elif name == "treepl":
                built, command = build_treepl(source, tmp, compiler)
            else:
                built = source / "source/trimal"
                command = ["make", "-C", str(source / "source"), "trimal", "-j2", f"CC={compiler}"]
                subprocess.run(command, check=True)
            if not built.is_file():
                raise ValueError(f"build did not produce {built}")
            with atomic_writer(target, "wb") as output, open(built, "rb") as handle:
                shutil.copyfileobj(handle, output)
            os.chmod(target, 0o755)
            write_json(manifest, {"created_at": now(), "source": spec, "url": url,
                       "command": command, "compiler": subprocess.check_output([compiler, "--version"], text=True),
                       "executable": file_record(target),
                       "patch": file_record(TREEPL_PATCH) if name == "treepl" else None})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--destination", default="resources/phylogeny_tools")
    parser.add_argument("--archives", help="offline directory with verified aster.tar.gz, treepl.tar.gz and trimal.tar.gz")
    prepare(**vars(parser.parse_args()))
