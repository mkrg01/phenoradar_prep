#!/usr/bin/env python3
"""Download NCBI taxonomy and build a snapshot only when it is missing."""
import argparse
import fcntl
import json
import os
import shutil
import subprocess
import sys
import tempfile
from importlib.metadata import version
from pathlib import Path
from urllib.request import urlopen

from common import file_record, write_json
from snapshot_taxonomy import snapshot

TAXDUMP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"


def download_taxdump(destination):
    print(f"Downloading taxonomy from {TAXDUMP_URL}", flush=True)
    with urlopen(TAXDUMP_URL, timeout=60) as response, open(destination, "wb") as handle:
        shutil.copyfileobj(response, handle)


def prepare(destination):
    destination = Path(destination).resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    with open(destination.parent / f".{destination.name}.prepare.lock", "w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if destination.is_file():
            print(f"Keeping existing taxonomy snapshot: {destination}")
            return
        if destination.exists():
            raise ValueError(f"taxonomy destination is not a file: {destination}")
        with tempfile.TemporaryDirectory(prefix=f".{destination.name}.preparing-",
                                         dir=destination.parent) as temporary:
            staging = Path(temporary)
            ready = staging / "snapshot.sqlite"
            taxdump = staging / "taxdump.tar.gz"
            download_taxdump(taxdump)
            generated = staging / "taxa.sqlite"
            # ETE writes taxa.tab, syn.tab, and merged.tab in its working
            # directory. Keep those files and the traversal cache isolated.
            subprocess.run([
                sys.executable, "-c",
                "import sys; from ete4 import NCBITaxa; "
                "ncbi = NCBITaxa(dbfile=sys.argv[1], taxdump_file=sys.argv[2], update=False); "
                "ncbi.db.close()",
                str(generated), str(taxdump),
            ], cwd=staging, check=True)
            snapshot(generated, ready)
            record = json.loads(Path(str(ready) + ".json").read_text())
            record["snapshot"]["path"] = str(destination)
            record["method"] = "ncbi_download"
            record["source"] = TAXDUMP_URL
            record["taxdump"] = {key: value for key, value in file_record(taxdump).items()
                                 if key != "path"}
            record["ete4_version"] = version("ete4")
            # Publish the database last: failures cannot leave a partial SQLite
            # file at the path that subsequent workflow runs will reuse.
            write_json(str(destination) + ".json", record)
            os.replace(ready, destination)
            print(f"Created taxonomy snapshot: {destination}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--destination", required=True)
    prepare(**vars(parser.parse_args()))
