import gzip
import json
import os
import shutil
import sqlite3
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))

from common import file_record, write_json, write_tsv


@pytest.fixture
def workflow_project(tmp_path):
    """Run the real Snakefile with fixed storage paths in an isolated project."""
    (tmp_path / "config").mkdir()
    shutil.copyfile(ROOT / "config/config.yaml", tmp_path / "config/config.yaml")
    return tmp_path


@pytest.fixture
def command_environment(tmp_path):
    """Expose real tools/test doubles under the workflow's fixed command names."""
    def build(commands):
        test_bin = tmp_path / "bin"
        test_bin.mkdir()
        for name, executable in commands.items():
            (test_bin / name).symlink_to(Path(executable).resolve())
        return {**os.environ, "PATH": str(test_bin) + os.pathsep + os.environ.get("PATH", ""),
                "XDG_CACHE_HOME": str(tmp_path / "cache")}
    return build


@pytest.fixture
def tiny_inputs(tmp_path):
    root = tmp_path / "inputs"
    root.mkdir()
    cds = root / "assembly" / "longest_cds"
    cds.mkdir(parents=True)
    quant = root / "assembly" / "quant"
    names = [("Alpha plant", "Alpha_plant", 42, 8), ("Beta sp-X", "Beta_sp-X", 43, 6),
             ("Gamma plant", "Gamma_plant", 44, 4)]
    meta, bus = [], []
    for name, species, taxid, complete in names:
        with gzip.open(cds / f"{species}_longestCDS.fa.gz", "wt") as handle:
            for i in [1, 2, 3]:
                handle.write(f">{species}_g{i}\nATGAAATAA\n")
        runs = ["A1", "A2"] if taxid == 42 else ["B1"] if taxid == 43 else ["G1"]
        for run in runs:
            meta.append({"scientific_name": name, "run": run, "taxid": taxid})
            path = quant / species / run / f"{run}_abundance.tsv"
            tpm = [20, 30, 50] if run != "A2" else [80, 10, 10]
            write_tsv(path, ["target_id", "tpm"],
                      [{"target_id": f"{species}_g{i}", "tpm": value} for i, value in zip([1, 2, 3], tpm)])
        bus.append({"Species": name, "busco_cds_single": complete - 1, "busco_cds_duplicated": 1,
                    "busco_cds_fragmented": 0, "busco_cds_missing": 10 - complete, "busco_cds_total": 10})
    write_tsv(root / "metadata.tsv", list(meta[0]), meta)
    write_tsv(root / "busco.tsv", list(bus[0]), bus)
    taxonomy = root / "taxa.sqlite"
    with sqlite3.connect(taxonomy) as db:
        db.executescript("""
            CREATE TABLE species (taxid INTEGER PRIMARY KEY, parent INTEGER, spname TEXT, common TEXT, rank TEXT, track TEXT);
            CREATE TABLE merged (taxid_old INTEGER, taxid_new INTEGER);
            CREATE TABLE stats (version INTEGER);
            INSERT INTO stats VALUES (2);
            INSERT INTO species VALUES (1, 1, 'root', '', 'no rank', '1');
            INSERT INTO species VALUES (2, 1, 'Plants', '', 'kingdom', '2,1');
        """)
        db.executemany("INSERT INTO species VALUES (?,2,?,'','species',?)",
                       [(tid, name, f"{tid},2,1") for name, _, tid, _ in names])
    return {"metadata": str(root / "metadata.tsv"), "busco": str(root / "busco.tsv"),
            "cds_dir": str(cds), "quant_dir": str(quant), "taxonomy_db": str(taxonomy)}


@pytest.fixture
def fake_odb(tmp_path):
    """A protocol double, confined to tests; never used in production configs."""
    command = tmp_path / "fake_odb"
    command.write_text('''#!/usr/bin/env python3
import json, os, sys, time
from pathlib import Path
root = Path(os.environ['ODBMAPPER_WORK']) / 'v12'
project = root / 'pipeline'
action = sys.argv[1]
if action == 'SETUP':
    for name in ['pipeline', 'etc', 'data/tarfiles']:
        (root / name).mkdir(parents=True, exist_ok=True)
    conf = project / 'orthologer_conf.sh'
    if not conf.exists():
        keys = ['SCHEDULER_LABEL','OP_NJOBMAX_BATCH','OP_NJOBMAX_LOCAL','OP_SAVE_JOBLOG','SKIP_REMAKE_CHECK','STEP_SLEEP','TMP_DIR_BASE']
        steps = ['PREPROC','MASKER','SELECT','STATS','FORMATDB','ALIGNMENT','MAKEBRH','MAKEINPAR','MAKEINPARSEL']
        conf.write_text(''.join(k + '=0\\n' for k in keys) + ''.join('OP_STEP_NPARALLEL[' + s + ']=1\\n' for s in steps))
elif action == 'CONFIG':
    print(str(project) if len(sys.argv) > 2 and sys.argv[2] == 'project' else 'v12')
elif action == 'MAP':
    label, manifest = sys.argv[2:4]
    if os.environ.get('FAKE_ODB_LOG'):
        with open(os.environ['FAKE_ODB_LOG'], 'a') as log:
            log.write(str(root) + '\\n')
    expected = int(os.environ.get('FAKE_ODB_BARRIER_COUNT', '0'))
    if expected:
        deadline = time.monotonic() + 30
        while len(Path(os.environ['FAKE_ODB_LOG']).read_text().splitlines()) < expected:
            if time.monotonic() >= deadline:
                raise SystemExit('mapping jobs did not overlap')
            time.sleep(0.05)
    flag = os.environ.get('FAKE_ODB_FAIL_ONCE')
    if flag and Path(flag).exists():
        Path(flag).unlink()
        sys.exit(23)
    result = project / 'Results'
    result.mkdir(exist_ok=True)
    lines = ['#query\\tODB_OG\\n']
    for name in Path(manifest).read_text().splitlines():
        for line in Path(name).read_text().splitlines():
            if line.startswith('>') and line.endswith(('_g1','_g2')):
                gene = line[1:].split()[0]
                lines.extend([gene + '\\tOG' + gene[-1] + '\\n'] * 2)
    (result / (label + '.og.annotations')).write_text(''.join(lines))
    (result / (label + '.og.hits')).write_text('#cluster_id\\tgene_id\\n')
    (result / (label + '.summary.txt')).write_text('fake ODB summary\\n')
elif action in ['REPORT', 'DBINFO']:
    print('fake ODB report')
else:
    raise SystemExit('unexpected fake ODB action: ' + action)
''')
    command.chmod(0o755)
    return command


@pytest.fixture
def frozen_reference(tmp_path):
    root = tmp_path / "reference"
    data = root / "odbmapper" / "v12" / "data"
    (data / "tarfiles").mkdir(parents=True)
    (data / "tiny.db").write_text("synthetic reference")
    inventory = [{"relative_path": "tiny.db", **file_record(data / "tiny.db")}]
    write_json(root / "files.json", inventory)
    write_json(root / "reference.json", {"version": "v12", "node": 3193, "data_dir": str(data),
                                         "inventory": file_record(root / "files.json")})
    for name in ["dbinfo.txt", "config.txt"]:
        (root / name).write_text("synthetic reference\n")
    return root
