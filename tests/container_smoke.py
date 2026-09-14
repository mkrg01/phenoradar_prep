#!/usr/bin/env python3
"""Exercise the built image without reference downloads or test-only packages."""
import argparse
import importlib
import json
import math
import os
from pathlib import Path
import shutil
import sqlite3
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow/scripts"))
MANIFEST = Path("/opt/phenoradar/container.json")


def dating(work):
    from common import write_json, write_tsv
    from date_phylogeny import date
    from infer_phylogeny import read_tree
    tree = work / "tree.nwk"
    tree.write_text("((A:0.03,B:0.07):0.04,(C:0.03,D:0.07):0.04);\n")
    provenance = work / "tree.json"
    write_json(provenance, {"branch_length_unit": "substitutions_per_site", "outgroup": "A",
                           "total_gene_sites": 32000})
    bounds = work / "bounds.tsv"
    write_tsv(bounds, ["taxa", "min_age_ma", "max_age_ma", "source"], [
        {"taxa": "A,D", "min_age_ma": 100, "max_age_ma": 100, "source": "smoke root"},
        {"taxa": "C,D", "min_age_ma": 50, "max_age_ma": 60, "source": "smoke internal"}])
    for variance in (0, 1, 2):
        output = work / f"dating-{variance}"
        date(tree, provenance, bounds, output, "lsd2", {"variance": variance})
        result = read_tree(output / "species_tree.dated.nwk", {"A", "B", "C", "D"})
        assert all(math.isclose(result.get_distance(result, leaf), 100, abs_tol=1e-8)
                   for leaf in result.leaves())
        age = 100 - result.get_distance(result, result.common_ancestor(["C", "D"]))
        assert 50 - 1e-8 <= age <= 60 + 1e-8
        report = json.loads((output / "provenance.json").read_text())
        assert report["topology_preserved"] and report["build"]["source_patch_applied"] is False


def monophy(work):
    from common import read_tsv, write_json, write_tsv
    from taxonomy_audit import audit
    names = ["A_one", "A_query", "A_three", "A_two", "B_five", "B_four",
             "B_one", "B_three", "B_two", "Outgroup"]
    taxonomy = work / "taxa.sqlite"
    with sqlite3.connect(taxonomy) as db:
        db.executescript("""
            CREATE TABLE species(taxid INTEGER PRIMARY KEY, spname TEXT, rank TEXT, track TEXT);
            CREATE TABLE merged(taxid_old INTEGER, taxid_new INTEGER);
            INSERT INTO species VALUES (1,'root','no rank','1');
            INSERT INTO species VALUES (10,'FamilyA','family','10,1');
            INSERT INTO species VALUES (20,'FamilyB','family','20,1');
            INSERT INTO species VALUES (30,'FamilyO','family','30,1');
        """)
        for index, name in enumerate(names):
            family = 10 if name.startswith("A_") else 20 if name.startswith("B_") else 30
            db.execute("INSERT INTO species VALUES (?,?,?,?)",
                       (100 + index, name, "species", f"{100 + index},{family},1"))
    samples = work / "samples.tsv"
    write_tsv(samples, ["species", "scientific_name", "taxid", "run"],
              [dict(species=name, scientific_name=name.replace("_", " "), taxid=100+i, run=f"SRR{i}")
               for i, name in enumerate(names)])
    tree = work / "tree.nwk"
    tree.write_text("(Outgroup:1,(((A_one:1,A_two:1):1,A_three:1):1,"
                    "((B_one:1,A_query:1):1,((B_two:1,B_three:1):1,(B_four:1,B_five:1):1):1):1):1);\n")
    qc = work / "tree.json"
    write_json(qc, {"outgroup": "Outgroup", "species": len(names)})
    output = work / "audit"
    result = audit(tree, qc, samples, taxonomy, output, {"ranks": ["family"]})
    assert result["method"] == "MonoPhy" and result["engine"]["version"] == "1.3.2"
    assert {(r["species"], r["role"], r["focal_taxon"]) for r in read_tsv(output / "candidates.tsv")} == {
        ("A_query", "outlier", "FamilyA"), ("A_query", "intruder", "FamilyB")}
    assert (output / "ranks/family/tree.pdf").stat().st_size > 1000
    assert (output / "ranks/family/tree.svg").stat().st_size > 1000


def phylogeny(work):
    from common import sha256
    from prepare_phylogeny_tools import prepare
    bundle = work / "tools"
    prepare(bundle)
    binary = bundle / "bin/astral4_int128"
    record = json.loads((bundle / "aster.json").read_text())
    assert record["executable"]["sha256"] == sha256(binary)
    assert "LARGE_DATA" in record["command"]
    trees = work / "genes.nwk"
    trees.write_text("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);\n" * 4)
    output = work / "species.nwk"
    subprocess.run([str(binary), "-i", str(trees), "-o", str(output), "-t", "1",
                    "--root", "A", "--seed", "1"], check=True)
    from infer_phylogeny import read_tree
    read_tree(output, {"A", "B", "C", "D"})


def check_environment(name):
    manifest = json.loads(MANIFEST.read_text())
    assert Path(sys.prefix) == Path(manifest["environments"][name]), sys.prefix
    imports = {"analysis": ["pandas", "ete4", "matplotlib"], "dating": ["Bio", "ete4"],
               "monophy": ["ete4"], "phylogeny": ["Bio", "ete4", "numpy", "cdskit"],
               "timetree": ["Bio", "nwkit"]}
    commands = {"alignment": ["famsa"], "kofam": ["exec_annotation", "hmmsearch"],
                "odb": ["ODB-mapper", "findmnt"], "seqkit": ["seqkit"],
                "phylogeny": ["famsa", "trimal", "VeryFastTree"],
                "dating": ["lsd2"], "monophy": ["Rscript"]}
    for module in imports.get(name, []):
        importlib.import_module(module)
    for command in commands.get(name, []):
        assert shutil.which(command), f"Missing {command} in {name}"
    with tempfile.TemporaryDirectory(prefix=f"phenoradar-{name}-") as temporary:
        work = Path(temporary)
        if name in {"dating", "monophy", "phylogeny"}:
            globals()[name](work)
        elif name == "seqkit":
            result = subprocess.run(["seqkit", "translate"], input=">test\nATGAAATAA\n",
                                    text=True, capture_output=True, check=True)
            assert "MK*" in result.stdout
        elif name == "alignment":
            source, target = work / "input.faa", work / "aligned.faa"
            source.write_text(">a\nMKWVTFISLLFLFSSAYS\n>b\nMKWVTFISLFLFSSAYS\n")
            subprocess.run(["famsa", "-t", "1", str(source), str(target)], check=True)
            assert target.read_text().count(">") == 2
    print(f"PASS {name}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--environment", help="check one environment after activation")
    args = parser.parse_args()
    if args.environment:
        check_environment(args.environment)
        return
    manifest = json.loads(MANIFEST.read_text())
    for name, prefix in sorted(manifest["environments"].items()):
        subprocess.run(["conda", "run", "--no-capture-output", "--prefix", prefix,
                        "python", str(Path(__file__).resolve()), "--environment", name],
                       check=True, env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"})
    print("PASS all image environments (no reference downloads)")


if __name__ == "__main__":
    main()
