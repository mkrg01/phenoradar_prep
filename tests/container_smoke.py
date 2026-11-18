#!/usr/bin/env python3
"""Exercise the built image without reference downloads or test-only packages."""
import argparse
import importlib
import json
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


def rooting(work):
    from ete4.ncbi_taxonomy.ncbiquery import DB_VERSION
    from common import read_tsv, write_json, write_tsv
    from contrast_pairs import skim
    from infer_phylogeny import read_tree
    from nwkit_runtime import nwkit_backend
    from phylogeny_root import ncbi_tree, prepare_root, root_species_tree

    work.mkdir()
    names = ["O1", "O2", "A", "B", "C"]
    scores = dict(zip(names, [95, 90, 98, 100, 99]))
    taxonomy = work / "taxa.sqlite"
    with sqlite3.connect(taxonomy) as db:
        db.executescript("""
            CREATE TABLE species(taxid INTEGER PRIMARY KEY, parent INTEGER, spname TEXT, common TEXT, rank TEXT, track TEXT);
            CREATE TABLE merged(taxid_old INTEGER, taxid_new INTEGER);
            CREATE TABLE stats(version INTEGER);
            INSERT INTO species VALUES (1,1,'root','','no rank','1');
            INSERT INTO species VALUES (2,1,'Plants','','kingdom','2,1');
            INSERT INTO species VALUES (3,2,'Outgroup','','family','3,2,1');
            INSERT INTO species VALUES (4,2,'Ingroup','','family','4,2,1');
        """)
        db.execute("INSERT INTO stats VALUES (?)", (DB_VERSION,))
        for taxid, name in enumerate(names, 100):
            parent = 3 if name.startswith("O") else 4
            db.execute("INSERT INTO species VALUES (?,?,?,'','species',?)",
                       (taxid, parent, f"Plant {name}", f"{taxid},{parent},2,1"))
    samples, metadata = work / "samples.tsv", work / "metadata.tsv"
    write_tsv(samples, ["species", "species_id", "scientific_name", "taxid", "cds"], [
        dict(species=n, species_id=f"Plant_{n}", scientific_name=f"Plant {n}", taxid=i, cds="unused")
        for i, n in enumerate(names, 100)])
    write_tsv(metadata, ["species", "busco_percent"],
              [dict(species=n, busco_percent=scores[n]) for n in names])
    guide, anchor, anchor_qc = work / "guide.nwk", work / "anchor.txt", work / "anchor.json"
    ncbi_tree(samples, taxonomy, guide, work / "taxids.tsv")
    prepare_root(samples, metadata, anchor, anchor_qc, tree=guide, for_anchor=True)
    assert anchor.read_text() == "O1\n"
    raw, raw_qc = work / "astral.nwk", work / "astral.json"
    raw.write_text("(O1:2,(O2:3,(A:4,(B:5,C:6)0.8:7)0.9:8):2);\n")
    write_json(raw_qc, {"astral_root_anchor": "O1", "outgroup": ["O1"], "branch_length_method": "CASTLES-II"})
    output, qc, outgroup, rooting_qc = [work / n for n in ["tree.nwk", "tree.json", "outgroup.txt", "rooting.json"]]
    root_species_tree(samples, raw, raw_qc, anchor_qc, output, qc, outgroup, rooting_qc, taxonomy_db=taxonomy)
    assert outgroup.read_text().splitlines() == ["O1", "O2"]
    assert json.loads(rooting_qc.read_text())["nwkit_version"] == nwkit_backend().__version__
    tree = read_tree(output, names)
    rows = [dict(leaf_name=n, trait="0" if n.startswith("O") else "1", busco_percent=scores[n]) for n in names]
    skim(tree, rows, work / "skim", 12345)
    assert {r["leaf_name"] for r in read_tsv(work / "skim.sampled.tsv")} == {"O1", "B"}


def dating(work):
    from unittest.mock import patch

    from angiocal_calibrations import prepare
    from common import read_tsv, write_json, write_tsv
    from date_phylogeny import date
    from infer_phylogeny import read_tree
    from timetree_calibrations import prepare as prepare_timetree

    rooting(work / "rooting")

    names = ("Nuphar_lutea", "Nymphaea_colorata", "Brasenia_schreberi", "Cabomba_caroliniana")
    tree = work / "tree.nwk"
    tree.write_text("((Nuphar_lutea:0.06,Nymphaea_colorata:0.06):0.06,"
                    "(Brasenia_schreberi:0.06,Cabomba_caroliniana:0.06):0.06);\n")
    qc, samples, taxonomy = work / "tree.json", work / "samples.tsv", work / "taxa.sqlite"
    write_json(qc, {"branch_length_unit": "substitutions_per_site", "outgroup": list(names[:2]),
                    "total_gene_sites": 32000})
    with sqlite3.connect(taxonomy) as db:
        db.executescript("""
            CREATE TABLE species(taxid INTEGER PRIMARY KEY, spname TEXT, rank TEXT, track TEXT);
            CREATE TABLE synonym(taxid INTEGER, spname TEXT);
            CREATE TABLE merged(taxid_old INTEGER, taxid_new INTEGER);
            INSERT INTO species VALUES (3398,'Magnoliopsida','class','3398');
            INSERT INTO species VALUES (10,'Nymphaeaceae','family','10,3398');
            INSERT INTO species VALUES (20,'Cabombaceae','family','20,3398');
        """)
        for index, name in enumerate(names):
            taxid, family = 42 + index, 10 if index < 2 else 20
            db.execute("INSERT INTO species VALUES (?,?,?,?)",
                       (taxid, name.replace("_", " "), "species", f"{taxid},{family},3398"))
    write_tsv(samples, ["species", "taxid"],
              [{"species": name, "taxid": 42 + index} for index, name in enumerate(names)])
    mapping = work / "mapping.tsv"
    write_tsv(mapping, ["fossil_id", "left_species", "right_species"],
              [{"fossil_id": 166, "left_species": names[0], "right_species": names[2]}])
    calibration = work / "angiocal"
    prepare(tree, samples, taxonomy, calibration, mapping)
    fossils = {row["fossil_id"]: row for row in read_tsv(calibration / "fossils.tsv")}
    assert fossils["166"]["mapping_method"] == "explicit_map" and fossils["166"]["status"] == "applied"
    assert fossils["165"]["mapping_method"] == "taxonomy_stem"
    output = work / "dating"
    date(tree, qc, calibration / "calibrations.tsv", output)
    result = read_tree(output / "species_tree.dated.nwk", names)
    assert all(113 - 2e-3 <= result.get_distance(result, leaf) <= 154 + 2e-3 for leaf in result.leaves())
    report = json.loads((output / "provenance.json").read_text())
    assert report["topology_preserved"] and report["build"]["source_patch_applied"] is False
    assert report["smoothing"] == 10 and report["root_max_age_ma"] == 154

    responses = {"42+43+44+45": (3398, 185, 170, 195), "42+43": (0, 90, 70, 105),
                 "44+45": (0, 100, 80, 120)}
    responses.update({str(i): (i, 0, 0, 0) for i in range(42, 46)})

    def fetch(url):
        row = responses.get(url.rsplit("/", 1)[-1])
        text = ("taxon_id,precomputed_age,precomputed_ci_low,precomputed_ci_high\r\n"
                + ",".join(map(str, row))) if row else "MRCA node not found"
        return {"error": None, "status_code": 200, "text": text, "elapsed": 0}

    secondary = work / "timetree"
    with patch("nwkit.mcmctree._fetch_timetree_url", side_effect=fetch):
        prepare_timetree(tree, samples, taxonomy, secondary)
        assert len(read_tsv(secondary / "calibrations.tsv")) == 3
        date(tree, qc, secondary / "calibrations.tsv", work / "dating-timetree", calibration_source="timetree")
        report = json.loads((work / "dating-timetree/provenance.json").read_text())
        assert report["calibration_source"] == "TimeTree" and 170 <= report["root_age_ma"] <= 195

        responses["42+43+44+45"] = responses["42+43"]
        prepare_timetree(tree, samples, taxonomy, secondary)
        assert len(read_tsv(secondary / "calibrations.tsv")) == 2


def monophy(work):
    from common import read_tsv, write_json, write_tsv
    from taxonomy_check import check
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
    output = work / "check"
    result = check(tree, qc, samples, taxonomy, output, {"ranks": ["family"]})
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
    imports = {"alignment": ["alignment_selection"],
               "analysis": ["pandas", "ete4", "matplotlib"], "dating": ["nwkit", "ete4"],
               "monophy": ["ete4"], "phylogeny": ["Bio", "ete4", "numpy", "cdskit"]}
    commands = {"alignment": ["famsa", "trimal"], "kofam": ["exec_annotation", "hmmsearch"],
                "odb": ["ODB-mapper"], "seqkit": ["seqkit"],
                "phylogeny": ["famsa", "trimal", "VeryFastTree"],
                "dating": ["treePL"], "monophy": ["Rscript"]}
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
            from alignment_trimming import trim_records
            from busco_phylogeny import fasta_records
            records, columns, _, _ = trim_records(list(fasta_records(target)), 'trimal', 'gappyout', work)
            assert len(records) == 2 and columns
        elif name == "odb":
            from check_odb_runtime import check
            assert check(work_dir=work)["checks"]["fasta_import"] == "passed"
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
