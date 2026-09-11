"""BUSCO reuse, missing-data safeguards, and real phylogeny tool integration."""
import csv
import gzip
import json
import os
from pathlib import Path
import random
import shutil
import sqlite3
import subprocess
import sys

import pytest
import yaml

from busco_phylogeny import busco_table, extract, fasta_records, plan, prepare_cds
from common import read_tsv, write_json, write_tsv
from date_phylogeny import date, calibration_rows
from infer_phylogeny import alignment_qc, astral, merge, read_tree, trim

ROOT = Path(__file__).resolve().parents[1]


def settings(**updates):
    cfg = yaml.safe_load((ROOT / "config/config.yaml").read_text())["phylogeny"]
    cfg.update(translation_table=1)
    cfg.update(updates)
    return cfg


def table(path, rows, lineage="embryophyta_odb12"):
    path.write_text(f"# BUSCO version is: 5.8.2\n# The lineage dataset is: {lineage} (test)\n"
                   "# Busco id\tStatus\tSequence\tScore\tLength\n" + "\n".join(rows) + "\n")


def test_busco_single_copy_and_lineage_validation(tmp_path):
    path = tmp_path / "busco.tsv"
    table(path, ["1at1\tComplete\tg1:0-299\t100\t100", "2at1\tDuplicated\tg2\t100\t100",
                 "2at1\tDuplicated\tg3\t100\t100", "3at1\tFragmented\tg4\t20\t20",
                 "4at1\tMissing", "5at1\tComplete\tg5:0-99\t50\t30",
                 "6at1\tComplete\tg5:101-199\t50\t30"])
    complete, universe = busco_table(path, "embryophyta_odb12")
    assert set(complete) == {"1at1"}
    assert len(universe) == 6
    with pytest.raises(ValueError, match="lineage/header"):
        busco_table(path, "embryophyta_odb10")
    path.write_text("# The lineage dataset is: embryophyta_odb12\n"
                    "# Busco id\tStatus\tSequence\tGene Start\tGene End\n")
    with pytest.raises(ValueError, match="genome BUSCO"):
        busco_table(path, "embryophyta_odb12")


@pytest.mark.parametrize("limit,expected", [(2, ["4at1", "6at1"]),
                                           (500, ["4at1", "6at1", "1at1", "3at1", "2at1", "5at1"])])
def test_plan_selects_highest_coverage_with_id_ties_independent_of_orders(tmp_path, limit, expected):
    # Thirty species in two orders (20 + 10): 1at1 is absent from the second
    # order. Its tied competitor 3at1 spans both orders and has longer hits.
    # Neither order representation nor match length may break the coverage tie.
    # A short BUSCO match (6at1) must remain eligible; extracted-sequence QC is separate.
    # Low overall coverage (5at1) is eligible too, but fewer than four species is not.
    species = [f"plant_{i:02d}" for i in range(30)]
    presence = {"1at1": set(species[:18]), "2at1": set(species[:17]),
                "3at1": set(species[:12] + species[20:26]), "4at1": set(species),
                "5at1": set(species[:9]), "6at1": set(species), "7at1": set(species[:3])}
    lengths = {"1at1": 100, "2at1": 2000, "3at1": 900, "4at1": 120, "5at1": 300, "6at1": 99, "7at1": 150}
    cds = tmp_path / "cds.fa"
    cds.write_text(">placeholder\n" + "ATG" * 100 + "\n")
    tables = tmp_path / "busco"
    tables.mkdir()
    for name in species:
        table(tables / f"{name}.busco.full.tsv",
              [f"{m}\tComplete\t{name}_g{m}\t100\t{lengths[m]}" if name in presence[m]
               else f"{m}\tMissing" for m in reversed(list(presence))])
    samples = tmp_path / "samples.tsv"
    rows = [{"species": name, "cds": str(cds), "order": "X" if i < 20 else "Y"}
            for i, name in enumerate(species)]
    cfg = settings(busco_full_dir=str(tables), outgroup=species[0], max_markers=limit)
    del cfg["min_protein_length"]  # Planning must not depend on sequence-length QC.
    selected = []
    for run in range(2):
        if run:
            # Changes to order labels, row order and replicate count must not
            # affect the species-level denominator or the selected loci.
            rows = [dict(r, order="") for r in reversed(rows)] + [dict(rows[0], order="Z")]
        write_tsv(samples, ["species", "cds", "order"], rows)
        out = tmp_path / f"plan_{run}"
        plan(samples, out, cfg)
        chosen = read_tsv(out / "markers.tsv")
        assert [r["marker"] for r in chosen] == expected
        assert [int(r["selection_rank"]) for r in chosen] == list(range(1, len(expected) + 1))
        stats = {r["marker"]: r for r in read_tsv(out / "marker_stats.tsv")}
        assert float(stats["1at1"]["occupancy"]) == pytest.approx(0.6)
        assert stats["1at1"]["selection_rank"] == "3"
        assert stats["3at1"]["selection_rank"] == "4"
        assert stats["5at1"]["eligible"] == "True"
        assert stats["5at1"]["selection_rank"] == "6"
        assert float(stats["5at1"]["occupancy"]) == pytest.approx(0.3)
        assert stats["7at1"]["eligible"] == "False"
        assert stats["7at1"]["selection_rank"] == ""
        assert stats["6at1"]["eligible"] == "True"
        assert stats["6at1"]["selection_rank"] == "2"
        assert float(stats["6at1"]["mean_busco_length"]) == 99
        assert json.loads((out / "provenance.json").read_text())["eligible_markers"] == 6
        selected.append(chosen)
    assert selected[0] == selected[1]


def test_extract_original_cds_audits_padding_masking_and_missing_ids(tmp_path):
    full, cds, markers = tmp_path / "busco.tsv", tmp_path / "cds.fa.gz", tmp_path / "markers.tsv"
    table(full, ["1at1\tComplete\tg1:3-11\t100\t3", "2at1\tComplete\tg2:0-11\t100\t3"])
    # Coordinates are MetaEuk hit metadata, not instructions to crop the original CDS.
    with gzip.open(cds, "wt") as handle:
        handle.write(">g1\nATGAAAGGGTAAA\n>g2\nATGGCNTARAAA\n")
    write_tsv(markers, ["marker"], [{"marker": "1at1"}, {"marker": "2at1"}])
    output, qc = tmp_path / "out.fa", tmp_path / "out.json"
    original = cds.read_bytes()
    extract("Plant_alpha", full, cds, markers, output, qc,
            settings(min_protein_length=2, max_unknown_fraction=0.5))
    assert list(fasta_records(output)) == [("1at1", "XERVX"), ("2at1", "MAXK")]
    assert cds.read_bytes() == original
    report = json.loads(qc.read_text())
    assert report["rejected"] == {}
    import cdskit
    assert report["cdskit"]["version"] == cdskit.__version__
    first, second = [r["preparation"] for r in report["records"]]
    assert first["head_padding_nt"] == first["tail_padding_nt"] == 1
    assert first["reading_frame_changed"] is True
    assert second["internal_stops_after_padding"] == 1
    assert second["mask_changed_codon_indices_0based"] == [2]
    extract("Plant_alpha", full, cds, markers, output, qc, settings(min_protein_length=2))
    assert json.loads(qc.read_text())["rejected"] == {"ambiguous_protein": 2}
    assert not output.read_text()
    table(full, ["1at1\tComplete\tabsent:0-11\t100\t3"])
    with pytest.raises(ValueError, match="absent from original"):
        extract("Plant_alpha", full, cds, markers, output, qc, settings(min_protein_length=2))


@pytest.mark.parametrize("sequence,code,expected,head,tail", [
    ("ATGAAAGGGTAA", 1, "MKGX", 0, 0),
    ("TGAAAGGGTAA", 1, "XKGX", 1, 0),
    ("GAAAGGGTAA", 1, "ERVX", 0, 2),
    ("ATGTAAAAAGGG", 1, "XVKRX", 1, 2),
    ("ATGGCNTARAAA", 1, "MAXK", 0, 0),
    ("ATGTGAAAATAA", 4, "MWKX", 0, 0),
    ("ATG---AAA", 1, "M-K", 0, 0),
    ("ATGA--AAA", 1, "MXK", 0, 0),
])
def test_cdskit_preparation_matches_real_cli(tmp_path, sequence, code, expected, head, tail):
    protein, audit = prepare_cds(sequence, "gene", code)
    assert protein == expected
    assert (audit["head_padding_nt"], audit["tail_padding_nt"]) == (head, tail)
    current = tmp_path / "input.fa"
    current.write_text(f">gene\n{sequence}\n")
    for action in ["pad", "mask", "translate"]:
        output = tmp_path / f"{action}.fa"
        subprocess.run([sys.executable, "-m", "cdskit.cli", action, "--seq_file", str(current),
                        "--out_file", str(output), "--codon_table", str(code)], check=True,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        current = output
    assert list(fasta_records(current)) == [("gene", protein)]


def test_protein_input_requires_known_residues_and_does_not_reframe(tmp_path):
    full, proteins, markers = tmp_path / "busco.tsv", tmp_path / "input.fa", tmp_path / "markers.tsv"
    table(full, ["1at1\tComplete\tg1\t100\t3", "2at1\tComplete\tg2\t100\t3",
                 "3at1\tComplete\tg3\t100\t3"])
    proteins.write_text(">g1\nMX*KG\n>g2\nMAXX\n>g3\nMKG*\n")
    write_tsv(markers, ["marker"], [{"marker": f"{i}at1"} for i in range(1,4)])
    out, qc = tmp_path / "out.fa", tmp_path / "qc"
    extract("plant", full, proteins, markers, out, qc,
            settings(sequence_mode="protein", min_protein_length=3, max_unknown_fraction=0.5))
    assert list(fasta_records(out)) == [("3at1", "MKG")]
    assert json.loads(qc.read_text())["rejected"] == {"internal_stop": 1, "short_protein": 1}


def trimal_binary():
    command = os.environ.get("TRIMAL_BIN") or shutil.which("trimal")
    if not command:
        pytest.skip("set TRIMAL_BIN for real trimming tests")
    return command


@pytest.mark.parametrize("mode", ["gappyout", "automated1"])
def test_real_trimal_column_map_and_missing_taxa(tmp_path, mode):
    raw, raw_qc = tmp_path / "raw.fa", tmp_path / "raw.json"
    raw.write_text(">A\nAAX-GGWKY\n>B\nAAX-GGWKY\n>C\nCCX-GGXKY\n>D\nCCX-GGWKY\n>E\nXXXXXXXKY\n")
    write_json(raw_qc, {"status": "retained", "sites": 9, "taxa": 5})
    out, qc, columns = tmp_path / "out.fa", tmp_path / "qc.json", tmp_path / "columns.tsv"
    trim(raw, raw_qc, out, qc, columns, trimal_binary(), mode,
         settings(min_protein_length=3))
    assert list(fasta_records(out)) == [("A", "AAGGKY"), ("B", "AAGGKY"), ("C", "CCGGKY"), ("D", "CCGGKY")]
    assert [int(r["famsa_column_1based"]) for r in read_tsv(columns)] == [1, 2, 5, 6, 8, 9]
    report = json.loads(qc.read_text())
    assert report["removed_species"] == ["E"]
    assert report["informative_sites"] == 2
    assert report["mode"] == mode
    assert report["raw_sites"] == 9


def test_real_trimal_retains_original_unknown_residues(tmp_path):
    raw, raw_qc = tmp_path / "raw.fa", tmp_path / "raw.json"
    # Uniformly 25% missing: gappyout retains these columns, including X.
    raw.write_text(">A\nXACC\n>B\nAXCC\n>C\nCCXA\n>D\nCCAX\n")
    write_json(raw_qc, {"status": "retained", "sites": 4, "taxa": 4})
    out, qc, columns = tmp_path / "out.fa", tmp_path / "qc.json", tmp_path / "columns.tsv"
    trim(raw, raw_qc, out, qc, columns, trimal_binary(), "gappyout",
         settings(min_protein_length=3))
    assert out.read_text() == raw.read_text()
    assert len(read_tsv(columns)) == 4


def test_trimal_rejects_inconsistent_column_map(tmp_path):
    raw, raw_qc = tmp_path / "raw.fa", tmp_path / "raw.json"
    raw.write_text("".join(f">{n}\nAAGG\n" for n in "ABCD"))
    write_json(raw_qc, {"status": "retained", "sites": 4, "taxa": 4})
    bad = tmp_path / "bad-trimal"
    bad.write_text("#!/bin/sh\nprintf '#ColumnsMap\\t0, 0, 1\\n'\n")
    bad.chmod(0o755)
    with pytest.raises(ValueError, match="duplicated, unordered or out of range"):
        trim(raw, raw_qc, tmp_path / "out", tmp_path / "qc", tmp_path / "cols", str(bad), "gappyout", settings())


def test_alignment_qc_drops_unusable_taxa_without_fabricating_signal():
    records = [("A", "AA--GG"), ("B", "AA--GG"), ("C", "CC--GG"),
               ("D", "CC--GG"), ("E", "------")]
    retained, report = alignment_qc(records, settings(min_protein_length=2))
    assert retained == [("A", "AAGG"), ("B", "AAGG"), ("C", "CCGG"), ("D", "CCGG")]
    assert report["informative_sites"] == 2
    assert report["removed_species"] == ["E"]
    assert report["status"] == "retained"
    _, report = alignment_qc([(n, "AAAA") for n in "ABCD"], settings(min_protein_length=2))
    assert report["status"] == "no_variable_sites"
    assert report["variable_sites"] == report["informative_sites"] == 0


@pytest.mark.parametrize("sequences,variable,informative", [
    (["A" * 100] * 2 + ["C" * 8 + "A" * 92] * 2, 8, 8),
    (["A" * 100] * 3 + ["C" + "A" * 99], 1, 0),
])
def test_alignment_qc_reports_low_information_without_a_count_cutoff(sequences, variable, informative):
    records = list(zip("ABCD", sequences))
    retained, report = alignment_qc(records, settings())
    assert retained == records
    assert report["status"] == "retained"
    assert report["sites"] == 100
    assert report["variable_sites"] == variable
    assert report["informative_sites"] == informative


def test_calibrations_require_absolute_ages_and_consistent_ancestors(tmp_path):
    from ete4 import Tree
    tree = Tree("((A:1,B:1)N2:2,(C:1,D:1)N3:2)N1;", parser=1)
    path = tmp_path / "calibrations.tsv"
    fields = ["taxa", "min_age_ma", "max_age_ma", "source"]
    write_tsv(path, fields, [])
    with pytest.raises(ValueError, match="at least one"):
        calibration_rows(tree, path)
    write_tsv(path, fields, [{"taxa": "A,C", "min_age_ma": 10, "max_age_ma": 20, "source": "test"},
                             {"taxa": "A,B", "min_age_ma": 30, "max_age_ma": 40, "source": "test"}])
    with pytest.raises(ValueError, match="inconsistent"):
        calibration_rows(tree, path)


def test_astral_rejects_non_int128_for_large_species_sets(tmp_path):
    manifest = tmp_path / "species.tsv"
    write_tsv(manifest, ["species"], [{"species": f"s{i}"} for i in range(5001)])
    with pytest.raises(ValueError, match="requires astral4_int128"):
        astral("unused", "unused", manifest, "unused", "unused", sys.executable, "s0", 1, 1)


def test_missing_species_diagnostics_survive_before_inference(tmp_path):
    manifest, markers = tmp_path / "species.tsv", tmp_path / "markers.tsv"
    write_tsv(manifest, ["species"], [{"species": n} for n in "ABCDE"])
    write_tsv(markers, ["marker"], [{"marker": "1at1"}])
    (tmp_path / "1at1.nwk").write_text("((A:1,B:1):1,(C:1,D:1):1);\n")
    write_json(tmp_path / "1at1.json", {"status": "retained", "sites": 250})
    output, coverage, qc = tmp_path / "genes.nwk", tmp_path / "coverage.tsv", tmp_path / "qc.json"
    merge(manifest, markers, tmp_path, output, coverage, qc)
    report = json.loads(qc.read_text())
    assert report["unrepresented_species"] == ["E"]
    assert report["gene_tree_count_distribution"] == [{"gene_trees": 0, "species": 1},
                                                        {"gene_trees": 1, "species": 4}]
    with pytest.raises(ValueError, match="insufficient gene-tree coverage"):
        astral(output, qc, manifest, tmp_path / "tree", tmp_path / "tree.json", sys.executable, "A", 1, 1)
    assert read_tsv(coverage)[-1]["gene_trees"] == "0"
    assert read_tsv(coverage)[-1]["represented"] == "False"


@pytest.mark.parametrize("present", [14, 4])
def test_post_alignment_coverage_is_diagnostic_without_overall_or_order_filter(tmp_path, present):
    manifest, markers = tmp_path / "species.tsv", tmp_path / "markers.tsv"
    species = [f"plant_{i:02d}" for i in range(24)]
    # Even if a legacy manifest has clade labels, absence of the ten-species
    # order from one gene must not exclude that gene. Four of 24 species is also
    # sufficient for retention. A second gene covers all selected species.
    write_tsv(manifest, ["species", "clade"],
              [{"species": n, "clade": "X" if i < 14 else "Y"} for i, n in enumerate(species)])
    write_tsv(markers, ["marker"], [{"marker": "1at1"}, {"marker": "2at1"}])
    for marker, names in [("1at1", species[:present]), ("2at1", species)]:
        (tmp_path / f"{marker}.nwk").write_text("(" + ",".join(f"{n}:1" for n in names) + ");\n")
        write_json(tmp_path / f"{marker}.json", {"status": "retained", "sites": 250})
    merge(manifest, markers, tmp_path, tmp_path / "trees", tmp_path / "coverage", tmp_path / "qc")
    report = json.loads((tmp_path / "qc").read_text())
    assert [r["marker"] for r in report["retained"]] == ["1at1", "2at1"]
    assert report["status"] == "retained"
    assert report["unrepresented_species"] == []
    assert all(r["represented"] == "True" for r in read_tsv(tmp_path / "coverage"))
    assert report["gene_tree_count_distribution"] == [
        {"gene_trees": 1, "species": 24 - present}, {"gene_trees": 2, "species": present}]
    assert report["excluded"] == []
    assert report["retained"][0]["occupancy"] == pytest.approx(present / 24)
    assert report["retained"][1]["occupancy"] == 1


def test_renamed_binary_cannot_pass_int128_check(tmp_path):
    binary = tmp_path / "bin/astral4_int128"
    binary.parent.mkdir()
    binary.write_text("#!/bin/sh\nexit 1\n")
    binary.chmod(0o755)
    manifest = tmp_path / "species.tsv"
    write_tsv(manifest, ["species"], [{"species": f"s{i}"} for i in range(5001)])
    with pytest.raises(ValueError, match="build provenance missing"):
        astral("unused", "unused", manifest, "unused", "unused", str(binary), "s0", 1, 1)


def test_real_lsd2_outputs_time_units_and_honors_calibration(tmp_path):
    lsd2 = os.environ.get("LSD2_BIN") or shutil.which("lsd2")
    if not lsd2:
        pytest.skip("set LSD2_BIN for real dating integration")
    tree, provenance = tmp_path / "input.nwk", tmp_path / "input.json"
    tree.write_text("(A:0.1,(B:0.08,(C:0.04,D:0.04):0.04):0.02);\n")
    write_json(provenance, {"branch_length_unit": "substitutions_per_site", "outgroup": "A", "mean_gene_length": 250, "total_gene_sites": 750})
    calibrations = tmp_path / "calibrations.tsv"
    write_tsv(calibrations, ["taxa", "min_age_ma", "max_age_ma", "source"],
              [{"taxa": "A,B", "min_age_ma": 100, "max_age_ma": 100, "source": "synthetic test only"}])
    date(tree, provenance, calibrations, tmp_path / "dated", lsd2, settings={"variance": 1})
    dated = read_tree(tmp_path / "dated/species_tree.dated.nwk")
    assert dated.get_distance("A", "B") == pytest.approx(200)
    assert dated.get_distance("C", "D") == pytest.approx(80, rel=0.01)
    report = json.loads((tmp_path / "dated/provenance.json").read_text())
    assert report["root_age_ma"] == pytest.approx(100)
    assert report["concatenation_used"] is False


def phylogeny_inputs(tmp_path):
    source = tmp_path / "input"
    source.mkdir()
    cds, busco = source / "cds", source / "busco"
    cds.mkdir(); busco.mkdir()
    rng = random.Random(73)
    amino = "ACDEFGHIKLMNPQRSTVWY"
    codons = ["GCT", "TGT", "GAT", "GAA", "TTT", "GGT", "CAT", "ATT", "AAA", "CTT",
              "ATG", "AAT", "CCT", "CAA", "CGT", "TCT", "ACT", "GTT", "TGG", "TAT"]
    code = dict(zip(amino, codons))
    ancestors = ["".join(rng.choice(amino) for _ in range(250)) for _ in range(3)]
    metadata, summaries, taxonomy = [], [], source / "taxa.sqlite"
    species = [f"Plant_{name}" for name in ["alpha", "beta", "gamma", "delta", "epsilon", "zeta"]]
    for i, name in enumerate(species):
        with gzip.open(cds / f"{name}_longestCDS.fa.gz", "wt") as handle:
            for marker, ancestor in enumerate(ancestors, 1):
                protein = list(ancestor)
                for j in range(i * 7, i * 7 + 35):
                    protein[j] = amino[(amino.index(protein[j]) + i + 1) % 20]
                handle.write(f'>{name}_g{marker}\n' + "".join(code[c] for c in protein) + "TAA\n")
        table(busco / f"{name}.busco.full.tsv",
              [f"{m}at1\tComplete\t{name}_g{m}:0-749\t200\t250" for m in range(1, 4)])
        metadata.append({"scientific_name": name.replace("_", " "), "run": f"R{i}", "taxid": str(42+i)})
        summaries.append({"Species": name.replace("_", " "), "busco_cds_single": 3, "busco_cds_duplicated": 0,
                          "busco_cds_fragmented": 0, "busco_cds_missing": 0, "busco_cds_total": 3})
        write_tsv(source / "quant" / name / f"R{i}" / f"R{i}_abundance.tsv", ["target_id", "tpm"],
                  [{"target_id": f"{name}_g1", "tpm": 100}])
    write_tsv(source / "metadata.tsv", list(metadata[0]), metadata)
    write_tsv(source / "busco.tsv", list(summaries[0]), summaries)
    with sqlite3.connect(taxonomy) as db:
        db.executescript("""CREATE TABLE species (taxid INTEGER PRIMARY KEY, parent INTEGER, spname TEXT, common TEXT, rank TEXT, track TEXT);
            CREATE TABLE merged (taxid_old INTEGER, taxid_new INTEGER); CREATE TABLE stats (version INTEGER);
            INSERT INTO stats VALUES (2); INSERT INTO species VALUES (1,1,'root','','no rank','1');
            INSERT INTO species VALUES (2,1,'Plants','','kingdom','2,1');""")
        db.executemany("INSERT INTO species VALUES (?,2,?,'','species',?)",
                       [(42+i, name.replace("_", " "), f"{42+i},2,1") for i, name in enumerate(species)])
    return source, species


def test_real_phylogeny_workflow_and_unchanged_rerun(tmp_path, command_environment, workflow_project):
    snakemake = os.environ.get("SNAKEMAKE_BIN") or shutil.which("snakemake")
    famsa = os.environ.get("FAMSA_BIN") or shutil.which("famsa")
    vft = os.environ.get("VERYFASTTREE_BIN") or shutil.which("VeryFastTree")
    astral4 = ROOT / "resources/phylogeny_tools/bin/astral4_int128"
    if not all([snakemake, famsa, vft]) or not astral4.is_file():
        pytest.skip("set phylogeny tool paths and run prepare_phylogeny_tools.py for workflow integration")
    trimal = trimal_binary()
    lsd2 = os.environ.get("LSD2_BIN") or shutil.which("lsd2")
    commands = {"python": sys.executable, "famsa": famsa, "trimal": trimal, "VeryFastTree": vft}
    if lsd2:
        commands["lsd2"] = lsd2
    env = command_environment(commands)
    source, species = phylogeny_inputs(tmp_path)
    cfg = {"analysis": "test", "inputs": {"metadata": str(source / "metadata.tsv"), "busco": str(source / "busco.tsv"),
           "cds_dir": str(source / "cds"), "quant_dir": str(source / "quant")},
           "taxonomy": {"source": str(source / "taxa.sqlite")},
           "phylogeny": {"busco_full_dir": str(source / "busco"), "outgroup": species[0],
               "max_markers": 3,
               "align_threads": 1, "tree_threads": 1, "astral_threads": 2, "astral_mem_gb": 4}}
    conda_prefix = os.environ.get("PHYLOGENY_CONDA_PREFIX")
    config = tmp_path / "config.yaml"
    config.write_text(yaml.safe_dump(cfg))
    argv = [snakemake, "--snakefile", str(ROOT / "workflow/Snakefile"), "--configfile", str(config),
            "--cores", "2", "--resources", "mem_mb=8000", "--", "phylogeny"]
    if conda_prefix:
        argv[1:1] = ["--use-conda", "--conda-prefix", conda_prefix]
    def run():
        result = subprocess.run(argv, cwd=workflow_project, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                env=env)
        assert result.returncode == 0, result.stdout + "\n" + "\n".join(
            p.read_text()[-4000:] for p in (tmp_path / "logs").rglob("*.log"))
        return result.stdout
    run()
    out = tmp_path / "results/test/phylogeny/all"
    read_tree(out / "species_tree.nwk", species)
    assert all(int(r["gene_trees"]) >= 1 for r in read_tsv(out / "species_coverage.tsv"))
    report = json.loads((out / "species_tree.json").read_text())
    assert report["branch_length_method"] == "CASTLES-II"
    assert report["concatenation_used"] is False
    timestamp = (out / "species_tree.nwk").stat().st_mtime_ns
    assert "Nothing to be done" in run()
    assert (out / "species_tree.nwk").stat().st_mtime_ns == timestamp
    assert not (tmp_path / "results/test/orthogroups/mapping").exists()
    assert not (tmp_path / "results/test/orthogroups/expression").exists()
    upstream = {p: p.stat().st_mtime_ns for folder in [out / "species", out / "alignments/raw"]
                for p in folder.iterdir()}
    cfg["phylogeny"]["trimal_mode"] = "automated1"
    config.write_text(yaml.safe_dump(cfg))
    run()
    assert all(p.stat().st_mtime_ns == mtime for p, mtime in upstream.items())
    assert all(json.loads(p.read_text())["mode"] == "automated1" for p in (out / "alignments").glob("*.json"))
    timestamp = (out / "species_tree.nwk").stat().st_mtime_ns
    assert "Nothing to be done" in run()
    # The dating target reuses the inferred species tree and does not restart loci.
    if lsd2:
        calibrations = tmp_path / "calibrations.tsv"
        write_tsv(calibrations, ["taxa", "min_age_ma", "max_age_ma", "source"],
                  [{"taxa": ",".join(species[:2]), "min_age_ma": 100, "max_age_ma": 100,
                    "source": "synthetic workflow test only"}])
        cfg["phylogeny"]["dating"] = {"calibrations": str(calibrations), "mem_gb": 4,
                                        "lsd2": {"variance": 1}}
        config.write_text(yaml.safe_dump(cfg))
        argv[-1] = "timetree"
        run()
        read_tree(out / "dating/species_tree.dated.nwk", species)
        assert (out / "species_tree.nwk").stat().st_mtime_ns == timestamp
        dated_timestamp = (out / "dating/species_tree.dated.nwk").stat().st_mtime_ns
        assert "Nothing to be done" in run()
        assert (out / "dating/species_tree.dated.nwk").stat().st_mtime_ns == dated_timestamp
        retained = {p: p.stat().st_mtime_ns for folder in [out / "gene_trees", out / "alignments/raw"]
                    for p in folder.iterdir()}
        (out / "dating/species_tree.dated.nwk").unlink()
        run()
        assert all(p.stat().st_mtime_ns == stamp for p, stamp in retained.items())
        cfg["phylogeny"]["dating"]["lsd2"]["variance"] = 0
        config.write_text(yaml.safe_dump(cfg))
        run()
        assert json.loads((out / "dating/provenance.json").read_text())["settings"]["variance"] == 0
        assert (out / "species_tree.nwk").stat().st_mtime_ns == timestamp
        # Switch from manual bounds to the automatic TimeTree branch. Exercise
        # the full rule graph offline with a recorded synthetic API response.
        from timetree_calibrations import cached_response
        from test_timetree_calibrations import fake_fetch, payload
        cache = workflow_project / "resources/timetree_cache"
        cached_response(range(42, 48), cache, delay=0,
                        backend=(fake_fetch(payload(range(42, 48))), {"synthetic": True}))
        cfg["phylogeny"]["dating"].update(calibration_source="timetree", timetree={
            "max_representatives": 6, "max_queries": 1,
            "offline": True})
        config.write_text(yaml.safe_dump(cfg))
        run()
        assert json.loads((out / "timetree/provenance.json").read_text())["status"] == "ready"
        assert len(read_tsv(out / "timetree/calibrations.tsv")) == 1
        assert (out / "species_tree.nwk").stat().st_mtime_ns == timestamp
        assert "Nothing to be done" in run()
