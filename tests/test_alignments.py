"""All-copy OG alignments, residue preservation and real workflow resumption."""
import json
import os
from pathlib import Path
import shutil
import sqlite3
import subprocess
import sys

import pytest
import yaml

from align_orthogroups import align, collect, fasta_records, finish
from common import read_tsv, write_tsv


ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def collected_inputs(tmp_path):
    samples = tmp_path / "samples.tsv"
    write_tsv(samples, ["species", "odb_species", "run"], [
        {"species": "Alpha", "odb_species": "Alpha", "run": "A1"},
        {"species": "Alpha", "odb_species": "Alpha", "run": "A2"},
        {"species": "Beta-X", "odb_species": "Beta_X", "run": "B1"},
    ])
    proteins = tmp_path / "proteins"
    proteins.mkdir()
    (proteins / "Alpha_protein.fa").write_text(
        ">a1\nMACDEFGHIKLMNPQRSTVWY*\n>a2\nMACDEFGHIKLMNPQRSTVWY*\n>unused\nMXXX*\n")
    (proteins / "Beta_X_protein.fa").write_text(
        ">b1\nMACDXXGHIKLMNPQRSTVWY*\n>b2\nMAUBZOJ*ACD\n")
    database = tmp_path / "mappings.sqlite"
    with sqlite3.connect(database) as db:
        db.executescript("""
            CREATE TABLE genes (query TEXT PRIMARY KEY, species TEXT NOT NULL);
            CREATE TABLE mappings (query TEXT, og TEXT, PRIMARY KEY (query, og));
            INSERT INTO genes VALUES ('a1','Alpha'), ('a2','Alpha'), ('unused','Alpha'),
                                     ('b1','Beta-X'), ('b2','Beta-X');
            INSERT INTO mappings VALUES ('a1','OG1'), ('a2','OG1'), ('b1','OG1'),
                                        ('b2','OG2'), ('a2','OG3');
        """)
    inputs = tmp_path / "inputs"
    collect(samples, database, proteins, inputs)
    return samples, database, proteins, inputs


def test_collect_preserves_all_copies_assignments_and_species(collected_inputs):
    _, _, _, inputs = collected_inputs
    assert set(p.name for p in inputs.glob("*.faa")) == {"OG1.faa", "OG2.faa", "OG3.faa"}
    og1 = dict(fasta_records(inputs / "OG1.faa"))
    assert set(og1) == {"a1", "a2", "b1"}
    assert og1["a1"] == og1["a2"]
    assert dict(fasta_records(inputs / "OG3.faa")) == {"a2": og1["a2"]}
    members = read_tsv(inputs / "members.tsv")
    assert len(members) == 5  # repeated runs do not duplicate sequences
    assert next(r for r in members if r["gene_id"] == "b2")["species"] == "Beta-X"
    assert {r["orthogroup"] for r in members if r["gene_id"] == "a2"} == {"OG1", "OG3"}
    assert "unused" not in {r["gene_id"] for r in members}


def test_collection_replaces_removed_ogs_and_handles_no_mappings(collected_inputs):
    samples, database, proteins, inputs = collected_inputs
    with sqlite3.connect(database) as db:
        db.execute("DELETE FROM mappings WHERE og = 'OG3'")
    collect(samples, database, proteins, inputs)
    assert not (inputs / "OG3.faa").exists()
    with sqlite3.connect(database) as db:
        db.execute("DELETE FROM mappings")
    collect(samples, database, proteins, inputs)
    assert not list(inputs.glob("*.faa"))
    assert read_tsv(inputs / "members.tsv") == []
    out = inputs.parent / "alignments"
    finish(inputs, out, inputs.parent / "reports")
    assert json.loads((out / "provenance.json").read_text())["alignments"] == []


@pytest.mark.parametrize("kind", ["missing", "duplicate", "invalid", "unsafe_og", "wrong_species"])
def test_collection_rejects_inconsistent_inputs_without_replacing_checkpoint(collected_inputs, kind):
    samples, database, proteins, inputs = collected_inputs
    previous = (inputs / "members.tsv").read_bytes()
    if kind == "missing":
        (proteins / "Beta_X_protein.fa").write_text(">b1\nMXX*\n")
    elif kind == "duplicate":
        with open(proteins / "Beta_X_protein.fa", "a") as handle:
            handle.write(">b2\nMAA*\n")
    elif kind == "invalid":
        (proteins / "Beta_X_protein.fa").write_text(">b1\nM-AA\n>b2\nMAA*\n")
    else:
        with sqlite3.connect(database) as db:
            db.execute("UPDATE mappings SET og = '../escape' WHERE og = 'OG3'" if kind == "unsafe_og"
                       else "UPDATE genes SET species = 'Other' WHERE query = 'b2'")
    with pytest.raises(ValueError):
        collect(samples, database, proteins, inputs)
    assert (inputs / "members.tsv").read_bytes() == previous


def test_singleton_is_saved_without_running_famsa(collected_inputs, tmp_path):
    *_, inputs = collected_inputs
    output, report = tmp_path / "out.faa", tmp_path / "out.json"
    align(inputs / "OG2.faa", output, report, command="nonexistent-famsa")
    assert output.read_bytes() == (inputs / "OG2.faa").read_bytes()
    assert json.loads(report.read_text())["method"] == "singleton"


def test_collection_reopens_many_og_files_without_losing_copies(collected_inputs):
    samples, database, proteins, inputs = collected_inputs
    groups = [f"many{i:03d}" for i in range(70)]
    with sqlite3.connect(database) as db:
        db.executemany("INSERT INTO mappings VALUES (?, ?)",
                       [(gene, og) for gene in ["a1", "a2"] for og in groups])
    collect(samples, database, proteins, inputs)
    for og in groups:
        assert [name for name, _ in fasta_records(inputs / f"{og}.faa")] == ["a1", "a2"]


@pytest.mark.parametrize("replacement, message", [
    (">a1\nMAA\n", "gene ID set"),
    (">a1\nMAA\n>a2\nMAAA\n", "unequal sequence lengths"),
    (">a1\nM-A\n>a2\nMAA\n", "changed input residues"),
])
def test_bad_aligner_output_is_not_published(tmp_path, replacement, message):
    fasta = tmp_path / "in.faa"
    fasta.write_text(">a1\nMAA\n>a2\nMAA\n")
    command = tmp_path / "fake_famsa"
    command.write_text(f"#!{sys.executable}\nfrom pathlib import Path\nimport sys\n"
                       f"Path(sys.argv[-1]).write_text({replacement!r})\n")
    command.chmod(0o755)
    output = tmp_path / "out.faa"
    output.write_text("previous result\n")
    with pytest.raises(ValueError, match=message):
        align(fasta, output, tmp_path / "out.json", command=str(command))
    assert output.read_text() == "previous result\n"
    assert not (tmp_path / "out.json").exists()


def famsa_binary():
    binary = os.environ.get("FAMSA_BIN") or shutil.which("famsa")
    if not binary:
        pytest.skip("set FAMSA_BIN for real alignment tests")
    return binary


def test_real_famsa_preserves_residues_and_finish_prunes_removed_ogs(collected_inputs, tmp_path):
    samples, database, proteins, inputs = collected_inputs
    output, reports = tmp_path / "alignments", tmp_path / "reports"
    binary = famsa_binary()
    # Different lengths, stop/ambiguity symbols, identical copies and singletons.
    with open(proteins / "Alpha_protein.fa", "a") as handle:
        handle.write(">special\nMAUBZOJ*ACD\n>unknown\nXXXXX\n")
    with sqlite3.connect(database) as db:
        db.executemany("INSERT INTO genes VALUES (?, 'Alpha')", [("special",), ("unknown",)])
        db.executemany("INSERT INTO mappings VALUES (?, 'OG1')", [("special",), ("unknown",)])
    collect(samples, database, proteins, inputs)
    for fasta in sorted(inputs.glob("*.faa")):
        align(fasta, output / fasta.name, reports / f"{fasta.stem}.json", command=binary)
        result = dict(fasta_records(output / fasta.name))
        assert len({len(seq) for seq in result.values()}) == 1
        assert {g: s.replace("-", "") for g, s in result.items()} == dict(fasta_records(fasta))
    (output / "obsolete.faa").write_text(">old\nMAA\n")
    (output / "notes.txt").write_text("keep\n")
    finish(inputs, output, reports)
    assert not (output / "obsolete.faa").exists()
    assert (output / "notes.txt").read_text() == "keep\n"
    assert (output / "members.tsv").read_bytes() == (inputs / "members.tsv").read_bytes()
    report = json.loads((output / "provenance.json").read_text())
    assert {r["orthogroup"] for r in report["alignments"]} == {"OG1", "OG2", "OG3"}
    with open(output / "OG1.faa", "a") as handle:
        handle.write(">unexpected\nMAA\n")
    with pytest.raises(ValueError, match="changed after its job completed"):
        finish(inputs, output, reports)


def test_real_alignment_workflow_resume_updates_and_opt_in(
    tiny_inputs, fake_odb, frozen_reference, tmp_path, command_environment, workflow_project,
):
    snakemake = os.environ.get("SNAKEMAKE_BIN") or shutil.which("snakemake")
    seqkit = os.environ.get("SEQKIT_BIN") or shutil.which("seqkit")
    famsa = famsa_binary()
    if not snakemake or not seqkit:
        pytest.skip("set SNAKEMAKE_BIN and SEQKIT_BIN for workflow integration")
    # Include copies, repeated OG assignments, a singleton and an unmapped gene.
    fake_odb.write_text(fake_odb.read_text().replace(
        "line.endswith(('_g1','_g2'))", "line.endswith(('_g1','_g2','_g3'))").replace(
        "lines.extend([gene + '\\tOG' + gene[-1] + '\\n'] * 2)",
        "groups = ['OG1'] if gene.endswith('_g1') else ['OG1', 'OG2'] if gene.endswith('_g2') "
        "else ['OG3'] if gene.startswith('Beta_') else []\n"
        "                for og in groups:\n"
        "                    lines.extend([gene + '\\t' + og + '\\n'] * 2)"))
    reference = workflow_project / "resources/orthodb/v12_3193"
    reference.parent.mkdir(parents=True)
    reference.symlink_to(frozen_reference, target_is_directory=True)
    config = {
        "analysis": "test", "inputs": {k: tiny_inputs[k] for k in ["metadata", "busco", "cds_dir", "quant_dir"]},
        "taxonomy": {"source": tiny_inputs["taxonomy_db"]},
        "odb": {"chunk_size": 1, "threads": 1, "batch_size": 1,
                "mem_gb": 3, "min_free_gb": 0, "allow_nonlocal": True},
        "alignment": {"enabled": False, "threads": 1, "mem_gb": 2},
    }
    configfile = tmp_path / "override.yaml"
    configfile.write_text(yaml.safe_dump(config))
    # Wrap the real aligner to record calls and simulate one interrupted OG job.
    wrapper = tmp_path / "famsa_wrapper"
    wrapper.write_text(f"#!{sys.executable}\nimport os, sys\nfrom pathlib import Path\n"
                       "og = Path(sys.argv[-2]).stem\n"
                       "with open(os.environ['FAMSA_EVENTS'], 'a') as f: f.write(og + '\\n')\n"
                       "flag = Path(os.environ['FAMSA_FAIL_ONCE'])\n"
                       "if og == 'OG2' and flag.exists():\n"
                       "    flag.unlink()\n    raise SystemExit(23)\n"
                       f"os.execv({str(famsa)!r}, [{str(famsa)!r}] + sys.argv[1:])\n")
    wrapper.chmod(0o755)
    events, fail_once, odb_events = tmp_path / "align_events", tmp_path / "fail_once", tmp_path / "odb_events"
    env = {**command_environment({"python": sys.executable, "seqkit": seqkit,
                                 "ODB-mapper": fake_odb, "famsa": wrapper}),
           "FAMSA_EVENTS": str(events), "FAMSA_FAIL_ONCE": str(fail_once), "FAKE_ODB_LOG": str(odb_events)}
    base = [snakemake, "--snakefile", str(ROOT / "workflow/Snakefile"), "--configfile", str(configfile),
            "--cores", "2", "--resources", "mem_mb=16000"]

    def execute(options=(), targets=("alignments",), fail=False):
        result = subprocess.run(base + list(options) + ["--"] + list(targets), cwd=workflow_project, env=env,
                                text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=180)
        if fail:
            assert result.returncode != 0, result.stdout
        elif result.returncode:
            logs = "\n".join(f"{p}:\n{p.read_text()}" for p in (tmp_path / "logs").rglob("*.log"))
            pytest.fail(result.stdout + "\n" + logs)
        return result.stdout

    execute()
    out = tmp_path / "results/test/alignments"
    assert {p.name for p in out.glob("*.faa")} == {"OG1.faa", "OG2.faa", "OG3.faa"}
    assert len(list(fasta_records(out / "OG1.faa"))) == 4
    assert len(list(fasta_records(out / "OG2.faa"))) == 2
    assert len(list(fasta_records(out / "OG3.faa"))) == 1
    assert len(read_tsv(out / "members.tsv")) == 7
    assert sorted(events.read_text().splitlines()) == ["OG1", "OG2"]
    assert len(odb_events.read_text().splitlines()) == 2
    assert not (out.parent / "tpm").exists()
    assert not (out.parent / "phylogeny").exists()
    times = {p.name: p.stat().st_mtime_ns for p in out.glob("*.faa")}
    assert "Nothing to be done" in execute()

    # Missing output + interrupted rerun: only the affected OG is recomputed.
    (out / "OG2.faa").unlink()
    fail_once.touch()
    execute(fail=True)
    assert (out / "OG1.faa").stat().st_mtime_ns == times["OG1.faa"]
    execute(options=("--rerun-incomplete",))
    assert events.read_text().splitlines().count("OG1") == 1
    assert events.read_text().splitlines().count("OG2") == 3

    # Quantification changes and TPM ambiguity policy do not affect alignments.
    abundance = Path(tiny_inputs["quant_dir"]) / "Alpha_plant/A1/A1_abundance.tsv"
    values = read_tsv(abundance)
    values[0]["tpm"] = 40
    write_tsv(abundance, list(values[0]), values)
    assert "Nothing to be done" in execute()
    config["alignment"]["enabled"] = True
    config["tpm"] = {"multimap": "drop"}
    configfile.write_text(yaml.safe_dump(config))
    times = {p.name: p.stat().st_mtime_ns for p in out.glob("*.faa")}
    (out / "OG2.faa").unlink()
    execute(targets=())
    assert (out.parent / "tpm/tpm.tsv").exists()
    assert (out / "OG2.faa").exists()
    assert all((out / name).stat().st_mtime_ns == t for name, t in times.items() if name != "OG2.faa")
    assert len(read_tsv(out / "members.tsv")) == 7  # multi-OG copies still present
    assert "Nothing to be done" in execute(targets=())

    # Species changes rebuild membership and remove no-longer-observed OGs.
    subset = tmp_path / "subset.txt"
    subset.write_text("Alpha_plant\n")
    config["selection"] = {"species_list": str(subset)}
    configfile.write_text(yaml.safe_dump(config))
    execute()
    assert {p.name for p in out.glob("*.faa")} == {"OG1.faa", "OG2.faa"}
    assert {r["species"] for r in read_tsv(out / "members.tsv")} == {"Alpha_plant"}
    assert len(list(fasta_records(out / "OG1.faa"))) == 2
    assert len(list(fasta_records(out / "OG2.faa"))) == 1
    assert "Nothing to be done" in execute()


@pytest.mark.parametrize("settings,message", [
    ({"enabled": "yes"}, "alignment.enabled must be true or false"),
    ({"threads": 0}, "alignment.threads must be a positive integer"),
    ({"mem_gb": True}, "alignment.mem_gb must be a positive integer"),
    ({"command": "custom"}, "fixed by the workflow"),
])
def test_invalid_alignment_config(tmp_path, workflow_project, settings, message):
    snakemake = os.environ.get("SNAKEMAKE_BIN") or shutil.which("snakemake")
    if not snakemake:
        pytest.skip("Snakemake is unavailable")
    config = tmp_path / "override.yaml"
    config.write_text(yaml.safe_dump({"alignment": settings}))
    result = subprocess.run([snakemake, "--snakefile", str(ROOT / "workflow/Snakefile"),
                             "--configfile", str(config), "--dry-run"], cwd=workflow_project,
                            text=True, capture_output=True, timeout=60,
                            env={**os.environ, "XDG_CACHE_HOME": str(tmp_path / "cache")})
    assert result.returncode != 0
    assert message in result.stdout + result.stderr
