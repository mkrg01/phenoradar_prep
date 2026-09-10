"""Exercise KEGG's real workflow graph with a test-only KofamScan substitute."""
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from common import read_tsv, write_tsv
from prepare_kegg_reference import prepare as prepare_reference
from publish_kegg_reference import publish
from verify_kegg_reference import verify


ROOT = Path(__file__).resolve().parents[1]


def small_reference(tmp_path):
    source = tmp_path / "kofam_source"
    source.mkdir()
    for ko in ["K00001", "K00002", "K00003"]:
        (source / f"{ko}.hmm").write_text(f"HMMER3/f\nNAME  {ko}\nALPH  amino\n//\n")
    ko_list = tmp_path / "ko_list"
    header = ["knum", "threshold", "score_type", "profile_type", "F-measure", "nseq", "nseq_used",
              "alen", "mlen", "eff_nseq", "re/pos", "definition"]
    ko_list.write_text("\t".join(header) + "\n" + "".join(
        f"{ko}\t10\tfull\tall\t1\t10\t10\t100\t100\t10\t1\ttest\n"
        for ko in ["K00001", "K00002", "K00003"]))
    modules, pathways = tmp_path / "module_links", tmp_path / "pathway_links"
    modules.write_text("ko:K00001\tmd:M00001\nko:K00002\tmd:M00001\n")
    pathways.write_text("ko:K00001\tpath:ko00010\nko:K00003\tpath:map00020\n")
    root = tmp_path / "kegg_reference"
    prepare_reference(source, ko_list, root, modules, pathways, release="test-fixture")
    return root


def fake_kofam_command(tmp_path):
    command = tmp_path / "fake_exec_annotation"
    command.write_text('''#!/usr/bin/env python3
import os, sys
from pathlib import Path
args = sys.argv[1:]
if '--help' in args:
    print('fake KofamScan for workflow tests')
    raise SystemExit(0)
output = Path(args[args.index('-o') + 1])
query = Path(args[-1])
genes = [line[1:].split()[0] for line in query.read_text().splitlines() if line.startswith('>')]
with open(os.environ['FAKE_KOFAM_LOG'], 'a') as log:
    log.write(genes[0] + '\\n')
lines = ['#\\tgene name\\tKO\\tthrshld\\tscore\\tE-value\\tKO definition\\n']
for gene in genes:
    if gene.endswith('_g1'):
        lines.append('*\\t' + gene + '\\tK00001\\t10.00\\t20.0\\t1e-8\\t"unique"\\n')
    elif gene.endswith('_g2'):
        for ko in ['K00002', 'K00003']:
            lines.append('*\\t' + gene + '\\t' + ko + '\\t10.00\\t20.0\\t1e-8\\t"ambiguous"\\n')
    else:
        lines.append('\\t' + gene + '\\tK00001\\t10.00\\t5.0\\t0.01\\t"weak"\\n')
output.write_text(''.join(lines))
''')
    command.chmod(0o755)
    return command


def test_publish_preserves_frozen_reference(tmp_path):
    reference = small_reference(tmp_path)
    original = verify(reference / "reference.json")["reference_id"]
    output = tmp_path / "published"
    publish(reference / "reference.json", output / "qc.json", output / "modules.tsv", output / "pathways.tsv")
    assert json.loads((output / "qc.json").read_text())["reference_id"] == original
    with pytest.raises(ValueError, match="frozen reference"):
        publish(reference / "reference.json", reference / "extra.json", output / "modules.tsv", output / "pathways.tsv")
    with pytest.raises(ValueError, match="distinct"):
        publish(reference / "reference.json", output / "same", output / "same", output / "pathways.tsv")
    assert verify(reference / "reference.json")["reference_id"] == original


def test_kegg_standalone_incremental_and_opt_in_full(tiny_inputs, fake_odb, frozen_reference, tmp_path, command_environment):
    snakemake = os.environ.get("SNAKEMAKE_BIN") or shutil.which("snakemake")
    seqkit = os.environ.get("SEQKIT_BIN") or shutil.which("seqkit")
    if not snakemake or not seqkit:
        pytest.skip("set SNAKEMAKE_BIN and SEQKIT_BIN for workflow integration")
    reference = small_reference(tmp_path)
    command = fake_kofam_command(tmp_path)
    config = {
        "analysis": "test", "inputs": {k: tiny_inputs[k] for k in ["metadata", "busco", "cds_dir", "quant_dir"]},
        "taxonomy": {"database": tiny_inputs["taxonomy_db"]},
        "paths": {"results": str(tmp_path / "results"), "work": str(tmp_path / "work"), "logs": str(tmp_path / "logs")},
        "odb": {"reference_dir": str(frozen_reference), "chunk_size": 1, "threads": 1, "batch_size": 1,
                "mem_gb": 3, "min_free_gb": 0, "allow_nonlocal": True},
        "kegg": {"enabled": False, "reference_dir": str(reference),
                 "threads": 1, "mem_gb": 2},
    }
    configfile = tmp_path / "config.yaml"
    configfile.write_text(yaml.safe_dump(config))
    events = tmp_path / "kofam_events.txt"
    odb_events = tmp_path / "odb_events.txt"
    env = {**command_environment({"python": sys.executable, "seqkit": seqkit,
                                 "ODB-mapper": fake_odb, "exec_annotation": command}),
           "FAKE_KOFAM_LOG": str(events), "FAKE_ODB_LOG": str(odb_events)}
    base = [snakemake, "--snakefile", str(ROOT / "workflow/Snakefile"), "--configfile", str(configfile),
            "--cores", "2", "--resources", "mem_mb=16000"]

    def execute(options=(), targets=("kegg",)):
        result = subprocess.run(base + list(options) + ["--"] + list(targets), cwd=ROOT, env=env,
                                text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if result.returncode:
            logs = "\n".join(f"{p}:\n{p.read_text()}" for p in (tmp_path / "logs").rglob("*.log"))
            pytest.fail(result.stdout + "\n" + logs)
        return result.stdout

    execute()
    out = tmp_path / "results/test/kegg"
    assert not odb_events.exists()  # A standalone KEGG target never maps to ODB.
    assert len(events.read_text().splitlines()) == 2  # two species, three runs
    benchmarks = list((out / "species").glob("*/benchmark.tsv"))
    assert len(benchmarks) == 2
    assert all(float(read_tsv(path)[0]["s"]) >= 0 for path in benchmarks)
    rows = read_tsv(out / "ko_tpm_sum.tsv")
    assert {(r["run"], r["ko"]): float(r["tpm_sum"]) for r in rows} == {
        (run, ko): value for run, unique, multi in [("A1", 20, 30), ("A2", 80, 10), ("B1", 20, 30)]
        for ko, value in [("K00001", unique), ("K00002", multi), ("K00003", multi)]}
    assert len(read_tsv(out / "genes.tsv")) == 6
    assert {r["assignment_status"] for r in read_tsv(out / "genes.tsv")} == {"unique", "ambiguous", "below_threshold"}
    qc = {r["run"]: r for r in read_tsv(out / "mapping_qc.tsv")}
    assert qc["A1"]["ambiguity"] == "duplicate"
    assert float(qc["A1"]["retained_tpm_fraction"]) == 0.5
    assert float(qc["A1"]["retained_tpm"]) == 50
    assert float(qc["A1"]["ko_tpm_sum"]) == 80
    assert len(read_tsv(out / "ko_support.tsv")) == 9
    assert read_tsv(out / "ko_modules.tsv")[0] == {"ko": "K00001", "module": "M00001"}
    assert "Nothing to be done" in execute(["--dry-run"])

    abundance = Path(tiny_inputs["quant_dir"]) / "Alpha_plant/A1/A1_abundance.tsv"
    values = read_tsv(abundance)
    values[0]["tpm"] = 40
    write_tsv(abundance, list(values[0]), values)
    dry = execute(["--dry-run"])
    assert "rule aggregate_ko_tpm:" in dry
    assert "rule annotate_kofam:" not in dry
    execute()
    assert len(events.read_text().splitlines()) == 2
    assert float(next(r for r in read_tsv(out / "ko_tpm_sum.tsv")
                      if r["run"] == "A1" and r["ko"] == "K00001")["tpm_sum"]) == 40

    # Aggregation policy changes reuse the species annotations.
    annotation_times = {p: p.stat().st_mtime_ns for p in (out / "species").glob("*/provenance.json")}
    for ambiguity, expected_kos in [("drop", {"K00001"}), ("duplicate", {"K00001", "K00002", "K00003"})]:
        config["kegg"]["ambiguity"] = ambiguity
        configfile.write_text(yaml.safe_dump(config))
        dry = execute(["--dry-run"])
        assert "rule aggregate_ko_tpm:" in dry
        assert "rule annotate_kofam:" not in dry
        execute()
        assert {r["ko"] for r in read_tsv(out / "ko_tpm_sum.tsv")} == expected_kos
        assert len(events.read_text().splitlines()) == 2
        assert all(p.stat().st_mtime_ns == timestamp for p, timestamp in annotation_times.items())

    subset = tmp_path / "subset.txt"
    subset.write_text("Beta_sp-X\n")
    config["selection"] = {"species_list": str(subset)}
    configfile.write_text(yaml.safe_dump(config))
    execute()
    assert [r["run"] for r in read_tsv(out / "ko_tpm_sum.tsv")] == ["B1"] * 3
    assert {r["species"] for r in read_tsv(out / "genes.tsv")} == {"Beta_sp-X"}
    assert len(events.read_text().splitlines()) == 2

    config["kegg"]["enabled"] = True
    configfile.write_text(yaml.safe_dump(config))
    execute(targets=())
    assert [r["run"] for r in read_tsv(out.parent / "tpm/tpm.tsv")] == ["B1", "B1"]
    assert len(odb_events.read_text().splitlines()) == 1
    assert "Nothing to be done" in execute(["--dry-run"], targets=())


@pytest.mark.parametrize("settings,message", [
    ({"enabled": "yes"}, "kegg.enabled must be true or false"),
    ({"threads": 0}, "kegg.threads must be a positive integer"),
    ({"mem_gb": True}, "kegg.mem_gb must be a positive integer"),
    ({"ambiguity": "split"}, "kegg.ambiguity must be duplicate, drop, or error"),
    ({"command": ""}, "fixed by the workflow"),
])
def test_invalid_kegg_config_is_rejected(tmp_path, settings, message):
    snakemake = os.environ.get("SNAKEMAKE_BIN") or shutil.which("snakemake")
    if not snakemake:
        pytest.skip("Snakemake is unavailable")
    configfile = tmp_path / "config.yaml"
    configfile.write_text(yaml.safe_dump({"kegg": settings}))
    result = subprocess.run([snakemake, "--snakefile", str(ROOT / "workflow/Snakefile"),
                             "--configfile", str(configfile), "--dry-run"], cwd=ROOT, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            env={**os.environ, "XDG_CACHE_HOME": str(tmp_path / "cache")})
    assert result.returncode != 0
    assert message in result.stdout
