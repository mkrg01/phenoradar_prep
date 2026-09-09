import json
from pathlib import Path
import subprocess

import pytest

from common import read_tsv, write_json
from prepare_kegg_reference import prepare
from run_kofam import normalize_protein, parse_detail, run


HEADER = "#\tgene name\tKO\tthrshld\tscore\tE-value\tKO definition\n#\t---------\t------\t-------\t------\t---------\t-------------\n"


def detail_row(gene, ko="K00001", threshold="10.00", score="11.0", evalue="1e-10", marker="*"):
    return "\t".join([marker, gene, ko, threshold, score, evalue, '"example"']) + "\n"


@pytest.fixture
def kofam_job(tmp_path, monkeypatch):
    source = tmp_path / "sources"
    profiles = source / "profiles"
    profiles.mkdir(parents=True)
    ko_list = source / "ko_list"
    ko_list.write_text("knum\tthreshold\tscore_type\tprofile_type\tF-measure\tnseq\tnseq_used\talen\tmlen\teff_nseq\tre/pos\tdefinition\n"
                       "K00001\t10\tfull\tall\t0.9\t3\t3\t2\t2\t3\t0.5\tfirst\n"
                       "K00002\t20\tdomain\tall\t0.9\t3\t3\t2\t2\t3\t0.5\tsecond\n"
                       "K00003\t-\t-\tall\t-\t1\t1\t2\t2\t1\t0.5\tno threshold\n")
    for ko in ["K00001", "K00002", "K00003"]:
        (profiles / f"{ko}.hmm").write_text(f"HMMER3/f\nNAME  {ko}\nALPH  amino\n//\n")
    modules, pathways = source / "modules.tsv", source / "pathways.tsv"
    modules.write_text("ko:K00001\tmd:M00001\n")
    pathways.write_text("ko:K00001\tpath:map00010\n")
    reference = prepare(profiles, ko_list, tmp_path / "reference", modules, pathways, release="fixture")
    command = tmp_path / "fake kofam scan"
    command.write_text('''#!/usr/bin/env python3
import json, os, sys
from pathlib import Path
args = sys.argv[1:]
assert args[args.index('-f') + 1] == 'detail-tsv'
assert args[args.index('-T') + 1] == '1'
assert '--no-report-unannotated' in args
config = json.loads(Path(args[args.index('-c') + 1]).read_text())
assert config['ko_list'] == args[args.index('-k') + 1]
protein = Path(args[-1]).read_text()
assert '*' not in protein
with open(os.environ['FAKE_KOFAM_EVENTS'], 'a') as log:
    log.write(json.dumps({'args': args, 'protein': protein}) + '\\n')
output = Path(args[args.index('-o') + 1])
output.write_text(Path(os.environ['FAKE_KOFAM_DETAIL']).read_text())
if os.environ.get('FAKE_KOFAM_FAIL') == '1':
    print('intentional fake failure', file=sys.stderr)
    sys.exit(29)
''')
    command.chmod(0o755)
    protein = tmp_path / "input proteins.fa"
    protein.write_text(">plant_g1 original header\nMK*\n")
    detail = tmp_path / "expected_detail.tsv"
    detail.write_text(HEADER + detail_row("plant_g1"))
    events = tmp_path / "events.jsonl"
    monkeypatch.setenv("FAKE_KOFAM_EVENTS", str(events))
    monkeypatch.setenv("FAKE_KOFAM_DETAIL", str(detail))
    monkeypatch.delenv("SLURM_CPUS_PER_TASK", raising=False)
    return {"args": dict(protein=protein, species="plant", reference=reference,
                         output_dir=tmp_path / "out", work_dir=tmp_path / "work",
                         command=str(command), threads=1),
            "detail": detail, "events": events}


def test_annotations_keep_all_gene_states_and_original_identifiers(kofam_job):
    args = kofam_job["args"]
    args["protein"].write_text("".join(f">plant_g{i} original header {i}\nMK{'*' if i == 1 else ''}\n"
                                     for i in range(1, 7)))
    accepted = detail_row("plant_g1")
    kofam_job["detail"].write_text(HEADER + accepted + accepted
                                  + detail_row("plant_g1", "K00002", "20.00", "2.0", marker="")
                                  + detail_row("plant_g2")
                                  + detail_row("plant_g2", "K00002", "20.00", "21.0")
                                  + detail_row("plant_g4", "K00003", "", "50.0", marker="")
                                  + detail_row("plant_g5", score="2.0", marker="")
                                  + detail_row("plant_g6", score="10.0", marker=""))
    result = run(**args)
    genes = read_tsv(args["output_dir"] / "genes.tsv")
    assert [gene["assignment_status"] for gene in genes] == [
        "unique", "ambiguous", "unannotated", "threshold_missing", "below_threshold", "below_threshold"]
    assert [gene["selected_ko"] for gene in genes] == ["K00001", "", "", "", "", ""]
    assert len(read_tsv(args["output_dir"] / "gene_kos.tsv")) == 7
    assert genes[0]["terminal_stop_stripped"] == "1"
    event = json.loads(kofam_job["events"].read_text())
    assert ">plant_g1 original header 1\nMK\n" in event["protein"]
    assert result["identity"]["species"] == "plant"
    assert result["terminal_stop_stripped_count"] == 1
    assert (args["output_dir"] / "detail.tsv").read_text() == kofam_job["detail"].read_text()


def test_zero_hits_is_valid_but_missing_header_is_not(kofam_job):
    kofam_job["detail"].write_text(HEADER)
    run(**kofam_job["args"])
    assert read_tsv(kofam_job["args"]["output_dir"] / "genes.tsv")[0]["assignment_status"] == "unannotated"
    kofam_job["detail"].write_text("")
    # A changed input invalidates the previously successful run.
    kofam_job["args"]["protein"].write_text(">plant_g1\nMKK\n")
    with pytest.raises(ValueError, match="missing.*header"):
        run(**kofam_job["args"])


@pytest.mark.parametrize("row,error", [
    (detail_row("foreign_gene"), "does not belong"),
    (detail_row("plant_g1", "K99999"), "absent from reference"),
    (detail_row("plant_g1", score="NaN"), "invalid score"),
    (detail_row("plant_g1", evalue="inf"), "invalid E-value"),
    (detail_row("plant_g1", evalue="-1"), "invalid E-value"),
    (detail_row("plant_g1", threshold="NaN"), "invalid threshold"),
    (detail_row("plant_g1", threshold="9"), "threshold differs"),
    (detail_row("plant_g1", score="2"), "marker disagrees"),
    (detail_row("plant_g1", score="20", marker=""), "marker disagrees"),
    (detail_row("plant_g1", "K00003", "", "20"), "has no threshold"),
    (detail_row("plant_g1") + detail_row("plant_g1", score="12"), "conflicting duplicate"),
])
def test_invalid_tool_output_never_publishes(kofam_job, row, error):
    kofam_job["detail"].write_text(HEADER + row)
    with pytest.raises(ValueError, match=error):
        run(**kofam_job["args"])
    assert not kofam_job["args"]["output_dir"].exists()
    status = json.loads((kofam_job["args"]["work_dir"] / "status.json").read_text())
    assert status["state"] == "failed"


def test_pre_rounding_acceptance_is_preserved(tmp_path):
    path = tmp_path / "detail.tsv"
    path.write_text(HEADER + detail_row("g1", threshold="10.02", score="10.0")
                    + detail_row("g2", threshold="10.02", score="10.0", marker=""))
    proteins = [{"gene_id": gene, "terminal_stop_stripped": 0} for gene in ["g1", "g2"]]
    _, genes = parse_detail(path, "plant", proteins, {"K00001": {"threshold": "10.024"}})
    assert [gene["assignment_status"] for gene in genes] == ["unique", "below_threshold"]


@pytest.mark.parametrize("content,error", [
    (">g1\nMK*K\n", "internal stop"),
    (">g1\nMK**\n", "internal stop"),
    (">g1\n*\n", "empty or invalid"),
    (">g1\nMK\n>g1\nMKK\n", "duplicate"),
    ("> g1\nMK\n", "header"),
    ("MK\n>g1\nMK\n", "precedes"),
])
def test_protein_validation_is_explicit(tmp_path, content, error):
    source = tmp_path / "input.fa"
    source.write_text(content)
    with pytest.raises(ValueError, match=error):
        normalize_protein(source, tmp_path / "normalized.faa")


def test_resume_verifies_outputs_and_input_tool_options(kofam_job):
    args = kofam_job["args"]
    first = run(**args)
    assert run(**args)["fingerprint"] == first["fingerprint"]
    assert len(kofam_job["events"].read_text().splitlines()) == 1
    assert json.loads((args["work_dir"] / "status.json").read_text())["reused"]
    # Even a same-size corruption must force annotation again.
    genes = args["output_dir"] / "genes.tsv"
    genes.write_text(genes.read_text().replace("K00001", "K00002"))
    run(**args)
    assert read_tsv(genes)[0]["selected_ko"] == "K00001"
    assert len(kofam_job["events"].read_text().splitlines()) == 2
    args["protein"].write_text(args["protein"].read_text().replace("MK", "MQ"))
    second = run(**args)
    assert second["fingerprint"] != first["fingerprint"]
    command = Path(args["command"])
    command.write_text(command.read_text() + "\n# tool revision\n")
    third = run(**args)
    assert third["fingerprint"] != second["fingerprint"]
    args["threads"] = 2
    fourth = run(**args)
    assert fourth["fingerprint"] != third["fingerprint"]
    reference = json.loads(args["reference"].read_text())
    reference["release"] = "fixture-revised"
    write_json(args["reference"], reference)
    assert run(**args)["fingerprint"] != fourth["fingerprint"]


def test_failed_rerun_keeps_previous_complete_output_and_retries_fresh(kofam_job, monkeypatch):
    args = kofam_job["args"]
    first = run(**args)
    args["protein"].write_text(">plant_g1\nMKQ*\n")
    monkeypatch.setenv("FAKE_KOFAM_FAIL", "1")
    with pytest.raises(subprocess.CalledProcessError):
        run(**args)
    saved = json.loads((args["output_dir"] / "provenance.json").read_text())
    assert saved["fingerprint"] == first["fingerprint"]
    failed = json.loads((args["work_dir"] / "status.json").read_text())
    assert failed["state"] == "failed"
    assert "intentional fake failure" in (Path(failed["work"]) / "result/stderr.log").read_text()
    monkeypatch.delenv("FAKE_KOFAM_FAIL")
    assert run(**args)["fingerprint"] != first["fingerprint"]
    events = [json.loads(line) for line in kofam_job["events"].read_text().splitlines()]
    assert len(events) == 3
    assert events[1]["args"][-1] != events[2]["args"][-1]  # partial HMMER work is not trusted
