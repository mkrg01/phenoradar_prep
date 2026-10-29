"""Per-copy abundances remain exact, bounded-memory, and independent of OG selection."""
import gc
import tracemalloc

import pytest

from common import read_tsv, sha256, write_tsv
from sequence_expression import FIELDS, merge


def inputs(tmp_path, compressed=False):
    rows = []
    for sample, run in [("Plant_A_R_1", "R_1"), ("Plant_A_R_2", "R_2")]:
        path = tmp_path / (run + (".tsv.gz" if compressed else ".tsv"))
        write_tsv(path, ["target_id", "length", "tpm"], [
            {"target_id": sample + "_g1", "length": "123", "tpm": "80.000000000000001"},
            {"target_id": sample + "_g2", "length": "99", "tpm": "2e1"},
            {"target_id": sample + "_g9", "length": "100", "tpm": "0.00000"},
        ])
        rows.append(dict(species=sample, run=run, abundance=str(path)))
    samples = tmp_path / "samples.tsv"
    write_tsv(samples, list(rows[0]), reversed(rows))
    return samples, rows


@pytest.mark.parametrize("compressed", [False, True])
def test_exact_original_values_all_copies_zero_and_unmapped_targets(tmp_path, compressed):
    samples, manifest = inputs(tmp_path, compressed)
    before = {r["abundance"]: sha256(r["abundance"]) for r in manifest}
    output = tmp_path / "expression/sequence_tpm.tsv.gz"
    assert merge(samples, output) == {"samples": 2, "sequences": 6}
    result = read_tsv(output)
    assert list(result[0]) == FIELDS
    assert result == [dict(species=s["species"], sequence_id=r["target_id"], tpm=r["tpm"])
                      for s in manifest for r in read_tsv(s["abundance"])]
    assert {r["abundance"]: sha256(r["abundance"]) for r in manifest} == before
    first = output.read_bytes()
    merge(samples, output)
    assert output.read_bytes() == first  # Reproducible gzip, including its header.
    assert not any(r["sequence_id"].endswith("_g3") for r in result)  # Missing is not zero.


@pytest.mark.parametrize("problem", ["duplicate", "empty_id", "foreign", "empty", "header",
    "duplicate_header", "short_row", "extra_column", "negative", "nan", "inf", "nonnumeric"])
def test_bad_later_sample_keeps_previous_output(tmp_path, problem):
    samples, manifest = inputs(tmp_path)
    path = manifest[-1]["abundance"]
    rows = read_tsv(path)
    fields = list(rows[0])
    if problem == "duplicate": rows.append(dict(rows[0]))
    elif problem == "empty_id": rows[0]["target_id"] = ""
    elif problem == "foreign": rows[0]["target_id"] = manifest[0]["species"] + "_g1"
    elif problem == "empty": rows = []
    elif problem == "header": fields = ["Name", "TPM"]; rows = []
    elif problem == "duplicate_header": fields = ["target_id", "tpm", "tpm"]; rows = []
    elif problem not in {"short_row", "extra_column"}:
        rows[0]["tpm"] = {"negative": "-1", "nan": "NaN", "inf": "inf", "nonnumeric": "NA"}[problem]
    write_tsv(path, fields, rows)
    if problem in {"short_row", "extra_column"}:
        with open(path, "a") as handle:
            handle.write("broken\t1\n" if problem == "short_row" else "broken\t1\t2\t3\n")
    output = tmp_path / "old.tsv.gz"
    output.write_bytes(b"previous completed output")
    with pytest.raises(ValueError): merge(samples, output)
    assert output.read_bytes() == b"previous completed output"
    assert not list(tmp_path.glob(".old.tsv.gz.*"))


@pytest.mark.parametrize("problem", ["duplicate_sample", "duplicate_run", "empty"])
def test_manifest_requires_unambiguous_sample_keys(tmp_path, problem):
    samples, rows = inputs(tmp_path)
    if problem == "duplicate_sample": rows[1]["species"] = rows[0]["species"]
    elif problem == "duplicate_run": rows[1]["run"] = rows[0]["run"]
    else: rows = []
    write_tsv(samples, ["species", "run", "abundance"], rows)
    output = tmp_path / "out.tsv.gz"
    with pytest.raises(ValueError): merge(samples, output)
    assert not output.exists()


def test_zero_expression_is_preserved_without_imputing_composition(tmp_path):
    samples, rows = inputs(tmp_path)
    for sample in rows:
        write_tsv(sample["abundance"], ["target_id", "tpm"],
                  [dict(target_id=sample["species"] + "_g1", tpm="0")])
    output = tmp_path / "out.tsv.gz"
    merge(samples, output)
    assert [r["tpm"] for r in read_tsv(output)] == ["0", "0"]


def test_memory_does_not_grow_with_all_sample_sequence_ids(tmp_path):
    samples = tmp_path / "samples.tsv"
    def peak(count):
        rows = []
        for i in range(count):
            sample = f"Plant_R{i}"
            path = tmp_path / f"{sample}.tsv.gz"
            write_tsv(path, ["target_id", "tpm"],
                      (dict(target_id=f"{sample}_g{j}", tpm="0.00000001") for j in range(4000)))
            rows.append(dict(species=sample, run=f"R{i}", abundance=str(path)))
        write_tsv(samples, list(rows[0]), rows)
        gc.collect()
        tracemalloc.start()
        try:
            merge(samples, tmp_path / "out.tsv.gz")
            return tracemalloc.get_traced_memory()[1]
        finally:
            tracemalloc.stop()
    assert peak(15) < peak(1) + 2_000_000
