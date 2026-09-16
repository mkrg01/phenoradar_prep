"""PhenoRadar TPM conversion must preserve values and reject implicit pooling."""
import pytest

from common import read_tsv, write_tsv
from export_species_tpm import export


@pytest.fixture
def tables(tmp_path):
    samples, source, output = [tmp_path / name for name in
                               ["samples.tsv", "tpm.tsv", "tpm_species.tsv"]]
    write_tsv(samples, ["species", "run"], [dict(species="Plant_A", run="A1"),
                                           dict(species="Plant_B", run="B1")])
    write_tsv(source, ["species", "run", "orthogroup", "tpm"], [
        dict(species="Plant_A", run="A1", orthogroup="OG1", tpm="0.0000012300"),
        dict(species="Plant_A", run="A1", orthogroup="OG2", tpm="0.00000"),
        dict(species="Plant_B", run="B1", orthogroup="OG3", tpm="1e+06"),
    ])
    return samples, source, output


def test_exact_three_columns_values_and_sparse_coordinates(tables):
    samples, source, output = tables
    before = source.read_bytes()
    export(*tables)
    assert output.read_text() == (
        "species\torthogroup\ttpm\n"
        "Plant_A\tOG1\t0.0000012300\n"
        "Plant_A\tOG2\t0.00000\n"
        "Plant_B\tOG3\t1e+06\n")
    assert source.read_bytes() == before
    export(*tables)  # Replacing a previous completed export is supported.
    assert len(read_tsv(output)) == 3


def test_multiple_runs_are_rejected_before_replacing_output(tables):
    samples, source, output = tables
    rows = read_tsv(samples)
    rows.append(dict(species="Plant_A", run="A2"))
    write_tsv(samples, ["species", "run"], rows)
    output.write_text("previous result\n")
    with pytest.raises(ValueError, match="multiple runs per species.*Plant_A"):
        export(*tables)
    assert output.read_text() == "previous result\n"


@pytest.mark.parametrize("problem", ["duplicate_og", "repeated_run", "missing_run", "wrong_species",
                                     "blank_og", "nan", "inf", "-1", "", "nonnumeric"])
def test_invalid_tables_do_not_replace_previous_export(tables, problem):
    samples, source, output = tables
    rows = read_tsv(source)
    if problem == "duplicate_og":
        rows[1]["orthogroup"] = rows[0]["orthogroup"]
    elif problem == "repeated_run":
        rows.append({**rows[0], "orthogroup": "OG4"})
    elif problem == "missing_run":
        rows.pop()
    elif problem == "wrong_species":
        rows[0]["species"] = "Plant_B"
    elif problem == "blank_og":
        rows[0]["orthogroup"] = " "
    else:
        rows[-1]["tpm"] = problem
    write_tsv(source, list(rows[0]), rows)
    output.write_text("previous result\n")
    with pytest.raises(ValueError):
        export(*tables)
    assert output.read_text() == "previous result\n"


def test_duplicate_manifest_run_is_rejected(tables):
    samples, source, output = tables
    rows = read_tsv(samples)
    write_tsv(samples, list(rows[0]), rows + [rows[0]])
    with pytest.raises(ValueError, match="duplicate run"):
        export(*tables)


def test_input_cannot_be_overwritten(tables):
    samples, source, output = tables
    before = source.read_bytes()
    with pytest.raises(ValueError, match="different files"):
        export(samples, source, source)
    assert source.read_bytes() == before
