"""Expression merging preserves its tables without retaining all row dictionaries."""
import gc
import tracemalloc
from pathlib import Path

import pytest

import merge_tpm
from common import read_tsv, write_json, write_tsv


TABLE_FIELDS = ["species", "run", "orthogroup", "tpm", "tpm_sum"]
QC_FIELDS = ["species", "run", "targets", "protein_genes", "quantified_proteins", "mapped_targets",
             "ambiguous_targets", "retained_targets", "orthogroups", "total_tpm", "mapped_tpm",
             "mapped_tpm_fraction", "retained_tpm", "retained_tpm_fraction", "multimap"]
OUTPUTS = ["tpm.tsv", "tpm_wide.tsv", "tpm_sum.tsv", "tpm_sum_wide.tsv", "mapping_qc.tsv"]


def save_run(directory, species, run, entries, compressed=False):
    rows = [dict(species=species, run=run, orthogroup=group, tpm=tpm, tpm_sum=summed)
            for group, tpm, summed in entries]
    path = directory / f"{run}.tsv"
    if compressed:
        path = Path(str(path) + ".gz")
    write_tsv(path, TABLE_FIELDS, rows)
    report = {field: index for index, field in enumerate(QC_FIELDS)}
    report.update(species=species, run=run, multimap="error")
    write_json(directory / f"{run}.qc.json", report)
    return rows, report, path


def test_exact_tables_single_input_pass_and_independent_runs(tmp_path, monkeypatch):
    samples = tmp_path / "samples.tsv"
    manifest = [{"species": "Beta", "run": "B1"}, {"species": "Alpha", "run": "A2"},
                {"species": "Alpha", "run": "A1"}]
    write_tsv(samples, ["species", "run"], manifest)
    runs = tmp_path / "runs"
    odd_group = 'OG4\t"雪\n'
    data = [
        ("Alpha", "A1", [("OG2", "1.2300", "0.0010"), ("OG1", "1e-12", "02.500")], True),
        ("Alpha", "A2", [("OG3", "0", "-0.0")], True),
        ("Beta", "B1", [("OG1", "-0", "9.999999999999999"), (odd_group, "2.0", "3e+10")], False),
    ]
    all_rows, reports, paths = [], [], []
    for species, run, entries, compressed in data:
        rows, report, path = save_run(runs, species, run, entries, compressed)
        all_rows.extend(rows)
        reports.append(report)
        paths.append(path)
    # A published gzip product takes precedence over a stale plain table.
    write_tsv(runs / "A1.tsv", TABLE_FIELDS, [])
    reads = []
    original_read = merge_tpm.read_tsv

    def tracked_read(path):
        reads.append(Path(path))
        return original_read(path)

    monkeypatch.setattr(merge_tpm, "read_tsv", tracked_read)
    actual = tmp_path / "actual"
    merge_tpm.merge(samples, runs, actual)
    assert reads == [samples, *paths]

    expected = tmp_path / "expected"
    ordered = sorted(manifest, key=lambda row: (row["species"], row["run"]))
    groups = ["OG1", "OG2", "OG3", odd_group]
    for value in ["tpm", "tpm_sum"]:
        fields = ["species", "run", "orthogroup", value]
        write_tsv(expected / f"{value}.tsv", fields,
                  [{key: row[key] for key in fields} for row in all_rows])
        wide = []
        for sample in ordered:
            row = dict(sample)
            row.update(dict.fromkeys(groups, 0))
            row.update({entry["orthogroup"]: entry[value] for entry in all_rows
                        if entry["run"] == sample["run"]})
            wide.append(row)
        write_tsv(expected / f"{value}_wide.tsv", ["species", "run", *groups], wide)
    write_tsv(expected / "mapping_qc.tsv", QC_FIELDS, reports)
    for name in OUTPUTS:
        assert (actual / name).read_bytes() == (expected / name).read_bytes(), name


@pytest.mark.parametrize("problem", ["empty", "species", "run", "duplicate", "qc"])
def test_bad_later_run_preserves_previous_outputs(tmp_path, problem):
    samples = tmp_path / "samples.tsv"
    write_tsv(samples, ["species", "run"], [{"species": "Alpha", "run": "A1"},
                                           {"species": "Beta", "run": "B1"}])
    runs = tmp_path / "runs"
    save_run(runs, "Alpha", "A1", [("OG1", "1.0", "2.0")])
    _, _, path = save_run(runs, "Beta", "B1", [("OG1", "3.0", "4.0")])
    rows = read_tsv(path)
    if problem == "empty":
        rows = []
    elif problem in {"species", "run"}:
        rows[0][problem] = "wrong"
    elif problem == "duplicate":
        rows.append(dict(rows[0]))
    else:
        write_json(runs / "B1.qc.json", {"species": "Beta", "run": "B1"})
    write_tsv(path, TABLE_FIELDS, rows)
    out = tmp_path / "out"
    out.mkdir()
    for name in OUTPUTS:
        (out / name).write_text("previous output\n")
    with pytest.raises(KeyError if problem == "qc" else ValueError):
        merge_tpm.merge(samples, runs, out)
    assert {p.name for p in out.iterdir()} == set(OUTPUTS)
    assert all((out / name).read_text() == "previous output\n" for name in OUTPUTS)


def test_many_runs_fit_compact_memory_budget(tmp_path):
    samples = tmp_path / "samples.tsv"
    manifest = [{"species": f"species_{i:03d}", "run": f"run_{i:03d}"} for i in range(96)]
    write_tsv(samples, ["species", "run"], manifest)
    runs = tmp_path / "runs"
    for sample in manifest:
        save_run(runs, sample["species"], sample["run"],
                 [(f"OG{i:05d}", f"{i}.0123456789", f"{i}.9876543210") for i in range(1200)])
    gc.collect()
    tracemalloc.start()
    try:
        merge_tpm.merge(samples, runs, tmp_path / "out")
        _, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    # About 115,000 rows: keeping their full dictionaries needs tens of MiB.
    assert peak < 12 * 1024 * 1024, f"peak allocation was {peak / 1024**2:.1f} MiB"


def test_empty_selection_writes_headers(tmp_path):
    samples = tmp_path / "samples.tsv"
    write_tsv(samples, ["species", "run"], [])
    out = tmp_path / "out"
    merge_tpm.merge(samples, tmp_path / "runs", out)
    assert (out / "tpm_wide.tsv").read_text() == "species\trun\n"
    assert all(read_tsv(out / name) == [] for name in OUTPUTS)
