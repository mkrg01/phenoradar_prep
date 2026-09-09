"""Biological mapping failures, persistent retrieval, and dating integration."""
import json
import os
from pathlib import Path
import sqlite3

import pytest
from ete4 import Tree

from common import read_tsv, write_tsv
from date_phylogeny import date
from infer_phylogeny import read_tree
from timetree_calibrations import (cached_response, interpret_response, nwkit_backend, prepare,
                                   representative_plan, select_calibrations, species_taxids)


def payload(ids, mrca="root", low=90, high=110):
    return {"found_ids": list(map(str, ids)), "missing_ids": [], "mrca_id": mrca,
            "precomputed_age": (low + high) / 2, "precomputed_ci_low": low, "precomputed_ci_high": high,
            "study_data": [{"title": f"Synthetic study {i}", "pub_year": 2000 + i,
                            "first_author": "Test", "tree_count": 1} for i in range(5)]}


def fake_fetch(body):
    return lambda url: {"url": url, "status_code": 200, "text": json.dumps(body), "error": None, "elapsed": 0}


def query():
    return {"node": "N1", "clade_taxa": 4, "query_taxa": ["A", "B", "C", "D"],
            "children": [["A", "B"], ["C", "D"]]}


def test_missing_taxa_must_leave_both_child_lineages_represented():
    mapping = dict(zip("ABCD", range(1, 5)))
    data = payload([1, 2])
    data["missing_ids"] = [3, 4]
    assert interpret_response(data, query(), mapping, 5)["reason"] == "missing_child_lineage"
    data = payload([1, 3])
    data["missing_ids"] = [2, 4]
    row = interpret_response(data, query(), mapping, 5)
    assert row["status"] == "candidate"
    assert row["used_taxa"] == ["A", "C"]


@pytest.mark.parametrize("change,reason", [
    ({"precomputed_ci_low": 0}, "invalid_or_point_only_bounds"),
    ({"precomputed_ci_low": 120}, "invalid_or_point_only_bounds"),
    ({"precomputed_age": float("nan")}, "invalid_or_point_only_bounds"),
    ({"precomputed_ci_low": 100, "precomputed_ci_high": 100}, "invalid_or_point_only_bounds"),
    ({"study_data": [{"title": "One study", "tree_count": 50}]}, "insufficient_studies"),
    ({"found_ids": [1, 2, 3, 999]}, "unexpected_taxon_response"),
    ({"mrca_id": None}, "missing_mrca_identifier"),
])
def test_invalid_calibration_candidates(change, reason):
    data = payload(range(1, 5))
    data.update(change)
    assert interpret_response(data, query(), dict(zip("ABCD", range(1, 5))), 5)["reason"] == reason


def test_cache_reuses_raw_response_offline_and_rejects_tampering(tmp_path):
    data = payload([1, 2])
    body, record, path = cached_response([2, 1], tmp_path, delay=0, backend=(fake_fetch(data), {"synthetic": True}))
    assert body == data
    assert cached_response([1, 2], tmp_path, offline=True)[1] == record
    with pytest.raises(ValueError, match="offline.*miss"):
        cached_response([1, 3], tmp_path, offline=True)
    record["text"] = "{}"
    path.write_text(json.dumps(record))
    with pytest.raises(ValueError, match="invalid.*cache"):
        cached_response([1, 2], tmp_path, offline=True)


def test_transport_failure_is_not_cached_as_missing_biology(tmp_path):
    fetch = lambda url: {"url": url, "status_code": 503, "text": "unavailable", "error": None, "elapsed": 0}
    with pytest.raises(ValueError, match="transport failed"):
        cached_response([1, 2], tmp_path, delay=0, backend=(fetch, {}))
    assert len(list(tmp_path.glob("*.failure.json"))) == 1
    with pytest.raises(ValueError, match="offline.*miss"):
        cached_response([1, 2], tmp_path, offline=True)


def test_real_nwkit_client_with_recorded_http_response(tmp_path, monkeypatch):
    pytest.importorskip("nwkit")
    import requests
    fetch, tool = nwkit_backend()
    data = payload([1, 2])
    calls = []
    def get(url, **kwargs):
        calls.append(url)
        assert kwargs == {"timeout": 30, "stream": True}
        response = requests.Response()
        response.status_code = 200
        response._content = json.dumps(data).encode()
        response._content_consumed = True
        return response
    monkeypatch.setattr(requests, "get", get)
    result = cached_response([1, 2], tmp_path, delay=0, backend=(fetch, tool))
    assert result[0] == data
    cached_response([2, 1], tmp_path, offline=True)
    assert len(calls) == 1
    assert result[1]["nwkit"]["version"] == "0.43.12"


def test_ambiguous_mrca_ids_and_conflicting_bounds_stop_automatic_dating(tmp_path):
    tree = Tree("(A:1,(B:1,(C:1,D:1)N3:1)N2:1)N1;", parser=1)
    mapping = dict(zip("ABCD", range(1, 5)))
    rows = []
    for taxa, children, node in [(["A", "B"], [["A"], ["B"]], "N1"),
                                 (["C", "D"], [["C"], ["D"]], "N3")]:
        row = interpret_response(payload([mapping[n] for n in taxa]),
              {"node": node, "query_taxa": taxa, "children": children}, mapping, 5)
        row.update(url="synthetic", retrieved_at="test")
        rows.append(row)
    output = tmp_path / "calibrations.tsv"
    assert select_calibrations(tree, rows, output)["status"] == "no_valid_calibrations"
    assert all(r["reason"] == "ambiguous_repeated_timetree_mrca" for r in rows)
    for row, mrca in zip(rows, ["root", "child"]):
        row.update(mrca_id=mrca, status="candidate", reason="")
    rows[1].update(min_age_ma=120, max_age_ma=150)
    assert select_calibrations(tree, rows, output)["status"] == "conflicting_calibrations"
    assert read_tsv(output) == []


def test_representatives_preserve_root_and_are_independent_of_child_order():
    selected_sets = []
    for text in ["(A:1,((B:1,C:1):1,(D:1,(E:1,F:1):1):1):1);",
                 "((((F:1,E:1):1,D:1):1,(C:1,B:1):1):1,A:1);"]:
        tree = Tree(text, parser=1)
        mapping = dict(zip("ABCDEF", range(1, 7)))
        selected, sizes, members = representative_plan(tree, mapping, {"C": 100, "F": 100}, 4)
        selected_sets.append(set(selected))
        assert len(selected) == 4
        assert all(members[c] for c in tree.children)
        assert sizes[tree] == 6
    assert selected_sets[0] == selected_sets[1]


def inputs(tmp_path):
    tree = tmp_path / "tree.nwk"
    tree.write_text("(A:0.1,(B:0.08,(C:0.04,D:0.04):0.04):0.02);\n")
    metadata, database, coverage = tmp_path / "metadata.tsv", tmp_path / "taxa.sqlite", tmp_path / "coverage.tsv"
    write_tsv(metadata, ["species", "taxid"], [{"species": n, "taxid": i} for i, n in enumerate("ABCD", 1)])
    write_tsv(coverage, ["species", "gene_trees"], [{"species": n, "gene_trees": 10} for n in "ABCD"])
    with sqlite3.connect(database) as db:
        db.executescript("CREATE TABLE species (taxid INTEGER, rank TEXT, track TEXT, spname TEXT);"
                        "CREATE TABLE merged (taxid_old INTEGER, taxid_new INTEGER);")
        db.executemany("INSERT INTO species VALUES (?, 'species', ?, ?)", [(i, str(i), n) for i, n in enumerate("ABCD", 1)])
    return tree, metadata, database, coverage


def test_taxid_resolution_rejects_aliases_mapping_to_same_species(tmp_path):
    _, metadata, database, _ = inputs(tmp_path)
    write_tsv(metadata, ["species", "taxid"], [{"species": "A", "taxid": 11}, {"species": "Alias", "taxid": 1}])
    with sqlite3.connect(database) as db:
        db.execute("INSERT INTO merged VALUES (11, 1)")
    mapping, report = species_taxids(metadata, database, {"A", "Alias"})
    assert not mapping
    assert {r["status"] for r in report} == {"duplicate_species_taxid"}


@pytest.mark.parametrize("maximum,expected_sizes", [(1, [4]), (2, [4, 3]), (32, [4, 3, 2])])
def test_small_clades_are_eligible_and_query_budget_still_applies(tmp_path, maximum, expected_sizes):
    tree, metadata, database, coverage = inputs(tmp_path)
    cache, out = tmp_path / "cache", tmp_path / "calibrations"
    # All nodes have fewer than 20 species. Distinct cached MRCA IDs and
    # consistent ages let us test selection through final calibration output.
    for ids, mrca, low, high in [(range(1, 5), "root", 90, 110),
                                (range(2, 5), "child", 60, 80), ([3, 4], "pair", 30, 50)]:
        cached_response(ids, cache, delay=0,
                        backend=(fake_fetch(payload(ids, mrca, low, high)), {"synthetic": True}))
    settings = {"max_representatives": 64, "max_queries": maximum,
                "min_studies": 5, "offline": True, "request_delay_seconds": 0}
    prepare(tree, metadata, database, coverage, out, cache, settings)
    candidates = read_tsv(out / "candidates.tsv")
    assert [int(row["clade_taxa"]) for row in candidates] == expected_sizes
    assert all(row["status"] == "retained" for row in candidates)
    assert len(read_tsv(out / "calibrations.tsv")) == len(expected_sizes)
    provenance = json.loads((out / "provenance.json").read_text())
    assert provenance["status"] == "ready"
    assert provenance["eligible_nodes"] == 3
    assert provenance["queries"] == len(expected_sizes)


def test_preparation_and_real_treepl_keep_topology_and_substitution_tree(tmp_path):
    tree, metadata, database, coverage = inputs(tmp_path)
    original = tree.read_bytes()
    settings = {"max_representatives": 4, "max_queries": 1,
                "min_studies": 5, "offline": True, "request_delay_seconds": 0}
    cache = tmp_path / "cache"
    cached_response(range(1, 5), cache, delay=0, backend=(fake_fetch(payload(range(1, 5))), {"synthetic": True}))
    out = tmp_path / "calibrations"
    prepare(tree, metadata, database, coverage, out, cache, settings)
    assert json.loads((out / "provenance.json").read_text())["status"] == "ready"
    assert len(read_tsv(out / "calibrations.tsv")) == 1
    treepl = os.environ.get("TREEPL_BIN")
    if treepl:
        provenance = tmp_path / "tree.json"
        provenance.write_text(json.dumps({"branch_length_unit": "substitutions_per_site", "outgroup": "A", "mean_gene_length": 250, "total_gene_sites": 750}))
        date(tree, provenance, out / "calibrations.tsv", tmp_path / "dated", treepl, settings={"smoothing": 10, "replicates": 2})
        dated = read_tree(tmp_path / "dated/species_tree.dated.nwk", "ABCD")
        assert 90 - 1e-3 <= dated.get_distance("A", "B") / 2 <= 110 + 1e-3
    assert tree.read_bytes() == original
