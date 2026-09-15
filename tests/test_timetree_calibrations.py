"""Biological mapping failures, persistent retrieval, and dating integration."""
import json
import os
from pathlib import Path
import sqlite3

import pytest
from ete4 import Tree

from common import file_record, read_tsv, write_tsv
from date_phylogeny import calibration_rows, date
from infer_phylogeny import read_tree
from timetree_calibrations import (cached_response, interpret_response, nwkit_backend, prepare,
                                   node_queries, select_calibrations, species_taxids)


@pytest.fixture(autouse=True)
def no_live_requests(monkeypatch):
    import requests
    monkeypatch.setattr(requests, "get", lambda *a, **k: pytest.fail("unexpected live TimeTree request"))
    monkeypatch.setattr("timetree_calibrations.time.sleep", lambda delay: None)


def payload(ids, mrca="root", low=90, high=110):
    return {"found_ids": list(map(str, ids)), "missing_ids": [], "mrca_id": mrca,
            "precomputed_age": (low + high) / 2, "precomputed_ci_low": low, "precomputed_ci_high": high,
            "study_data": [{"title": f"Synthetic study {i}", "pub_year": 2000 + i,
                            "first_author": "Test", "tree_count": 1} for i in range(5)]}


def fake_fetch(body):
    return lambda url: {"url": url, "status_code": 200, "text": json.dumps(body), "error": None, "elapsed": 0}


def cache_missing_estimates(tree_path, mapping, cache):
    """Fill uncached internal-node responses with synthetic insufficient studies.

    Tests preload their valid calibrations first. This covers the remaining
    clades without requiring live HTTP in Snakemake subprocesses.
    """
    for node in read_tree(tree_path).traverse():
        if not node.is_leaf:
            ids = sorted(mapping[s] for s in node.leaf_names())
            data = payload(ids, mrca="+".join(map(str, ids)))
            data["study_data"] = []
            cached_response(ids, cache, delay=0, backend=(fake_fetch(data), {"synthetic": True}))


def query():
    return {"node": "N1", "clade_taxa": 4, "query_taxa": ["A", "B", "C", "D"],
            "children": [["A", "B"], ["C", "D"]]}


def test_missing_taxa_must_leave_both_child_lineages_represented():
    mapping = dict(zip("ABCD", range(1, 5)))
    data = payload([1, 2])
    data["missing_ids"] = [3, 4]
    assert interpret_response(data, query(), mapping)["reason"] == "missing_child_lineage"
    data = payload([1, 3])
    data["missing_ids"] = [2, 4]
    row = interpret_response(data, query(), mapping)
    assert row["status"] == "candidate"
    assert row["used_taxa"] == ["A", "C"]
    assert row["missing_taxids"] == [2, 4]


@pytest.mark.parametrize("change,reason", [
    ({"precomputed_ci_low": 0}, "invalid_or_point_only_bounds"),
    ({"precomputed_ci_low": 120}, "invalid_or_point_only_bounds"),
    ({"precomputed_age": float("nan")}, "invalid_or_point_only_bounds"),
    ({"precomputed_ci_low": 100, "precomputed_ci_high": 100}, "invalid_or_point_only_bounds"),
    ({"study_data": [{"title": "One study", "tree_count": 50}]}, "insufficient_studies"),
    ({"found_ids": [1, 2, 3, 999]}, "unexpected_taxon_response"),
    ({"found_ids": [1, 3]}, "unexpected_taxon_response"),
    ({"missing_ids": [1]}, "unexpected_taxon_response"),
    ({"study_data": {}}, "unexpected_study_response"),
    ({"mrca_id": None}, "missing_mrca_identifier"),
])
def test_invalid_calibration_candidates(change, reason):
    data = payload(range(1, 5))
    data.update(change)
    assert interpret_response(data, query(), dict(zip("ABCD", range(1, 5))))["reason"] == reason


def test_five_studies_are_fixed_and_duplicate_trees_do_not_inflate_count():
    data = payload(range(1, 5))
    mapping = dict(zip("ABCD", range(1, 5)))
    data["study_data"].pop()
    data["study_data"].append(dict(data["study_data"][0], title="  SYNTHETIC   STUDY 0 ", tree_count=50))
    row = interpret_response(data, query(), mapping)
    assert row["reason"] == "insufficient_studies" and row["studies"] == 4
    assert row["min_age_ma"] == 90 and row["max_age_ma"] == 110
    assert len(row["study_records"]) == 4
    data["study_data"].append({"title": "Fifth study", "pub_year": 2020, "first_author": "Test"})
    assert interpret_response(data, query(), mapping)["status"] == "candidate"


def test_cache_reuses_raw_response_fetches_only_missing_and_rejects_tampering(tmp_path, monkeypatch):
    data = payload([1, 2])
    calls = []
    def fetch(url):
        calls.append(url)
        return fake_fetch(data)(url)
    monkeypatch.setattr("timetree_calibrations.nwkit_backend", lambda: (fetch, {"synthetic": True}))
    body, record, path = cached_response([2, 1], tmp_path, delay=0)
    assert body == data
    original, timestamp = path.read_bytes(), path.stat().st_mtime_ns
    assert cached_response([1, 2], tmp_path)[1] == record
    assert calls == [record["url"]]
    data = payload([1, 3])
    body, new_record, new_path = cached_response([1, 3], tmp_path, delay=0)
    assert body == data and new_path != path
    assert calls == [record["url"], new_record["url"]]
    assert path.read_bytes() == original and path.stat().st_mtime_ns == timestamp
    record["text"] = "{}"
    path.write_text(json.dumps(record))
    with pytest.raises(ValueError, match="invalid.*cache"):
        cached_response([1, 2], tmp_path)
    assert len(calls) == 2


def test_transport_failure_is_not_cached_as_missing_biology(tmp_path):
    fetch = lambda url: {"url": url, "status_code": 503, "text": "unavailable", "error": None, "elapsed": 0}
    with pytest.raises(ValueError, match="transport failed"):
        cached_response([1, 2], tmp_path, delay=0, backend=(fetch, {}))
    assert len(list(tmp_path.glob("*.failure.json"))) == 1
    # A failed request must be retried and then become reusable after success.
    data = payload([1, 2])
    body, record, path = cached_response([1, 2], tmp_path, delay=0, backend=(fake_fetch(data), {}))
    assert body == data and path.is_file()
    assert cached_response([1, 2], tmp_path)[1] == record


def test_real_nwkit_client_with_recorded_http_response(tmp_path, monkeypatch):
    pytest.importorskip("nwkit")
    import requests
    fetch, tool = nwkit_backend()
    data = payload([1, 2])
    calls = []
    def get(url, **kwargs):
        calls.append(url)
        assert kwargs == ({"timeout": 30} if tool["version"] == "0.27.0" else {"timeout": 30, "stream": True})
        response = requests.Response()
        response.status_code = 200
        response._content = json.dumps(data).encode()
        response._content_consumed = True
        return response
    monkeypatch.setattr(requests, "get", get)
    result = cached_response([1, 2], tmp_path, delay=0, backend=(fetch, tool))
    assert result[0] == data
    cached_response([2, 1], tmp_path)
    assert len(calls) == 1
    assert result[1]["nwkit"]["version"] == tool["version"]


def test_ambiguous_mrca_ids_and_conflicting_bounds_stop_automatic_dating(tmp_path):
    tree = Tree("(A:1,(B:1,(C:1,D:1)N3:1)N2:1)N1;", parser=1)
    mapping = dict(zip("ABCD", range(1, 5)))
    rows = []
    for taxa, children, node in [(["A", "B"], [["A"], ["B"]], "N1"),
                                 (["C", "D"], [["C"], ["D"]], "N3")]:
        row = interpret_response(payload([mapping[n] for n in taxa]),
              {"node": node, "query_taxa": taxa, "children": children}, mapping)
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


def test_all_nodes_and_labels_are_independent_of_child_order():
    plans = []
    for text in ["(A:1,((B:1,C:1):1,(D:1,(E:1,F:1):1):1):1);",
                 "((((F:1,E:1):1,D:1):1,(C:1,B:1):1):1,A:1);"]:
        tree = Tree(text, parser=1)
        mapping = dict(zip("ABCDEF", range(1, 7)))
        queries = [query for node, query in node_queries(tree, mapping)]
        plans.append(queries)
        assert len(queries) == 5
        assert queries[0]["query_taxa"] == list("ABCDEF")
        assert {tuple(q["query_taxa"]) for q in queries} == {
            tuple("ABCDEF"), tuple("BCDEF"), tuple("BC"), tuple("DEF"), tuple("EF")}
    assert plans[0] == plans[1]


def inputs(tmp_path):
    tree = tmp_path / "tree.nwk"
    tree.write_text("(A:0.1,(B:0.08,(C:0.04,D:0.04):0.04):0.02);\n")
    metadata, database = tmp_path / "metadata.tsv", tmp_path / "taxa.sqlite"
    write_tsv(metadata, ["species", "taxid"], [{"species": n, "taxid": i} for i, n in enumerate("ABCD", 1)])
    with sqlite3.connect(database) as db:
        db.executescript("CREATE TABLE species (taxid INTEGER, rank TEXT, track TEXT, spname TEXT);"
                        "CREATE TABLE merged (taxid_old INTEGER, taxid_new INTEGER);")
        db.executemany("INSERT INTO species VALUES (?, 'species', ?, ?)", [(i, str(i), n) for i, n in enumerate("ABCD", 1)])
    return tree, metadata, database


def test_taxid_resolution_rejects_aliases_mapping_to_same_species(tmp_path):
    _, metadata, database = inputs(tmp_path)
    write_tsv(metadata, ["species", "taxid"], [{"species": "A", "taxid": 11}, {"species": "Alias", "taxid": 1}])
    with sqlite3.connect(database) as db:
        db.execute("INSERT INTO merged VALUES (11, 1)")
    mapping, report = species_taxids(metadata, database, {"A", "Alias"})
    assert not mapping
    assert {r["status"] for r in report} == {"duplicate_species_taxid"}


def test_all_clades_are_reported_with_citations_and_reused_responses(tmp_path):
    tree, metadata, database = inputs(tmp_path)
    cache, out = tmp_path / "cache", tmp_path / "calibrations"
    # All nodes have fewer than 20 species. Distinct cached MRCA IDs and
    # consistent ages let us test selection through final calibration output.
    for ids, mrca, low, high in [(range(1, 5), "root", 90, 110),
                                (range(2, 5), "child", 60, 80), ([3, 4], "pair", 30, 50)]:
        cached_response(ids, cache, delay=0,
                        backend=(fake_fetch(payload(ids, mrca, low, high)), {"synthetic": True}))
    prepare(tree, metadata, database, out, cache)
    candidates = read_tsv(out / "candidates.tsv")
    assert [int(row["clade_taxa"]) for row in candidates] == [4, 3, 2]
    assert all(row["status"] == "retained" for row in candidates)
    assert len(read_tsv(out / "calibrations.tsv")) == 3
    provenance = json.loads((out / "provenance.json").read_text())
    assert provenance["status"] == "ready"
    assert provenance["eligible_nodes"] == 3
    assert provenance["internal_nodes"] == provenance["queries"] == 3
    assert provenance["policy"]["min_studies"] == 5
    assert len(read_tsv(out / "studies.tsv")) == 15
    for row in candidates:
        assert row["retrieved_at"] and row["query_taxids"] and row["mrca_id"]
        studies = [s for s in read_tsv(out / "studies.tsv") if s["node"] == row["node"]]
        assert len(studies) == 5
        assert all(json.loads(s["record_json"])["title"] == s["title"] for s in studies)
    for name, record in provenance["outputs"].items():
        assert record == file_record(out / name)
    nodes = read_tree(out / "nodes.nwk", "ABCD")
    assert len(calibration_rows(nodes, out / "calibrations.tsv")) == 3
    assert {n.name for n in nodes.traverse() if not n.is_leaf} == {r["node"] for r in candidates}
    original = read_tree(tree)
    assert all(nodes.get_distance(a, b) == original.get_distance(a, b) for a in "ABCD" for b in "ABCD")
    # A rerun preserves the original cache timestamps and all result data.
    results = {p: p.read_bytes() for p in out.iterdir() if p.name != "provenance.json"}
    cached = {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in cache.iterdir()}
    for name in ["representatives.txt", "representatives.nwk"]:
        (out / name).write_text("obsolete selection")
    prepare(tree, metadata, database, out, cache)
    assert all(p.read_bytes() == content for p, content in results.items())
    assert all((p.read_bytes(), p.stat().st_mtime_ns) == old for p, old in cached.items())
    assert not list(out.glob("representatives.*"))


def test_more_than_64_species_and_32_queries_include_small_clades(tmp_path, monkeypatch):
    tree, metadata, database = inputs(tmp_path)
    mapping = {f"Species_{i:03d}": 100 + i for i in range(70)}
    def balanced(names):
        if len(names) == 1:
            return names[0] + ":0.1"
        mid = len(names) // 2
        return f"({balanced(names[:mid])},{balanced(names[mid:])}):0.1"
    tree.write_text(balanced(sorted(mapping)) + ";\n")
    write_tsv(metadata, ["species", "taxid"], [{"species": s, "taxid": t} for s, t in mapping.items()])
    with sqlite3.connect(database) as db:
        db.executemany("INSERT INTO species VALUES (?, 'species', ?, ?)",
                       [(t, str(t), s) for s, t in mapping.items()])
    calls = []
    def fetch(url):
        ids = list(map(int, url.rsplit("/", 2)[1].split("+")))
        calls.append(tuple(ids))
        return fake_fetch(payload(ids, mrca="+".join(map(str, ids)), low=1, high=200))(url)
    monkeypatch.setattr("timetree_calibrations.nwkit_backend", lambda: (fetch, {"synthetic": True}))
    out = tmp_path / "out"
    prepare(tree, metadata, database, out, tmp_path / "cache")
    expected = {tuple(sorted(mapping[s] for s in n.leaf_names()))
                for n in read_tree(tree).traverse() if not n.is_leaf}
    assert set(calls) == expected and len(calls) == 69
    assert len(calls[0]) == 70 and min(map(len, calls)) == 2
    assert len(read_tsv(out / "calibrations.tsv")) == 69


@pytest.mark.parametrize("resolved", ["AC", ""])
def test_unresolved_child_lineages_are_reported_without_requests(tmp_path, resolved):
    tree, metadata, database = inputs(tmp_path)
    write_tsv(metadata, ["species", "taxid"],
              [{"species": s, "taxid": i} for i, s in enumerate("ABCD", 1) if s in resolved])
    cache, out = tmp_path / "cache", tmp_path / "out"
    if resolved:
        cached_response([1, 3], cache, delay=0, backend=(fake_fetch(payload([1, 3])), {}))
    prepare(tree, metadata, database, out, cache)
    report = json.loads((out / "provenance.json").read_text())
    assert report["internal_nodes"] == 3
    assert report["queries"] == report["retained"] == (1 if resolved else 0)
    rows = read_tsv(out / "candidates.tsv")
    assert len(rows) == 3
    assert all(r["reason"] == "unresolved_child_lineage" and not r["url"]
               for r in rows if r["status"] == "excluded")
    assert set(read_tree(out / "nodes.nwk").leaf_names()) == set("ABCD")


def test_conflict_preserves_bounds_and_sources_and_prevents_dating(tmp_path):
    tree, metadata, database = inputs(tmp_path)
    cache, out = tmp_path / "cache", tmp_path / "out"
    for ids, low, high in [(range(1, 5), 90, 110), (range(2, 5), 120, 150), ([3, 4], 30, 50)]:
        cached_response(ids, cache, delay=0,
                        backend=(fake_fetch(payload(ids, mrca=str(list(ids)), low=low, high=high)), {}))
    prepare(tree, metadata, database, out, cache)
    report = json.loads((out / "provenance.json").read_text())
    assert report["status"] == "conflicting_calibrations"
    assert "inconsistent ancestor/descendant" in report["error"]
    assert report["retained"] == 0 and not read_tsv(out / "calibrations.tsv")
    rows = read_tsv(out / "candidates.tsv")
    assert all(r["reason"] == "conflicting_calibration_set" for r in rows)
    assert [(float(r["min_age_ma"]), float(r["max_age_ma"])) for r in rows] == [(90, 110), (120, 150), (30, 50)]
    assert len(read_tsv(out / "studies.tsv")) == 15
    provenance = tmp_path / "tree.json"
    provenance.write_text(json.dumps({"branch_length_unit": "substitutions_per_site", "outgroup": "A",
                                     "total_gene_sites": 750}))
    # Dating must stop before resolving or launching an executable.
    with pytest.raises(ValueError, match="at least one explicit calibration"):
        date(tree, provenance, out / "calibrations.tsv", tmp_path / "dated", "must-not-run-lsd2")
    assert not (tmp_path / "dated").exists()


def test_interrupted_retrieval_resumes_with_the_completed_response(tmp_path, monkeypatch):
    tree, metadata, database = inputs(tmp_path)
    cache, out = tmp_path / "cache", tmp_path / "out"
    calls = []
    def fetch(url):
        ids = list(map(int, url.rsplit("/", 2)[1].split("+")))
        calls.append(ids)
        if len(calls) == 2:
            return {"url": url, "status_code": 503, "text": "unavailable", "error": None, "elapsed": 0}
        return fake_fetch(payload(ids, mrca=str(ids), low=1, high=200))(url)
    monkeypatch.setattr("timetree_calibrations.nwkit_backend", lambda: (fetch, {"synthetic": True}))
    with pytest.raises(ValueError, match="transport failed"):
        prepare(tree, metadata, database, out, cache)
    assert not (out / "calibrations.tsv").exists()
    completed = next(p for p in cache.glob("*.json") if not p.name.endswith(".failure.json"))
    record, timestamp = completed.read_bytes(), completed.stat().st_mtime_ns
    prepare(tree, metadata, database, out, cache)
    assert calls == [[1, 2, 3, 4], [2, 3, 4], [2, 3, 4], [3, 4]]
    assert completed.read_bytes() == record and completed.stat().st_mtime_ns == timestamp
    assert len(read_tsv(out / "calibrations.tsv")) == 3


def test_preparation_and_real_lsd2_keep_topology_and_substitution_tree(tmp_path):
    tree, metadata, database = inputs(tmp_path)
    original = tree.read_bytes()
    cache = tmp_path / "cache"
    cached_response(range(1, 5), cache, delay=0, backend=(fake_fetch(payload(range(1, 5))), {"synthetic": True}))
    cache_missing_estimates(tree, dict(zip("ABCD", range(1, 5))), cache)
    out = tmp_path / "calibrations"
    prepare(tree, metadata, database, out, cache)
    assert json.loads((out / "provenance.json").read_text())["status"] == "ready"
    assert len(read_tsv(out / "calibrations.tsv")) == 1
    lsd2 = os.environ.get("LSD2_BIN")
    if lsd2:
        provenance = tmp_path / "tree.json"
        provenance.write_text(json.dumps({"branch_length_unit": "substitutions_per_site", "outgroup": "A", "mean_gene_length": 250, "total_gene_sites": 750}))
        date(tree, provenance, out / "calibrations.tsv", tmp_path / "dated", lsd2)
        dated = read_tree(tmp_path / "dated/species_tree.dated.nwk", "ABCD")
        assert 90 - 1e-3 <= dated.get_distance("A", "B") / 2 <= 110 + 1e-3
        result = json.loads((tmp_path / "dated/provenance.json").read_text())
        assert result["native_report"]["unique_scale"] is False
        assert result["native_report"]["rate_substitutions_per_site_per_ma"] is None
        assert result["native_report"]["root_date_interval"] == [-110, -90]
        assert result["confidence_intervals"] is False
        assert result["review_status"] == "needs_review"
        assert "midpoint" in result["scale_selection"]
    assert tree.read_bytes() == original
