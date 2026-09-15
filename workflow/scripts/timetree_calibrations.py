#!/usr/bin/env python3
"""Fetch auditable secondary calibrations with nwkit's TimeTree HTTP client.

The pinned nwkit mcmctree client is used directly: the JSON endpoint preserves
MRCA identifiers, missing taxa and studies that its calibration Newick omits.
No nwkit internals are modified and no alignment or branch length is optimized.
"""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import sqlite3
import time

from common import atomic_writer, file_record, now, write_json, write_tsv
from date_phylogeny import calibration_rows, newick_text
from infer_phylogeny import read_tree


NWKIT_COMMIT = "db5b8a32c7608248db9f2b7b8aed16376779c5fb"
NWKIT_MODULE_SHA256 = "56dd8efd5367ca50143ef9df783441d9693da4fa87a6639b5704eef575d05b79"
NWKIT_CONDA_MODULE_SHA256 = "cbcd120846bdbac1c596056444c7c8622d5386ede7f63277f581bd634c5c9cf2"
ENDPOINT = "https://timetree.org/api/mrca/id/"
FIELDS = ["taxa", "min_age_ma", "max_age_ma", "source"]
MIN_STUDIES = 5
REQUEST_DELAY_SECONDS = 1.0


def nwkit_backend():
    import nwkit
    from nwkit import mcmctree
    # This is deliberately a narrow, pinned private-API boundary. A changed
    # implementation must be reviewed/tested before it is used for retrieval.
    record = file_record(mcmctree.__file__)
    if record["sha256"] not in {NWKIT_MODULE_SHA256, NWKIT_CONDA_MODULE_SHA256}:
        raise ValueError("unsupported nwkit TimeTree client; use workflow/envs/timetree.yaml")
    return mcmctree._fetch_timetree_url, {"version": nwkit.__version__,
                                        "source": "bioconda::nwkit=0.27.0" if record["sha256"] == NWKIT_CONDA_MODULE_SHA256 else NWKIT_COMMIT,
                                        "module": record}


def cached_response(taxids, cache_dir, delay=REQUEST_DELAY_SECONDS, backend=None):
    url = ENDPOINT + "+".join(str(t) for t in sorted(set(taxids))) + "/json"
    key = hashlib.sha256(url.encode()).hexdigest()
    path = Path(cache_dir) / (key + ".json")
    if path.exists():
        record = json.loads(path.read_text())
        if (record.get("url") != url or record.get("schema") != 1
                or record.get("response_sha256") != hashlib.sha256(record["text"].encode()).hexdigest()
                or record.get("status_code") != 200):
            raise ValueError(f"invalid TimeTree cache record: {path}")
    else:
        fetch, tool = backend if backend is not None else nwkit_backend()
        time.sleep(delay)  # Serial retrieval; avoid bursts against the public API.
        response = fetch(url)
        record = {**response, "error": str(response["error"]) if response["error"] else None,
                  "schema": 1, "retrieved_at": now(), "nwkit": tool,
                  "response_sha256": hashlib.sha256(response["text"].encode()).hexdigest()}
        if response["status_code"] != 200 or response["error"]:
            write_json(path.with_suffix(".failure.json"), record)
            raise ValueError(f"TimeTree transport failed; retry later: {url}")
        # Do not freeze a proxy/HTML error as a successful empty biological result.
        payload = json.loads(response["text"])
        if not isinstance(payload, dict) or not {"found_ids", "missing_ids"} <= payload.keys():
            raise ValueError(f"unexpected TimeTree JSON schema: {url}")
        write_json(path, record)
    return json.loads(record["text"]), record, path


def species_taxids(metadata, database, names):
    """Resolve input/merged/subspecies taxids to species using the local snapshot."""
    identifiers = defaultdict(set)
    with open(metadata, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("species") in names and row.get("taxid", "").isdigit():
                identifiers[row["species"]].add(int(row["taxid"]))
    mapping, report = {}, []
    with sqlite3.connect(Path(database).resolve().as_uri() + "?mode=ro", uri=True) as db:
        for name in sorted(names):
            resolved = set()
            for taxid in identifiers[name]:
                row = db.execute("SELECT taxid, rank, track, spname FROM species WHERE taxid=?", (taxid,)).fetchone()
                if row is None:
                    merged = db.execute("SELECT taxid_new FROM merged WHERE taxid_old=?", (taxid,)).fetchone()
                    if merged:
                        row = db.execute("SELECT taxid, rank, track, spname FROM species WHERE taxid=?", merged).fetchone()
                if row is None:
                    continue
                lineage = [row[0]] + [int(t) for t in row[2].split(",") if t.isdigit()]
                for ancestor in lineage:
                    sp = db.execute("SELECT taxid, rank, spname FROM species WHERE taxid=?", (ancestor,)).fetchone()
                    if sp and sp[1] == "species":
                        resolved.add((sp[0], sp[2]))
                        break
            if len(resolved) == 1:
                mapping[name] = next(iter(resolved))[0]
            report.append({"species": name, "input_taxids": ",".join(map(str, sorted(identifiers[name]))),
                           "species_taxid": mapping.get(name, ""),
                           "ncbi_name": next(iter(resolved))[1] if len(resolved) == 1 else "",
                           "status": "resolved" if len(resolved) == 1 else "unresolved_or_conflicting"})
    duplicates = {taxid for taxid, count in Counter(mapping.values()).items() if count > 1}
    for row in report:
        if row["species_taxid"] in duplicates:
            mapping.pop(row["species"], None)
            row["status"] = "duplicate_species_taxid"
    return mapping, report


def node_queries(tree, mapping):
    """Visit every internal node, using all resolved descendant species.

    Canonical child order makes node labels and requests independent of Newick
    child order. Only the in-memory tree is relabelled; lengths/root are retained.
    """
    counts, first_tip = {}, {}
    for node in tree.traverse("postorder"):
        if node.is_leaf:
            counts[node] = 1
            first_tip[node] = node.name
        else:
            counts[node] = sum(counts[c] for c in node.children)
            first_tip[node] = min(first_tip[c] for c in node.children)
            node.children.sort(key=lambda c: first_tip[c])
    names = set(tree.leaf_names())
    internal = [n for n in tree.traverse("postorder") if not n.is_leaf]
    for i, node in enumerate(internal, 1):
        node.name = f"T{i:06d}"
        while node.name in names:
            node.name = "_" + node.name
    for node in tree.traverse("preorder"):
        if node.is_leaf:
            continue
        children = [sorted(s for s in c.leaf_names() if s in mapping) for c in node.children]
        yield node, {"node": node.name, "clade_taxa": counts[node],
                     "query_taxa": sorted(s for child in children for s in child), "children": children}


def study_records(payload):
    """Deduplicate named bibliographic records, not trees within a paper."""
    studies = payload.get("study_data", [])
    if not isinstance(studies, list) or any(not isinstance(s, dict) for s in studies):
        return None
    unique = {}
    for study in studies:
        if not isinstance(study.get("title"), str) or not study["title"].strip():
            continue
        key = tuple(" ".join(str(study.get(k) or "").split()).casefold()
                    for k in ["title", "pub_year", "first_author"])
        unique.setdefault(key, study)
    return [unique[key] for key in sorted(unique)]


def interpret_response(payload, query, mapping):
    """Require observed taxa on every child lineage, finite bounds and studies."""
    row = {**query, "status": "excluded", "reason": "", "mrca_id": "", "studies": 0,
           "query_taxids": sorted(mapping[s] for s in query["query_taxa"]),
           "used_taxa": [], "used_taxids": [], "missing_taxids": [], "study_records": [],
           "min_age_ma": None, "age_ma": None, "max_age_ma": None,
           "url": "", "retrieved_at": "", "cache": None}
    if payload is None:
        row["reason"] = "unresolved_child_lineage"
        return row
    # Preserve available ages and citations even when a candidate is excluded.
    for field, key in [("min_age_ma", "precomputed_ci_low"), ("age_ma", "precomputed_age"),
                       ("max_age_ma", "precomputed_ci_high")]:
        try:
            value = float(payload[key])
            if math.isfinite(value):
                row[field] = value
        except (KeyError, TypeError, ValueError):
            pass
    records = study_records(payload)
    row["study_records"] = records or []
    row["studies"] = len(row["study_records"])
    mrca = payload.get("mrca_id", payload.get("mrca_ttid"))
    if mrca is not None and str(mrca) not in {"", "-1", "0"}:
        row["mrca_id"] = str(mrca)
    requested = set(map(str, row["query_taxids"]))
    found_raw, missing_raw = payload.get("found_ids"), payload.get("missing_ids")
    if not isinstance(found_raw, list) or not isinstance(missing_raw, list):
        row["reason"] = "unexpected_taxon_response"
        return row
    found = {str(t) for t in found_raw if t is not None}
    missing = {str(t) for t in missing_raw if t is not None}
    if found | missing != requested or found & missing:
        row["reason"] = "unexpected_taxon_response"
        return row
    row["used_taxids"] = sorted(map(int, found))
    row["missing_taxids"] = sorted(map(int, missing))
    row["used_taxa"] = sorted(s for s in query["query_taxa"] if str(mapping[s]) in found)
    if any(not any(str(mapping[s]) in found for s in child) for child in query["children"]):
        row["reason"] = "missing_child_lineage"
        return row
    if payload.get("error"):
        row["reason"] = "no_timetree_estimate"
        return row
    if not row["mrca_id"]:
        row["reason"] = "missing_mrca_identifier"
        return row
    if records is None:
        row["reason"] = "unexpected_study_response"
        return row
    if row["studies"] < MIN_STUDIES:
        row["reason"] = "insufficient_studies"
        return row
    low, age, high = (row[k] for k in ["min_age_ma", "age_ma", "max_age_ma"])
    if None in (low, age, high) or not 0 < low <= age <= high or low == high:
        row["reason"] = "invalid_or_point_only_bounds"
        return row
    row.update(status="candidate", min_age_ma=low, max_age_ma=high)
    return row


def select_calibrations(tree, candidates, output):
    # Equal ages alone are not evidence of a duplicate. Use the returned MRCA ID.
    assignments = defaultdict(list)
    for row in candidates:
        if row["mrca_id"] and row["reason"] not in {
                "unexpected_taxon_response", "missing_child_lineage", "no_timetree_estimate",
                "unresolved_child_lineage", "input_mrca_mismatch"}:
            assignments[row["mrca_id"]].append(row)
    for rows in assignments.values():
        if len(rows) > 1:
            for row in rows:
                row.update(status="excluded", reason="ambiguous_repeated_timetree_mrca")
    chosen = []
    for row in candidates:
        if row["status"] == "candidate":
            row["status"] = "retained"
            chosen.append({"taxa": ",".join(row["used_taxa"]), "min_age_ma": row["min_age_ma"],
                           "max_age_ma": row["max_age_ma"],
                           "source": f'TimeTree secondary calibration; MRCA={row["mrca_id"]}; '
                                     f'{row["url"]}; retrieved={row["retrieved_at"]}'})
    write_tsv(output, FIELDS, chosen)
    try:
        calibration_rows(tree, output)
    except ValueError as exc:
        # Preserve the candidate report, but prevent an invalid automatic dating.
        write_tsv(output, FIELDS, [])
        for row in candidates:
            if row["status"] == "retained":
                row.update(status="excluded", reason="conflicting_calibration_set")
        return {"status": "no_valid_calibrations" if not chosen else "conflicting_calibrations",
                "error": str(exc), "retained": 0}
    return {"status": "ready", "retained": len(chosen)}


def prepare(tree, metadata, taxonomy_db, outdir, cache_dir):
    source = tree
    tree = read_tree(tree)
    if len(tree.children) != 2:
        raise ValueError("TimeTree calibration requires the rooted species tree")
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    mapping, taxa_report = species_taxids(metadata, taxonomy_db, set(tree.leaf_names()))
    candidates, queries = [], 0
    for node, query in node_queries(tree, mapping):
        if not all(query["children"]):
            candidates.append(interpret_response(None, query, mapping))
            continue
        payload, raw, path = cached_response([mapping[s] for s in query["query_taxa"]], cache_dir)
        queries += 1
        row = interpret_response(payload, query, mapping)
        # The actually observed species must still identify this exact input node.
        if row["status"] == "candidate" and tree.common_ancestor(row["used_taxa"]) is not node:
            row.update(status="excluded", reason="input_mrca_mismatch")
        row.update(url=raw["url"], retrieved_at=raw["retrieved_at"], cache=file_record(path))
        candidates.append(row)
    status = select_calibrations(tree, candidates, out / "calibrations.tsv")
    write_tsv(out / "taxa.tsv", list(taxa_report[0]), taxa_report)
    with atomic_writer(out / "nodes.nwk") as handle:
        handle.write(newick_text(tree))
    write_json(out / "candidates.json", candidates)
    columns = ["node", "clade_taxa", "query_taxa", "query_taxids", "used_taxa", "used_taxids", "missing_taxids",
               "mrca_id", "studies", "min_age_ma", "age_ma", "max_age_ma", "status", "reason", "url", "retrieved_at"]
    write_tsv(out / "candidates.tsv", columns,
              [{k: ",".join(map(str, row[k])) if isinstance(row[k], list) else row[k] for k in columns}
               for row in candidates])
    fields = ["node", "mrca_id", "title", "pub_year", "first_author", "record_json"]
    write_tsv(out / "studies.tsv", fields,
              ({"node": row["node"], "mrca_id": row["mrca_id"],
                **{k: study.get(k, "") for k in ["title", "pub_year", "first_author"]},
                "record_json": json.dumps(study, ensure_ascii=False, sort_keys=True)}
               for row in candidates for study in row["study_records"]))
    write_json(out / "provenance.json", {**status, "created_at": now(),
        "policy": {"candidate_nodes": "all internal nodes", "min_studies": MIN_STUDIES,
                   "study_count": "distinct named bibliographic records (title, pub_year, first_author)",
                   "duplicate_mrca": "exclude all ambiguous assignments", "conflicting_bounds": "stop dating"},
        "request_delay_seconds": REQUEST_DELAY_SECONDS, "code": file_record(__file__),
        "tree": file_record(source), "metadata": file_record(metadata), "taxonomy": file_record(taxonomy_db),
        "nwkit_environment": "workflow/envs/timetree.yaml", "internal_nodes": len(candidates),
        "eligible_nodes": queries, "queries": queries,
        "node_eligibility": "internal node with resolved species on every child lineage",
        "query_taxa": "all resolved descendant species", "query_order": "canonical preorder by first tip label",
        "calibration_type": "TimeTree secondary bounds interpreted as hard LSD2 bounds",
        "topology_check": "sampled child-lineage coverage and distinct MRCA IDs; not proof of full topological concordance",
        "outputs": {name: file_record(out / name) for name in
                    ["calibrations.tsv", "candidates.tsv", "candidates.json", "taxa.tsv", "nodes.nwk", "studies.tsv"]},
        "confidence_intervals_propagated": False, "cache_dir": str(Path(cache_dir).resolve())})
    for name in ["representatives.txt", "representatives.nwk"]:
        (out / name).unlink(missing_ok=True)
    print(f"TimeTree: {queries} queries, {status['retained']} calibrations, {status['status']}", flush=True)
    if status.get("error"):
        print(status["error"], flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["tree", "metadata", "taxonomy-db", "outdir", "cache-dir"]:
        parser.add_argument("--" + name, required=True)
    args = vars(parser.parse_args())
    prepare(**args)
