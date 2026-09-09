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

from common import atomic_writer, file_record, now, read_tsv, sha256, write_json, write_tsv
from date_phylogeny import calibration_rows
from infer_phylogeny import read_tree


NWKIT_COMMIT = "db5b8a32c7608248db9f2b7b8aed16376779c5fb"
NWKIT_MODULE_SHA256 = "56dd8efd5367ca50143ef9df783441d9693da4fa87a6639b5704eef575d05b79"
ENDPOINT = "https://timetree.org/api/mrca/id/"
FIELDS = ["taxa", "min_age_ma", "max_age_ma", "source"]


def nwkit_backend():
    import nwkit
    from nwkit import mcmctree
    # This is deliberately a narrow, pinned private-API boundary. A changed
    # implementation must be reviewed/tested before it is used for retrieval.
    record = file_record(mcmctree.__file__)
    if record["sha256"] != NWKIT_MODULE_SHA256:
        raise ValueError("unsupported nwkit TimeTree client; use workflow/envs/timetree.yaml")
    return mcmctree._fetch_timetree_url, {"version": nwkit.__version__, "commit": NWKIT_COMMIT,
                                        "module": record}


def cached_response(taxids, cache_dir, offline=False, delay=1.0, backend=None):
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
        if offline:
            raise ValueError(f"offline TimeTree cache miss: {url}")
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


def representative_plan(tree, mapping, coverage, maximum, representatives=None):
    """Partition by the largest clade, then take its best-covered resolved tip.

    Uses O(n * maximum) stored representative membership, not all descendant
    sets for every node. The original tree, rooting and lengths are preserved.
    """
    counts, best = {}, {}
    for i, node in enumerate(tree.traverse("postorder"), 1):
        if node.is_leaf:
            counts[node] = 1
            best[node] = node.name if node.name in mapping else None
        else:
            node.name = f"T{i:06d}"
            counts[node] = sum(counts[c] for c in node.children)
            options = [best[c] for c in node.children if best[c] is not None]
            best[node] = min(options, key=lambda n: (-coverage.get(n, 0), n)) if options else None
    if representatives:
        selected = [s.strip() for s in Path(representatives).read_text().splitlines() if s.strip()]
        if len(selected) != len(set(selected)) or not 2 <= len(selected) <= maximum:
            raise ValueError("representatives must contain 2..max_representatives unique species")
        if set(selected) - mapping.keys():
            raise ValueError("representatives include absent or unresolved species taxids")
    else:
        frontier = [tree]
        while True:
            options = [n for n in frontier if not n.is_leaf and
                       len(frontier) + len(n.children) - 1 <= maximum]
            if not options:
                break
            node = min(options, key=lambda n: (-counts[n], best[n] or "", n.name))
            frontier.remove(node)
            frontier.extend(node.children)
        selected = sorted(best[n] for n in frontier if best[n] is not None)
    membership = {}
    for node in tree.traverse("postorder"):
        membership[node] = ([node.name] if node.name in selected else []) if node.is_leaf else sorted(
            s for c in node.children for s in membership[c])
    if any(not membership[c] for c in tree.children):
        raise ValueError("representatives must include resolved species on both sides of the root")
    return selected, counts, membership


def interpret_response(payload, query, mapping, min_studies):
    """Require observed taxa on every child lineage, finite bounds and studies."""
    row = {**query, "status": "excluded", "reason": "", "mrca_id": "", "studies": 0,
           "used_taxa": [], "min_age_ma": None, "max_age_ma": None}
    requested = {str(mapping[s]) for s in query["query_taxa"]}
    found_raw, missing_raw = payload.get("found_ids"), payload.get("missing_ids")
    if not isinstance(found_raw, list) or not isinstance(missing_raw, list):
        row["reason"] = "unexpected_taxon_response"
        return row
    found = {str(t) for t in found_raw if t is not None}
    missing = {str(t) for t in missing_raw if t is not None}
    if found - requested or found & missing:
        row["reason"] = "unexpected_taxon_response"
        return row
    row["used_taxa"] = sorted(s for s in query["query_taxa"] if str(mapping[s]) in found)
    if any(not any(str(mapping[s]) in found for s in child) for child in query["children"]):
        row["reason"] = "missing_child_lineage"
        return row
    if payload.get("error"):
        row["reason"] = "no_timetree_estimate"
        return row
    mrca = payload.get("mrca_id", payload.get("mrca_ttid"))
    if mrca is None or str(mrca) in {"", "-1", "0"}:
        row["reason"] = "missing_mrca_identifier"
        return row
    row["mrca_id"] = str(mrca)
    studies = payload.get("study_data", [])
    if not isinstance(studies, list) or any(not isinstance(s, dict) for s in studies):
        row["reason"] = "unexpected_study_response"
        return row
    # Count distinct bibliographic records, not multiple trees from one paper.
    row["studies"] = len({(s.get("title"), s.get("pub_year"), s.get("first_author")) for s in studies if s.get("title")})
    if row["studies"] < min_studies:
        row["reason"] = "insufficient_studies"
        return row
    try:
        low, age, high = (float(payload[k]) for k in ["precomputed_ci_low", "precomputed_age", "precomputed_ci_high"])
        if not all(math.isfinite(v) for v in (low, age, high)) or not 0 < low <= age <= high or low == high:
            raise ValueError("invalid bounds")
    except (KeyError, TypeError, ValueError):
        row["reason"] = "invalid_or_point_only_bounds"
        return row
    row.update(status="candidate", min_age_ma=low, max_age_ma=high)
    return row


def select_calibrations(tree, candidates, output):
    # Equal ages alone are not evidence of a duplicate. Use the returned MRCA ID.
    assignments = defaultdict(list)
    for row in candidates:
        if row["mrca_id"]:
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


def prepare(tree, metadata, taxonomy_db, coverage, outdir, cache_dir, settings, representatives=None):
    settings = dict(settings)
    if "min_clade_taxa" in settings:
        raise ValueError("TimeTree min_clade_taxa was removed; remove it from settings. Query limits control workload")
    for key in ["max_representatives", "max_queries", "min_studies"]:
        if type(settings[key]) is not int or settings[key] < 1:
            raise ValueError(f"TimeTree {key} must be a positive integer")
    if settings["max_representatives"] < 2 or type(settings["offline"]) is not bool:
        raise ValueError("TimeTree requires max_representatives >= 2 and a boolean offline setting")
    if not math.isfinite(settings["request_delay_seconds"]) or settings["request_delay_seconds"] < 0:
        raise ValueError("TimeTree request delay must be finite and nonnegative")
    source = tree
    tree = read_tree(tree)
    if len(tree.children) != 2:
        raise ValueError("TimeTree calibration requires the rooted species tree")
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    mapping, taxa_report = species_taxids(metadata, taxonomy_db, set(tree.leaf_names()))
    counts = {r["species"]: int(r["gene_trees"]) for r in read_tsv(coverage)}
    selected, sizes, members = representative_plan(tree, mapping, counts, settings["max_representatives"], representatives)
    write_tsv(out / "taxa.tsv", list(taxa_report[0]), taxa_report)
    with atomic_writer(out / "representatives.txt") as handle:
        handle.write("\n".join(selected) + "\n")
    skeleton = tree.copy()
    skeleton.prune(selected, preserve_branch_length=True)
    with atomic_writer(out / "representatives.nwk") as handle:
        handle.write(skeleton.write(parser=1, format_root_node=True) + "\n")
    eligible = [n for n in tree.traverse() if not n.is_leaf and all(members[c] for c in n.children)]
    eligible.sort(key=lambda n: (-sizes[n], tuple(members[n])))
    candidates = []
    for node in eligible[:settings["max_queries"]]:
        query = {"node": node.name, "clade_taxa": sizes[node], "query_taxa": members[node],
                 "children": [members[c] for c in node.children]}
        payload, raw, path = cached_response([mapping[s] for s in members[node]], cache_dir,
            offline=settings["offline"], delay=settings["request_delay_seconds"])
        row = interpret_response(payload, query, mapping, settings["min_studies"])
        # The actually observed species must still identify this exact input node.
        if row["status"] == "candidate" and tree.common_ancestor(row["used_taxa"]) is not node:
            row.update(status="excluded", reason="input_mrca_mismatch")
        row.update(url=raw["url"], retrieved_at=raw["retrieved_at"], cache=file_record(path))
        candidates.append(row)
    status = select_calibrations(tree, candidates, out / "calibrations.tsv")
    write_json(out / "candidates.json", candidates)
    write_tsv(out / "candidates.tsv", ["node", "clade_taxa", "query_taxa", "used_taxa", "mrca_id", "studies",
              "min_age_ma", "max_age_ma", "status", "reason", "url"],
              [{k: ",".join(row[k]) if k in {"query_taxa", "used_taxa"} else row[k] for k in
                ["node", "clade_taxa", "query_taxa", "used_taxa", "mrca_id", "studies", "min_age_ma", "max_age_ma",
                 "status", "reason", "url"]} for row in candidates])
    write_json(out / "provenance.json", {**status, "created_at": now(), "settings": settings,
        "tree": file_record(source), "metadata": file_record(metadata), "taxonomy": file_record(taxonomy_db),
        "coverage": file_record(coverage), "representatives_file": file_record(representatives) if representatives else None,
        "nwkit_commit": NWKIT_COMMIT, "eligible_nodes": len(eligible), "queries": len(candidates),
        "node_eligibility": "internal node with representatives on every child lineage",
        "query_ranking": ["clade_taxa descending", "representative species labels ascending (lexicographic)"],
        "calibration_type": "TimeTree secondary bounds interpreted as hard treePL bounds",
        "topology_check": "sampled child-lineage coverage and distinct MRCA IDs; not proof of full topological concordance",
        "confidence_intervals_propagated": False, "cache_dir": str(Path(cache_dir).resolve())})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["tree", "metadata", "taxonomy-db", "coverage", "outdir", "cache-dir", "settings"]:
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--representatives")
    args = vars(parser.parse_args())
    args["settings"] = json.loads(args["settings"])
    prepare(**args)
