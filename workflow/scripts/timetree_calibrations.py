#!/usr/bin/env python3
"""Collect secondary TimeTree bounds on eligible nodes of the rooted species tree."""
import argparse
from contextlib import closing
from csv import DictReader
import math
from pathlib import Path

from common import file_record, now, read_tsv, write_json, write_tsv
from infer_phylogeny import read_tree
from taxonomy_check import Taxonomy

ENDPOINT = "https://timetree.org/api"
MIN_CLADE_PROP = 0.2
NWKIT_VERSION = "0.46.11"
FIELDS = ["left_species", "right_species", "min_age_ma", "max_age_ma", "clade_tips", "query_url"]


def query_mrca(taxids, cache):
    from nwkit.mcmctree import _fetch_timetree_url_cached

    query_ids = "+".join(map(str, sorted(set(taxids))))
    url = f"{ENDPOINT}/mrca/id/{query_ids}"
    response = _fetch_timetree_url_cached(url, response_cache=cache)
    if response["error"] is not None or response["status_code"] in {429, 500, 502, 503, 504}:
        raise ValueError(f"TimeTree request failed after retries: {url}")
    return url, response


def mrca_key(response):
    if response["status_code"] != 200:
        return None
    lines = response["text"].replace("<br>", "\n").splitlines()
    row = next(DictReader(lines), {})
    taxid = (row.get("taxon_id") or "").strip()
    if taxid.isdecimal() and int(taxid) > 0:
        return int(taxid)
    if taxid == "0":  # Unnamed TimeTree nodes share this ID.
        try:
            ages = tuple(float(row[key]) for key in ("precomputed_age", "precomputed_ci_low", "precomputed_ci_high"))
        except (KeyError, TypeError, ValueError):
            return None
        if all(math.isfinite(age) and age >= 0 for age in ages) and ages[0] > 0:
            return ages
    return None


def collect_bounds(tree, tip_ranks):
    from nwkit.mcmctree import (
        SEARCH_RANKS, _build_timetree_rank_attempt, _constraint_from_timetree_response,
        remove_constraint_equal_upper,
    )
    from nwkit.time_tree import parse_mcmctree_calibration
    from nwkit.util import get_subtree_leaf_name_sets

    leaf_sets = get_subtree_leaf_name_sets(tree)
    threshold = MIN_CLADE_PROP * len(leaf_sets[tree])
    candidates = []
    for node in tree.traverse():
        if node.is_leaf:
            continue
        node.name = ""
        if node.is_root or len(leaf_sets[node]) >= threshold:
            candidates.append(node)

    # Use the same clade size condition and rank search as NWKIT, applying the
    # size condition early so that nodes below it do not need an API query.
    args = argparse.Namespace(timetree="ci", lower_tail_prob="0.025", upper_tail_prob="0.025")
    cache, urls = {}, {}
    for node in candidates:
        tips = sorted(leaf_sets[node])
        lineage_taxids = [tip_ranks[tip] for tip in tips]
        context = dict(node=node, leaf_names=tips, lineage_taxids=lineage_taxids,
                       leaf_rank_pairs=list(zip(tips, lineage_taxids)))
        for rank in SEARCH_RANKS:
            attempt = _build_timetree_rank_attempt(context, rank, ENDPOINT, leaf_sets)
            if attempt is None:
                continue
            tips_by_taxid = attempt["taxid_to_species_labels"]
            url, response = query_mrca(tips_by_taxid, cache)
            try:
                label = _constraint_from_timetree_response(attempt, response, None, leaf_sets, args)
            except ValueError as exc:
                raise ValueError(f"Could not read TimeTree response from {url}: {exc}") from exc
            if not label:
                continue
            # TimeTree can ignore absent TaxIDs. Require the combined query to
            # reach a different node from the reference found for either side.
            combined = mrca_key(response)
            sides = []
            for child in node.children:
                ids = [taxid for taxid, names in tips_by_taxid.items() if not leaf_sets[child].isdisjoint(names)]
                _, side_response = query_mrca(ids, cache)
                sides.append(mrca_key(side_response))
            if combined is None or all(side is None for side in sides) or combined in sides:
                print(f"Skipping TimeTree interval at {url}: the reference did not resolve both sides of the split", flush=True)
                continue
            node.name = label
            urls[node] = url
            break

    remove_constraint_equal_upper(tree)
    rows = []
    for node in candidates:
        interval = parse_mcmctree_calibration(node.name)
        if interval is None:
            continue
        left, right = sorted(min(leaf_sets[child]) for child in node.children)
        rows.append(dict(left_species=left, right_species=right, min_age_ma=interval["lower"],
                         max_age_ma=interval["upper"], clade_tips=len(leaf_sets[node]), query_url=urls[node]))
    return rows, len(candidates)


def prepare(tree, samples, taxonomy, outdir):
    import nwkit

    if nwkit.__version__ != NWKIT_VERSION:
        raise ValueError(f"TimeTree calibrations require pinned NWKIT {NWKIT_VERSION}")
    source, tree = tree, read_tree(tree)
    if any(len(node.children) != 2 for node in tree.traverse() if not node.is_leaf):
        raise ValueError("TimeTree calibration requires a rooted bifurcating tree")
    taxids = {}
    for row in read_tsv(samples):
        tip, taxid = row["species"], int(row["taxid"])
        if taxids.setdefault(tip, taxid) != taxid:
            raise ValueError(f"conflicting frozen NCBI TaxIDs for dating tip {tip}")
    ranks = {}
    database = Taxonomy(taxonomy)
    with closing(database.db):
        for tip in tree.leaf_names():
            lineage = database.lineage(taxids[tip]) if tip in taxids else {}
            if not lineage:
                raise ValueError(f"no frozen NCBI taxonomy for dating tip {tip}")
            ranks[tip] = {rank: taxon[0] for rank, taxon in lineage.items()}

    rows, eligible = collect_bounds(tree, ranks)
    outdir = Path(outdir)
    write_tsv(outdir / "calibrations.tsv", FIELDS, rows)
    write_json(outdir / "provenance.json", {
        "created_at": now(), "status": "ready" if rows else "no_applicable_bounds",
        "tree": file_record(source), "samples": file_record(samples), "taxonomy": file_record(taxonomy),
        "source": ENDPOINT, "nwkit": nwkit.__version__, "min_clade_prop": MIN_CLADE_PROP,
        "eligible_nodes": eligible, "calibrated_nodes": len(rows),
        "calibrated_root": any(row["clade_tips"] == len(ranks) for row in rows),
    })


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["tree", "samples", "taxonomy", "outdir"]:
        parser.add_argument("--" + name, required=True)
    prepare(**vars(parser.parse_args()))
