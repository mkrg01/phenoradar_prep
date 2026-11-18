#!/usr/bin/env python3
"""Place AngioCal v1.0 fossil minima using the frozen NCBI TaxIDs."""
import argparse
from collections import Counter, defaultdict
from dataclasses import asdict
from pathlib import Path
import sys
from types import SimpleNamespace

from common import file_record, read_tsv, write_json, write_tsv
from infer_phylogeny import read_tree
from taxonomy_check import Taxonomy
from nwkit_runtime import nwkit_backend

ANGIOSPERMS = 3398
CALIBRATION_FIELDS = ["left_species", "right_species", "min_age_ma", "fossil_ids"]


class FossilTaxonomy(Taxonomy):
    """Supply the name lookups used by NWKIT from the same read-only snapshot."""

    def get_name_translator(self, names):
        matches = {}
        for name in names:
            ids = [row[0] for row in self.db.execute(
                "SELECT taxid FROM species WHERE spname=? COLLATE NOCASE", (name,))]
            if not ids:
                ids = [row[0] for row in self.db.execute(
                    "SELECT taxid FROM synonym WHERE spname=? COLLATE NOCASE", (name,))]
            matches[name] = sorted({self.canonical(taxid) for taxid in ids})
        return matches


def tip_lineages(tree, samples, taxonomy):
    taxids = {}
    for row in read_tsv(samples):
        name, taxid = row["species"], int(row["taxid"])
        if taxids.setdefault(name, taxid) != taxid:
            raise ValueError(f"conflicting frozen NCBI TaxIDs for dating tip {name}")
    lineages = {}
    for tip in tree.leaf_names():
        row = taxonomy.node(taxids[tip]) if tip in taxids else None
        if row is None:
            raise ValueError(f"no frozen NCBI taxonomy for dating tip {tip}")
        lineage = {row[0], *map(int, row[3].split(","))}
        if ANGIOSPERMS not in lineage:
            raise ValueError(f"AngioCal dating currently supports angiosperms only: {tip}")
        lineages[tip] = lineage
    return lineages


def prepare(tree, samples, taxonomy, outdir, mapping=None, workbook=None):
    nwkit = nwkit_backend()
    from nwkit.angiocal import read_angiocal
    from nwkit.angiocal_constraints import (
        _taxonomy_stem_target, _taxonomy_tip_clades, read_calibration_map,
    )
    from nwkit.util import get_subtree_leaf_name_sets

    workbook = workbook or Path(sys.prefix) / "share/angiocal/v1.0.xls"
    dataset = read_angiocal(workbook, official=True)
    source, tree = tree, read_tree(tree)
    if any(len(node.children) != 2 for node in tree.traverse() if not node.is_leaf):
        raise ValueError("AngioCal dating requires a rooted bifurcating tree")
    explicit = read_calibration_map(
        mapping, dataset.records, tree, SimpleNamespace(species_parser="taxonomic"))
    leaf_sets = get_subtree_leaf_name_sets(tree)
    db = FossilTaxonomy(taxonomy)
    rows, placements = [], defaultdict(list)
    try:
        lineages = tip_lineages(tree, samples, db)
        taxon_clades = _taxonomy_tip_clades(lineages)
        conflicts = {}
        for fossil in dataset.records:
            row = asdict(fossil)
            row.update(status="excluded", reason="", mapping_method="", left_species="", right_species="")
            rows.append(row)
            if fossil.fossil_id in explicit:
                node, reason, method = explicit[fossil.fossil_id], "", "explicit_map"
            elif fossil.placement == "stem":
                node, reason, method = _taxonomy_stem_target(
                    tree, fossil, db, (None, lineages, ""), leaf_sets, taxon_clades, conflicts)
            else:
                node, reason, method = None, "crown_requires_reviewed_anchors", ""
            row.update(reason=reason, mapping_method=method)
            if node is not None:
                sides = sorted(min(leaf_sets[child]) for child in node.children)
                row.update(status="supporting", left_species=sides[0], right_species=sides[1])
                placements[node].append(row)
    finally:
        db.db.close()

    calibrations = []
    for fossils in placements.values():
        minimum = max(row["minimum_age_ma"] for row in fossils)
        for row in fossils:
            if row["minimum_age_ma"] == minimum:
                row["status"] = "applied"
        calibrations.append(dict(left_species=fossils[0]["left_species"],
                                 right_species=fossils[0]["right_species"], min_age_ma=minimum,
                                 fossil_ids=",".join(row["fossil_id"] for row in fossils)))
    calibrations.sort(key=lambda row: (row["left_species"], row["right_species"]))
    outdir = Path(outdir)
    write_tsv(outdir / "fossils.tsv", list(rows[0]), rows)
    write_tsv(outdir / "calibrations.tsv", CALIBRATION_FIELDS, calibrations)
    write_json(outdir / "provenance.json", {
        "status": "ready" if calibrations else "no_applicable_fossils",
        "tree": file_record(source), "samples": file_record(samples), "taxonomy": file_record(taxonomy),
        "dataset": {"version": "v1.0", "source": dataset.source, "sha256": dataset.sha256},
        "nwkit": nwkit.__version__, "map": file_record(mapping) if mapping else None,
        "calibrated_nodes": len(calibrations), "fossil_status": dict(Counter(row["status"] for row in rows)),
        "exclusion_reasons": dict(Counter(row["reason"] for row in rows if row["status"] == "excluded")),
    })


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["tree", "samples", "taxonomy", "outdir"]:
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--mapping")
    parser.add_argument("--workbook")
    prepare(**vars(parser.parse_args()))
