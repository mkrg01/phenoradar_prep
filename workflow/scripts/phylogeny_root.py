#!/usr/bin/env python3
"""Prepare ASTRAL anchors and root inferred trees."""
import argparse
import json
import sqlite3
import tempfile
from pathlib import Path
from types import SimpleNamespace

from common import atomic_writer, file_record, now, read_tsv, write_json, write_tsv
from sample_identity import species_id
from phylogeny_outgroup import outgroup_ids, root_on_outgroup, validate_outgroup, validate_root
from nwkit_runtime import nwkit_backend


def species_rows(samples):
    result = {}
    for row in read_tsv(samples):
        name = row["species"]
        if name in result and any(row[k] != result[name][k] for k in ["taxid", "cds"]):
            raise ValueError(f"conflicting species inputs: {name}")
        result[name] = row
    if not result:
        raise ValueError("no selected species")
    return dict(sorted(result.items()))


def ncbi_guide(rows, taxonomy_db):
    """Build an NWKIT guide directly from the frozen database, without cache discovery."""
    nwkit_backend()
    import pandas as pd
    from ete4 import NCBITaxa
    from ete4.ncbi_taxonomy.ncbiquery import DB_VERSION
    from nwkit.constrain import get_lineages_from_taxid, get_taxid_counts, taxid2tree
    from nwkit.rooting_state import set_rooting_info
    from nwkit.util import remove_singleton

    database = Path(taxonomy_db).resolve()
    with sqlite3.connect(database.as_uri() + "?mode=ro", uri=True) as db:
        if db.execute("SELECT version FROM stats").fetchone()[0] != DB_VERSION:
            raise ValueError("taxonomy snapshot schema differs from ETE; prepare a compatible snapshot")
    groups = {}
    for name, row in rows.items():
        groups.setdefault(row["taxid"], []).append(name)
    unique = pd.DataFrame([{"leaf_name": names[0], "taxid": tid} for tid, names in groups.items()])
    ncbi = NCBITaxa(dbfile=str(database), update=False)
    try:
        lineages = get_lineages_from_taxid(unique, rank="no", ncbi=ncbi)
        tree = taxid2tree(lineages, get_taxid_counts(lineages), ncbi=ncbi)
    finally:
        ncbi.db.close()
    tree = remove_singleton(tree, verbose=False, preserve_branch_length=False)
    for leaf in list(tree.leaves()):
        names = groups[rows[leaf.name]["taxid"]]
        if len(names) > 1:
            leaf.name = ""
            for name in names:
                leaf.add_child(name=name)
    if set(tree.leaf_names()) != set(rows):
        raise ValueError("NCBI guide tree lost or added species; check taxids")
    set_rooting_info(tree, True, source="taxonomy")
    return tree


def ncbi_tree(samples, taxonomy_db, output, taxids):
    """Build the guide for ASTRAL's anchor from frozen NCBI lineages."""
    rows = species_rows(samples)
    tree = ncbi_guide(rows, taxonomy_db)
    write_tsv(taxids, ["leaf_name", "taxid"],
              [{"leaf_name": name, "taxid": row["taxid"]} for name, row in rows.items()])
    with atomic_writer(output) as handle:
        handle.write(tree.write(parser=9) + "\n")


def anchor_clade(tree, scores, taxa):
    """Choose the smallest basal lineage in the NCBI guide for ASTRAL."""
    while len(tree.children) == 1:
        tree = tree.children[0]
    children = [set(child.leaf_names()) for child in tree.children]
    if len(children) < 2:
        raise ValueError("NCBI root has fewer than two basal lineages")

    def rank(tips):
        return len({taxa[n] for n in tips}), min((-scores[n], n) for n in tips)

    detail = {"source": "ncbi", "basal_clade_sizes": sorted(map(len, children)),
              "basal_species_counts": sorted(len({taxa[n] for n in tips}) for tips in children)}
    return sorted(min(children, key=rank)), detail


def prepare_root(samples, metadata, output, qc, outgroup="auto", tree=None, for_anchor=False):
    if not (isinstance(outgroup, str) and outgroup.strip() or
            isinstance(outgroup, list) and outgroup and all(isinstance(n, str) and n.strip() for n in outgroup)):
        raise ValueError("outgroup must be auto, a species/sample ID or a nonempty list of IDs")
    rows = species_rows(samples)
    taxa = {n: species_id(r) for n, r in rows.items()}
    requested = outgroup
    if outgroup != "auto":
        requested_ids = [outgroup] if isinstance(outgroup, str) else outgroup
        selected = set()
        for value in requested_ids:
            taxon = taxa.get(value, value.replace(" ", "_"))
            found = {n for n in rows if taxa[n] == taxon}
            if not found:
                raise ValueError(f"outgroup is absent from inference species; no samples will be added: {value}")
            selected.update(found)
        outgroup = sorted(selected)
        detail = {"source": "explicit"}
    else:
        from ete4 import Tree
        scores = {r["species"]: float(r["busco_percent"]) for r in read_tsv(metadata)}
        guide = Tree(Path(tree).read_text(), parser=9)
        if set(guide.leaf_names()) != set(rows):
            raise ValueError("rooting guide species differ from selected dataset")
        outgroup, detail = anchor_clade(guide, scores, taxa)
        detail["guide"] = file_record(tree)
    outgroup = outgroup_ids(outgroup)
    validate_outgroup(rows, outgroup)
    labels = outgroup
    if for_anchor:
        anchor = min(outgroup, key=lambda n: (-scores[n], n)) if requested == "auto" else outgroup[0]
        labels = [anchor]
        detail["astral_root_anchor"] = anchor
    with atomic_writer(output) as handle:
        handle.write("\n".join(labels) + "\n")
    write_json(qc, {"outgroup": outgroup, "outgroup_species": sorted({taxa[n] for n in outgroup}),
               "requested": requested, **detail, "samples": file_record(samples),
                    "selection_scope": "inference_species", "candidate_species_count": len(set(taxa.values())), "candidate_sample_count": len(rows),
                    "created_at": now()})


def taxonomy_split(tree, rows, taxonomy_db, directory):
    """Use nwkit's default taxonomy sources to find the root split."""
    nwkit = nwkit_backend()
    from nwkit.root import DEFAULT_TAXONOMY_SOURCE_CHAIN, taxonomy_rooting, transfer_root

    Path(directory).mkdir(parents=True, exist_ok=True)
    sources = DEFAULT_TAXONOMY_SOURCE_CHAIN.split(",")
    with tempfile.TemporaryDirectory(prefix=".rooting-", dir=directory) as tmp:
        names = Path(tmp) / "species_map.tsv"
        write_tsv(names, ["leaf_name", "species_label", "taxonomy_query"],
                  [{"leaf_name": name, "species_label": species_id(row),
                    "taxonomy_query": row.get("scientific_name") or species_id(row).replace("_", " ")}
                   for name, row in rows.items()])
        errors = []
        for source in sources:
            try:
                if source == "ncbi":
                    rooted = transfer_root(tree.copy(), ncbi_guide(rows, taxonomy_db), verbose=True,
                                           redistribute_root_length=False)
                else:
                    rooted = taxonomy_rooting(tree.copy(), taxonomy_source=source, verbose=True,
                                              args=SimpleNamespace(download_dir=tmp, species_map_tsv=str(names)))
                break
            except ValueError as exc:
                errors.append(f"{source}: {exc}")
        else:
            raise ValueError("Automatic rooting failed. Specify phylogeny.outgroup explicitly. "
                             + " | ".join(errors))

    outgroup = min((sorted(child.leaf_names()) for child in rooted.children), key=len)
    detail = {"source": "nwkit", "method": "nwkit_taxonomy", "nwkit_version": nwkit.__version__,
              "taxonomy_sources": sources, "taxonomy_source": source,
              "taxonomy_db": file_record(taxonomy_db)}
    return outgroup, detail


def root_species_tree(samples, tree, tree_qc, anchor_qc, output, qc, outgroup_file,
                      rooting_qc, taxonomy_db=None):
    from ete4.parser.newick import make_parser
    from infer_phylogeny import read_tree

    rows = species_rows(samples)
    inferred = read_tree(tree, rows)
    selection = json.loads(Path(anchor_qc).read_text())
    report = json.loads(Path(tree_qc).read_text())
    anchor = selection["astral_root_anchor"]
    if report.get("astral_root_anchor") != anchor:
        raise ValueError("ASTRAL tree and anchor selection disagree")
    validate_root(inferred, [anchor])
    if selection["requested"] == "auto":
        outgroup, detail = taxonomy_split(inferred, rows, taxonomy_db, Path(output).parent)
    else:
        outgroup = outgroup_ids(selection["outgroup"])
        detail = {"source": "explicit", "method": "exact_outgroup"}

    # Apply the split to the original with numerical supports as edge attributes.
    for node in inferred.traverse():
        if not node.is_leaf and not node.is_root and node.name:
            node.support = float(node.name)
            node.name = ""
    root_on_outgroup(inferred, outgroup)
    for node in inferred.traverse():
        if not node.is_leaf:
            node.name = "" if node.is_root or node.support is None else str(node.support)
    with atomic_writer(output) as handle:
        handle.write(inferred.write(parser=make_parser(1, dist="%.17g")) + "\n")
    with atomic_writer(outgroup_file) as handle:
        handle.write("\n".join(outgroup) + "\n")
    rooting = {"created_at": now(), **detail, "outgroup": outgroup,
               "outgroup_species": sorted({species_id(rows[n]) for n in outgroup}),
               "requested": selection["requested"], "astral_root_anchor": anchor,
               "anchor_selection": file_record(anchor_qc), "input_tree": file_record(tree)}
    write_json(rooting_qc, rooting)
    report.update(created_at=now(), outgroup=outgroup, rooting=file_record(rooting_qc),
                  root_edge="final outgroup split",
                  multi_sample_basal_branch_length_validation=(
                      "not independently benchmarked" if len(outgroup) > 1 else "not applicable"))
    write_json(qc, report)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    p = sub.add_parser("ncbi_tree")
    for name in ["samples", "taxonomy-db", "output", "taxids"]:
        p.add_argument("--" + name, required=True)
    p = sub.add_parser("prepare_root")
    for name in ["samples", "metadata", "output", "qc"]:
        p.add_argument("--" + name, required=True)
    p.add_argument("--outgroup", default="auto")
    p.add_argument("--outgroup-json", dest="outgroup", type=json.loads, default=argparse.SUPPRESS)
    p.add_argument("--tree")
    p.add_argument("--for-anchor", action="store_true")
    p = sub.add_parser("root_species_tree")
    for name in ["samples", "tree", "tree-qc", "anchor-qc", "output", "qc", "outgroup-file", "rooting-qc"]:
        p.add_argument("--" + name, required=True)
    p.add_argument("--taxonomy-db")
    args = vars(parser.parse_args())
    globals()[args.pop("action")](**args)
