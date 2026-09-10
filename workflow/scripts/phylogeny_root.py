#!/usr/bin/env python3
"""NCBI guide trees and conservative single-outgroup selection for both runs."""
import argparse
import tempfile
from pathlib import Path
from types import SimpleNamespace

from common import atomic_writer, file_record, now, read_tsv, write_json, write_tsv


def nwkit_backend():
    import nwkit
    if nwkit.__version__ != "0.27.0":
        raise ValueError("use the pinned nwkit 0.27.0 in workflow/envs/timetree.yaml")
    return nwkit


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


def ncbi_tree(samples, taxonomy_db, output, taxids):
    """Run nwkit constrain against the existing SQLite snapshot, without downloads.

    In pinned nwkit 0.27.0, constrain uses lineage queries, not ETE's traversal
    pickle. An isolated cache link supplies the frozen DB to its standard API.
    """
    nwkit_backend()
    from ete4.ncbi_taxonomy.ncbiquery import DB_VERSION
    import sqlite3
    from nwkit.constrain import constrain_main
    database = Path(taxonomy_db).resolve()
    with sqlite3.connect(database.as_uri() + "?mode=ro", uri=True) as db:
        if db.execute("SELECT version FROM stats").fetchone()[0] != DB_VERSION:
            raise ValueError("taxonomy snapshot schema differs from ETE; prepare a compatible snapshot")
    rows = species_rows(samples)
    write_tsv(taxids, ["leaf_name", "taxid"],
              [{"leaf_name": name, "taxid": row["taxid"]} for name, row in rows.items()])
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".ncbi-", dir=Path(output).parent) as tmp:
        cache = Path(tmp) / "ete4"
        cache.mkdir()
        (cache / "taxa.sqlite").symlink_to(database)
        args = SimpleNamespace(taxid_tsv=str(taxids), species_list=None, backbone="ncbi", rank="no",
                               download_dir=tmp, outfile=str(Path(tmp) / "tree.nwk"), collapse=False,
                               outformat=9, quoted_node_names=False)
        constrain_main(args)
        from ete4 import Tree
        tree = Tree(Path(args.outfile).read_text(), parser=9)
        if set(tree.leaf_names()) != set(rows):
            raise ValueError("NCBI guide tree lost or added species; check taxids")
        with atomic_writer(output) as handle:
            handle.write(tree.write(parser=9) + "\n")


def apg_reference(representatives, taxids, taxonomy_db):
    """Expand nwkit's bundled APG IV order tree using exact frozen NCBI taxids."""
    from importlib import resources
    from ete4 import Tree, NCBITaxa
    from nwkit.constrain import (initialize_tree, delete_nomatch_leaves,
                                polytomize_one2many_matches)
    from nwkit.util import remove_singleton
    source = resources.files("nwkit").joinpath("data_tree/apgiv.nwk")
    tree = initialize_tree(Tree(source.read_text(), parser=0))
    by_order = {leaf.name: leaf for leaf in tree.leaves()}
    ncbi = NCBITaxa(dbfile=str(Path(taxonomy_db).resolve()), update=False)
    try:
        for name in representatives:
            lineage_names = set(ncbi.get_taxid_translator(ncbi.get_lineage(int(taxids[name]))).values())
            matches = set(by_order) & lineage_names
            if len(matches) != 1:
                raise ValueError(f"APG IV cannot resolve basal representative {name}; specify phylogeny.outgroup")
            leaf = by_order[matches.pop()]
            leaf.props["has_taxon"] = True
            leaf.props["taxon_names"].append(name)
    finally:
        ncbi.db.close()
    tree = delete_nomatch_leaves(tree)
    tree = polytomize_one2many_matches(tree)
    tree = remove_singleton(tree, verbose=False, preserve_branch_length=False)
    return tree, file_record(source)


def resolve_outgroup(tree, scores, reference=None):
    """Accept only a singleton basal lineage, never one tip of a larger clade.

    If the NCBI root is unresolved, a reference callback receives one species
    per basal lineage. The production callback uses the bundled APG IV tree.
    """
    nwkit_backend()
    while len(tree.children) == 1:
        tree = tree.children[0]
    children = [set(child.leaf_names()) for child in tree.children]
    if len(children) < 2:
        raise ValueError("at least two basal lineages are required for automatic rooting")
    singletons = {next(iter(names)) for names in children if len(names) == 1}
    detail = {"source": "ncbi", "basal_clade_sizes": sorted(map(len, children))}
    if len(children) == 2:
        if not singletons:
            raise ValueError("NCBI root has two multi-species clades; specify phylogeny.outgroup explicitly")
        return min(singletons, key=lambda n: (-scores[n], n)), detail
    if not singletons:
        raise ValueError("no singleton basal lineage; specify phylogeny.outgroup explicitly")
    representatives = sorted(min(names, key=lambda n: (-scores[n], n)) for names in children)
    if reference is None:
        raise ValueError("NCBI root is unresolved; a rooted reference or explicit phylogeny.outgroup is required")
    guide, source = reference(representatives)
    while len(guide.children) == 1:
        guide = guide.children[0]
    if set(guide.leaf_names()) != set(representatives) or len(guide.children) != 2:
        raise ValueError("reference root is unresolved or missing species; specify phylogeny.outgroup")
    candidates = [c.name for c in guide.children if c.is_leaf and c.name in singletons]
    if len(candidates) != 1:
        raise ValueError("reference does not define a single-species outgroup; specify phylogeny.outgroup")
    detail.update(source="nwkit_apgiv", reference=source, representatives=representatives,
                  reference_tree=guide.write(parser=9))
    return candidates[0], detail


def prepare_root(samples, metadata, output, qc, outgroup="auto", tree=None,
                 taxonomy_db=None):
    if not isinstance(outgroup, str) or not outgroup.strip():
        raise ValueError("outgroup must be auto or an exact inference species label; use auto instead of null")
    rows = species_rows(samples)
    if outgroup != "auto":
        if outgroup not in rows:
            raise ValueError(f"outgroup is absent from inference species; no species will be added: {outgroup}")
        detail = {"source": "explicit"}
    else:
        from ete4 import Tree
        scores = {r["species"]: float(r["busco_percent"]) for r in read_tsv(metadata)}
        guide = Tree(Path(tree).read_text(), parser=9)
        if set(guide.leaf_names()) != set(rows):
            raise ValueError("rooting guide species differ from selected dataset")
        outgroup, detail = resolve_outgroup(guide, scores,
            reference=lambda names: apg_reference(names, {n: r["taxid"] for n, r in rows.items()}, taxonomy_db))
        detail["guide"] = file_record(tree)
    with atomic_writer(output) as handle:
        handle.write(outgroup + "\n")
    write_json(qc, {"outgroup": outgroup, **detail, "samples": file_record(samples),
                    "selection_scope": "inference_species", "candidate_species_count": len(rows),
                    "created_at": now()})


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
    p.add_argument("--tree")
    p.add_argument("--taxonomy-db")
    args = vars(parser.parse_args())
    globals()[args.pop("action")](**args)
