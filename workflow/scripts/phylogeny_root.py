#!/usr/bin/env python3
"""NCBI guide trees and conservative outgroup-species selection for both runs."""
import argparse
import json
import tempfile
from pathlib import Path
from types import SimpleNamespace

from common import atomic_writer, file_record, now, read_tsv, write_json, write_tsv
from sample_identity import species_id
from phylogeny_outgroup import outgroup_ids, root_on_outgroup, validate_outgroup, validate_root


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
    groups = {}
    for name, row in rows.items():
        groups.setdefault(row["taxid"], []).append(name)
    write_tsv(taxids, ["leaf_name", "taxid"],
              [{"leaf_name": name, "taxid": row["taxid"]} for name, row in rows.items()])
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".ncbi-", dir=Path(output).parent) as tmp:
        cache = Path(tmp) / "ete4"
        cache.mkdir()
        (cache / "taxa.sqlite").symlink_to(database)
        unique = Path(tmp) / "unique_taxids.tsv"
        write_tsv(unique, ["leaf_name", "taxid"], [{"leaf_name": names[0], "taxid": tid} for tid, names in groups.items()])
        args = SimpleNamespace(taxid_tsv=str(unique), species_list=None, backbone="ncbi", rank="no",
                               download_dir=tmp, outfile=str(Path(tmp) / "tree.nwk"), collapse=False,
                               outformat=9, quoted_node_names=False)
        constrain_main(args)
        from ete4 import Tree
        tree = Tree(Path(args.outfile).read_text(), parser=9)
        for leaf in list(tree.leaves()):
            names = groups[rows[leaf.name]["taxid"]]
            if len(names) > 1:
                leaf.name = ""
                for name in names:
                    leaf.add_child(name=name)
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


def resolve_outgroup(tree, scores, reference=None, taxa=None, for_anchor=False):
    """Select a basal lineage from the guide.

    Reference resolution receives one representative per basal lineage.
    Multi-species candidates are allowed when choosing an inference anchor.
    """
    nwkit_backend()
    while len(tree.children) == 1:
        tree = tree.children[0]
    names = set(tree.leaf_names())
    taxa = taxa or {n: n for n in names}
    children = [set(child.leaf_names()) for child in tree.children]
    if len(children) < 2:
        raise ValueError("NCBI root has fewer than two basal lineages")
    rank = lambda n: (-scores[n], n)
    singletons = {}
    for tips in children:
        species = {taxa[n] for n in tips}
        if len(species) == 1 and tips == {n for n in names if taxa[n] in species}:
            singletons[min(tips, key=rank)] = sorted(tips)
    detail = {"source": "ncbi", "basal_clade_sizes": sorted(map(len, children)),
              "basal_species_counts": sorted(len({taxa[n] for n in tips}) for tips in children)}
    candidate = lambda tips: (len({taxa[n] for n in tips}), rank(min(tips, key=rank)))
    if len(children) == 2:
        if singletons:
            return singletons[min(singletons, key=rank)], detail
        if for_anchor:
            return sorted(min(children, key=candidate)), detail
        raise ValueError("NCBI root has two multi-species clades; specify phylogeny.outgroup explicitly")
    if not singletons and not for_anchor:
        raise ValueError("no complete single-species basal lineage; specify phylogeny.outgroup explicitly")
    lineages = {min(tips, key=rank): sorted(tips) for tips in children}
    representatives = sorted(lineages)
    if reference is None:
        raise ValueError("NCBI root is unresolved; a rooted reference or explicit phylogeny.outgroup is required")
    try:
        guide, source = reference(representatives)
        while len(guide.children) == 1:
            guide = guide.children[0]
        if set(guide.leaf_names()) != set(representatives) or len(guide.children) != 2:
            raise ValueError("reference root is unresolved or missing species; specify phylogeny.outgroup")
        candidates = [c.name for c in guide.children if c.is_leaf and (for_anchor or c.name in singletons)]
        if len(candidates) != 1:
            raise ValueError("reference does not define a single basal lineage; specify phylogeny.outgroup")
    except ValueError as exc:
        if not for_anchor:
            raise
        detail.update(source="ncbi_smallest_basal_candidate", reference_error=str(exc))
        return sorted(min(children, key=candidate)), detail
    detail.update(source="nwkit_apgiv", reference=source, representatives=representatives,
                  reference_tree=guide.write(parser=9))
    return (lineages if for_anchor else singletons)[candidates[0]], detail


def prepare_root(samples, metadata, output, qc, outgroup="auto", tree=None,
                 taxonomy_db=None, for_anchor=False):
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
        outgroup, detail = resolve_outgroup(guide, scores,
            reference=lambda names: apg_reference(names, {n: r["taxid"] for n, r in rows.items()}, taxonomy_db),
            taxa=taxa, for_anchor=for_anchor)
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


def taxonomy_split(tree, guide):
    """Use NWKIT's taxonomy root transfer to identify an exact molecular split."""
    from nwkit.root import transfer_root
    while len(guide.children) == 1:
        guide = guide.children[0]
    oriented = transfer_root(tree.copy(), guide)
    # NWKIT defines the smaller reference side as outgroup, taking the first on ties.
    outgroup = min((sorted(c.leaf_names()) for c in guide.children), key=len)
    validate_root(oriented, outgroup)
    return outgroup


def root_species_tree(samples, tree, tree_qc, anchor_qc, output, qc, outgroup_file,
                      rooting_qc, guide=None, taxonomy_db=None):
    from ete4 import Tree
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
        nwkit_backend()
        try:
            outgroup = taxonomy_split(inferred, Tree(Path(guide).read_text(), parser=9))
            detail = {"source": "ncbi", "reference": file_record(guide)}
        except ValueError as ncbi_error:
            try:
                reference, source = apg_reference(rows, {n: r["taxid"] for n, r in rows.items()}, taxonomy_db)
                outgroup = taxonomy_split(inferred, reference)
            except ValueError as apg_error:
                raise ValueError(f"cannot root the inferred tree by taxonomy; specify phylogeny.outgroup "
                                 f"(NCBI: {ncbi_error}; APG IV: {apg_error})") from apg_error
            detail = {"source": "nwkit_apgiv", "reference": source, "ncbi_error": str(ncbi_error)}
        detail["method"] = "nwkit_taxonomy"
    else:
        outgroup = outgroup_ids(selection["outgroup"])
        detail = {"source": "explicit", "method": "exact_outgroup"}

    # NWKIT 0.27 transfers a root on a copy. Apply the verified split to the
    # original with numerical supports treated as edge attributes by ETE.
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
    p.add_argument("--taxonomy-db")
    p.add_argument("--for-anchor", action="store_true")
    p = sub.add_parser("root_species_tree")
    for name in ["samples", "tree", "tree-qc", "anchor-qc", "output", "qc", "outgroup-file", "rooting-qc"]:
        p.add_argument("--" + name, required=True)
    p.add_argument("--guide")
    p.add_argument("--taxonomy-db")
    args = vars(parser.parse_args())
    globals()[args.pop("action")](**args)
