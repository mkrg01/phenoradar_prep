"""Pruned derivatives of completed all-species phylogenies, without inference."""
from collections import Counter
import json
import math

from common import write_json, write_tsv
from layout import PHYLOGENY_BRANCHES

ALL_PHYLOGENY = PHYLOGENY_BRANCHES["all"]


def read_tree(text, allowed, exact=False):
    from ete4 import Tree
    if text.count(";") != 1:
        raise ValueError("expected one Newick tree")
    tree = Tree(text.strip(), parser=1)
    names = list(tree.leaf_names())
    if not names or len(names) != len(set(names)) or set(names) - allowed:
        raise ValueError("duplicate or unexpected tree species")
    if exact and set(names) != allowed:
        raise ValueError("species tree differs from source sample manifest")
    if any(n.dist is None or not math.isfinite(n.dist) or n.dist < 0
           for n in tree.traverse() if not n.is_root):
        raise ValueError("tree requires finite nonnegative branch lengths")
    return tree


def prune(tree, keep):
    """Keep path lengths; supports cannot be reinterpreted on contracted edges."""
    selected = set(tree.leaf_names()) & keep
    if len(selected) < 2:
        return None
    tree.prune(sorted(selected), preserve_branch_length=True)
    for node in tree.traverse():
        if not node.is_leaf:
            node.name = ""
    # Parser 5 writes leaf names and all branch lengths, with no internal labels.
    from ete4.parser.newick import make_parser
    return tree.write(parser=make_parser(5, dist="%.17g"), format_root_node=False) + "\n"


def filter_fasta(source, destination, keep, allowed):
    """Filter species-labeled marker FASTA without changing columns or residues."""
    from filter_species import fasta_records, safe_name
    seen, retained, width = set(), [], None
    for header, sequence in fasta_records(source):
        name = header.split()[0]
        safe_name(name)
        if name in seen or name not in allowed:
            raise ValueError(f"unexpected/duplicate species in marker FASTA: {source}: {name}")
        seen.add(name)
        if width is None:
            width = len(sequence)
        if width != len(sequence):
            raise ValueError(f"unequal alignment lengths: {source}")
        if name in keep:
            retained.append((header, sequence))
    if not seen:
        raise ValueError(f"empty marker alignment: {source}")
    if retained:
        destination.parent.mkdir(parents=True, exist_ok=True)
        with destination.open("w") as handle:
            for header, sequence in retained:
                handle.write(f">{header}\n{sequence}\n")
    return len(seen), len(retained), width


def export_phylogeny(job):
    from filter_species import safe_name
    source, out = job.source / ALL_PHYLOGENY, job.stage / ALL_PHYLOGENY
    out.mkdir(parents=True)
    original = json.loads(job.input(source / "species_tree.json").read_text())
    tree = read_tree(job.input(source / "species_tree.nwk").read_text(), job.species, exact=True)
    if original.get("species") != len(job.species):
        raise ValueError("species-tree QC differs from sample manifest")
    genes = job.input(source / "gene_trees.nwk")
    expected_hash = original.get("input", {}).get("sha256")
    if expected_hash and expected_hash != job.records[str(genes)]["sha256"]:
        raise ValueError("gene trees differ from the species-tree inference input")
    report = json.loads(job.input(source / "gene_trees.json").read_text())
    loci = report.get("retained", [])
    if report.get("status") != "retained" or not loci:
        raise ValueError("source gene-tree QC is not a completed inference input")
    markers = [safe_name(r["marker"]) for r in loci]
    if len(markers) != len(set(markers)):
        raise ValueError("duplicate gene-tree markers")
    source_coverage, coverage, details = Counter(), Counter(), []
    outgroup = original.get("outgroup")
    rooting = "original_outgroup_retained" if outgroup in job.keep else "original_root_unverified"
    if outgroup and outgroup not in job.keep:
        rooting = "original_outgroup_removed_requires_review"
    tree_text = prune(tree, job.keep)
    if tree_text:
        (out / "species_tree.pruned.nwk").write_text(tree_text)
    with genes.open() as handle, (out / "gene_trees.pruned.nwk").open("w") as merged:
        used = 0
        for line in handle:
            if not line.strip():
                continue
            if used >= len(markers):
                raise ValueError("more gene trees than QC records")
            marker = markers[used]
            gene = read_tree(line, job.species)
            present = set(gene.leaf_names())
            source_coverage.update(present)
            retained = present & job.keep
            text = prune(gene, job.keep)
            path = f"gene_trees/{marker}.pruned.nwk"
            if text:
                (out / "gene_trees").mkdir(exist_ok=True)
                (out / path).write_text(text)
                merged.write(text)
                coverage.update(retained)
            details.append({"marker": marker, "source_species": len(present), "retained_species": len(retained),
                            "status": "pruned" if text else "fewer_than_two_tips", "tree": path if text else "",
                            "at_least_four_tips": len(retained) >= 4})
            used += 1
    if used != len(markers) or set(source_coverage) != job.species:
        raise ValueError("gene-tree coverage/count differs from completed source inference")
    for marker in markers:
        for folder in ["alignments", "alignments/raw"]:
            path = source / folder / f"{marker}.faa"
            if path.is_file():
                before, after, columns = filter_fasta(job.input(path), out / folder / path.name,
                                                       job.keep, job.species)
                job.counts[f"{ALL_PHYLOGENY}/{folder}/{path.name}"] = dict(before=before, after=after, columns=columns)
                column_map = path.with_suffix(".columns.tsv")
                if column_map.is_file() and after:
                    job.copy(column_map, out / folder / column_map.name)
    write_tsv(out / "gene_trees.tsv", ["marker", "source_species", "retained_species", "status", "tree", "at_least_four_tips"], details)
    write_tsv(out / "species_coverage.tsv", ["species", "gene_trees", "represented"],
              [dict(species=s, gene_trees=coverage[s], represented=coverage[s] > 0) for s in sorted(job.keep)])
    # Old node IDs, calibration bounds and fit statistics do not describe a new
    # inference. Export only a pruned dated tree and clearly marked provenance.
    dated_source = source / "dating/species_tree.dated.nwk"
    dated = False
    if dated_source.is_file():
        dated_tree = read_tree(job.input(dated_source).read_text(), job.species, exact=True)
        dated_text = prune(dated_tree, job.keep)
        if dated_text:
            (out / "dating").mkdir()
            (out / "dating/species_tree.dated.pruned.nwk").write_text(dated_text)
            dated = True
    write_json(out / "pruning.json", {
        "operation": "prune_only", "reestimated": False, "source_species": len(job.species),
        "retained_species": len(job.keep), "species_tree_written": bool(tree_text),
        "branch_length_unit": original.get("branch_length_unit", "unspecified"),
        "branch_lengths": "sums along surviving paths; not reestimated",
        "internal_labels": "removed; source branch supports are not recomputed supports",
        "source_outgroup": outgroup, "rooting": rooting, "dated_tree_written": dated,
        "dated_tree_note": "Source time-tree path lengths retained; calibrations and node ages were not refitted",
        "gene_trees": details,
        "limitations": ["Original inference can still reflect excluded sequences.",
                        "Pruned trees are not new ASTRAL/CASTLES-II/LSD2 estimates.",
                        "Trees with fewer than four tips are not suitable ASTRAL inputs.",
                        "No old calibration, support, node-age, or inference-QC tables are relabeled as filtered results."]})
    job.counts["phylogeny"] = dict(before=len(job.species), after=len(job.keep), gene_trees=used)
