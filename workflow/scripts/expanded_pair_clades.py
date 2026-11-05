#!/usr/bin/env python3
"""Expand fixed contrast pairs to maximal clades containing no other pair."""
import argparse
import json
import os
from pathlib import Path
import shutil
import tempfile

from common import file_record, now, read_tsv, sha256, write_json, write_tsv

COLUMN = "expanded_pair_clade_id"
POLICY = "maximal rooted clade containing one selected contrast pair and no other selected pair"
TABLE_FIELDS = [COLUMN, "contrast_pair_id", "state_a", "state_b", "n_representatives",
                "n_species", "n_samples", "n_species_a", "n_species_b", "n_samples_a", "n_samples_b"]


def clade_roots(tree, rows):
    """Return disjoint maximal clades; IDs retain their anchoring pair IDs."""
    tips = list(tree.leaves())
    if len({tip.name for tip in tips}) != len(tips) or {tip.name for tip in tips} != set(rows):
        raise ValueError("expanded-pair labels differ from summary-tree tips")
    pairs = {}
    for tip in tips:
        pair = str(rows[tip.name].get("contrast_pair_id") or "")
        if pair:
            pairs.setdefault(pair, []).append(tip)
    pair_sets = {}
    for node in tree.traverse("postorder"):
        pair_sets[node] = (set().union(*(pair_sets[c] for c in node.children)) if node.children
                           else ({str(rows[node.name]["contrast_pair_id"])} if rows[node.name].get("contrast_pair_id") else set()))
    roots = {}
    for pair, members in sorted(pairs.items()):
        if len(members) != 2:
            raise ValueError(f"contrast pair {pair} must have exactly two representative tips")
        node = tree.common_ancestor(members)
        if set(node.leaf_names()) != {tip.name for tip in members}:
            raise ValueError(f"contrast pair {pair} is not a summary-tree clade")
        while node.up is not None and pair_sets[node.up] == {pair}:
            node = node.up
        roots[pair] = node
    return roots


def annotate(tree, metadata, trait):
    """Compose summary-tip clades with the existing full sample membership."""
    representatives = {r["species"]: r for r in metadata if str(r["is_representative"]) == "1"}
    if len({r["species"] for r in metadata}) != len(metadata):
        raise ValueError("duplicate species in contrast metadata")
    roots = clade_roots(tree, representatives) if tree is not None else {}
    if tree is None and representatives:
        raise ValueError("empty summary tree has nonempty representatives")
    membership = {}
    for pair, root in roots.items():
        for tip in root.leaves():
            if tip.name in membership:
                raise ValueError("expanded pair clades overlap")
            membership[tip.name] = pair
    result = []
    for row in metadata:
        representative = row.get("representative", "")
        if representative and representative not in representatives:
            raise ValueError(f"unknown contrast representative: {representative}")
        expanded = membership.get(representative, "")
        if row.get("contrast_pair_id") and expanded != str(row["contrast_pair_id"]):
            raise ValueError("expanded clade must retain the original contrast pair")
        if expanded and not row[trait]:
            raise ValueError("missing-trait species cannot inherit an expanded pair clade")
        result.append({**row, COLUMN: expanded})
    clades = []
    for pair, root in roots.items():
        members = [r for r in result if r[COLUMN] == pair]
        states = sorted({r[trait] for r in representatives.values() if str(r.get("contrast_pair_id") or "") == pair})
        if len(states) != 2:
            raise ValueError(f"contrast pair {pair} requires two opposite trait states")
        def species(rows):
            return len({r.get("species_id") or r["species"] for r in rows})
        a, b = ([r for r in members if r[trait] == state] for state in states)
        clades.append({COLUMN: pair, "contrast_pair_id": pair, "state_a": states[0], "state_b": states[1],
                       "n_representatives": len(list(root.leaves())), "n_species": species(members),
                       "n_samples": len(members), "n_species_a": species(a), "n_species_b": species(b),
                       "n_samples_a": len(a), "n_samples_b": len(b)})
    assigned = [r for r in result if r[COLUMN]]
    report = {"expanded_pair_clades": len(clades),
              "expanded_pair_species": len({r.get("species_id") or r["species"] for r in assigned}),
              "expanded_pair_samples": len(assigned), "expanded_pair_policy": POLICY}
    return result, clades, report


def validate_result(folder, trait):
    """Recompute memberships and counts for native-output integrity checks."""
    from ete4 import Tree
    folder = Path(folder)
    text = (folder / "summary_tree.nwk").read_text()
    saved = read_tsv(folder / "species_metadata.tsv")
    rows, clades, report = annotate(Tree(text, parser=0) if text.strip() else None, saved, trait)
    expected = [{key: str(value) for key, value in row.items()} for row in clades]
    summary = json.loads((folder / "summary.json").read_text())
    if rows != saved or expected != read_tsv(folder / "expanded_pair_clades.tsv") or any(summary.get(k) != v for k, v in report.items()):
        raise ValueError("expanded-pair outputs differ from the saved summary tree")
    pairs = read_tsv(folder / "contrast_pairs.tsv")
    by_name = {r["species"]: r for r in saved}
    if len(pairs) != len(clades) or {r["contrast_pair_id"] for r in pairs} != {r[COLUMN] for r in clades}:
        raise ValueError("expanded-pair anchors differ from saved contrast pairs")
    for pair in pairs:
        for key in ["representative_a", "representative_b"]:
            if by_name.get(pair[key], {}).get("contrast_pair_id") != pair["contrast_pair_id"]:
                raise ValueError("expanded-pair representative differs from saved contrast pairs")
    return report


def export_analysis(analysis, outdir=None):
    """Publish a derivative of saved pairs, without changing a frozen analysis."""
    from ete4 import Tree
    from consumer_metadata import BRANCHES, create, validate
    from plot_contrast_tree import plot
    source = Path(analysis).resolve()
    requested = Path(outdir) if outdir else source / "expanded_pair_clades"
    if requested.is_symlink():
        raise ValueError("expanded-pair output must not be a symlink")
    out = requested.resolve()
    if out.exists():
        raise ValueError("expanded-pair output already exists; choose a new directory")
    inputs = {}
    def consume(path):
        path = Path(path).resolve()
        if path == out or path.is_relative_to(out):
            raise ValueError("expanded-pair output overlaps source inputs")
        inputs[path] = file_record(path)
        return path
    manifest = json.loads(consume(source / "analysis.json").read_text())
    trait = manifest["pipeline"]["trait"]
    phy = manifest["pipeline"]["phylogeny"]
    branches = phy.get("trees", []) if phy.get("contrast_pairs", {}).get("enabled") else []
    if not branches or set(branches) - BRANCHES or len(branches) != len(set(branches)):
        raise ValueError("analysis requires configured contrast branches")
    if (source / "completed.json").is_file():
        consume(source / "completed.json")
    samples = consume(source / "metadata/samples.tsv")
    base = consume(source / "metadata/species_metadata.tsv")
    out.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=".expanded-pairs-", dir=out.parent))
    reports, pair_tables = {}, {}
    try:
        (stage / "metadata").mkdir()
        shutil.copyfile(base, stage / "metadata/species_metadata_base.tsv")
        shutil.copyfile(samples, stage / "metadata/samples.tsv")
        for branch in branches:
            original = source / "phylogeny" / branch / "contrast"
            target = stage / "phylogeny" / branch / "contrast"
            target.mkdir(parents=True)
            text = consume(original / "summary_tree.nwk").read_text()
            tree = Tree(text, parser=0) if text.strip() else None
            rows, clades, report = annotate(tree, read_tsv(consume(original / "species_metadata.tsv")), trait)
            summary = json.loads(consume(original / "summary.json").read_text())
            if summary["trait"] != trait or summary["contrast_pairs"] != len(clades):
                raise ValueError("saved contrast summary differs from expanded-pair input")
            for name in ["summary_tree.nwk", "contrast_pairs.tsv"]:
                shutil.copyfile(consume(original / name), target / name)
            write_tsv(target / "species_metadata.tsv", list(rows[0]) if rows else ["species", trait, COLUMN], rows)
            write_tsv(target / "expanded_pair_clades.tsv", TABLE_FIELDS, clades)
            write_json(target / "summary.json", {**summary, **report})
            validate_result(target, trait)
            plot(target / "summary_tree.nwk", target / "species_metadata.tsv", target / "summary.json",
                 target, samples=stage / "metadata/samples.tsv", expanded=True)
            pair_tables[branch] = target / "species_metadata.tsv"
            reports[branch] = report
        create(stage / "metadata/species_metadata_base.tsv", stage / "metadata/species_metadata_base.tsv",
               stage / "metadata/species_metadata.tsv", pair_tables, trait)
        validate(stage / "metadata/species_metadata.tsv", trait)
        for path, record in inputs.items():
            if sha256(path) != record["sha256"]:
                raise ValueError(f"expanded-pair input changed during export: {path}")
        code = Path(__file__).resolve().parent
        outputs = [dict(file_record(p), path=str(p.relative_to(stage))) for p in sorted(stage.rglob("*")) if p.is_file()]
        write_json(stage / "provenance.json", {
            "schema_version": 1, "created_at": now(), "source_analysis": str(source),
            "version": (code.parent.parent / "VERSION").read_text().strip(),
            "reestimated": False, "pairs_reselected": False, "branches": reports,
            "sources": list(inputs.values()), "outputs": outputs,
            "implementation": [file_record(code / name) for name in [
                "expanded_pair_clades.py", "plot_contrast_tree.py", "consumer_metadata.py", "phenoradar_metadata.py", "common.py"]],
        })
        os.rename(stage, out)
        return {"output": str(out), "branches": reports}
    finally:
        if stage.exists():
            shutil.rmtree(stage)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis", required=True, help="Analysis with completed, saved contrast outputs")
    parser.add_argument("--outdir", help="New derivative directory; default: ANALYSIS/expanded_pair_clades")
    args = parser.parse_args()
    print(json.dumps(export_analysis(**vars(args)), indent=2))
