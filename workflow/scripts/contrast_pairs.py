#!/usr/bin/env python3
"""Trait-guided nwkit skims, representative manifests, pair tables and figures."""
import argparse
import json
import math
from collections import Counter
from pathlib import Path
from types import SimpleNamespace

from common import atomic_writer, file_record, now, read_tsv, write_json, write_tsv
from phylogeny_root import nwkit_backend, species_rows
from species_traits import read_species_traits


def skim(tree, rows, prefix, seed, contrastive=False, topology_only=False):
    """Use nwkit's grouping and selection functions, including for zero pairs."""
    nwkit_backend()
    import pandas as pd
    from nwkit.skim import (mark_traits_to_nodes, add_group_ids,
                            add_contrastive_clade_ids, sample_from_groups)
    args = SimpleNamespace(group_by="trait", filter_by="busco_percent", filter_mode="descending",
                           retain_per_clade=1, prioritize_non_missing=True, seed=seed)
    names = set(tree.leaf_names())
    frame = pd.DataFrame(sorted(rows, key=lambda r: r["leaf_name"]))
    if set(frame["leaf_name"]) != names or frame["leaf_name"].duplicated().any():
        raise ValueError("skim trait rows must match tree tips exactly")
    # No missing traits enter grouping.
    if frame["trait"].eq("").any():
        raise ValueError("skim requires observed traits")
    tree = mark_traits_to_nodes(tree, frame, args)
    frame = add_group_ids(frame, tree)
    if contrastive:
        frame = add_contrastive_clade_ids(frame, tree)
    sampled = sample_from_groups(frame, args)
    if contrastive:
        sampled = sampled.loc[sampled["contrastive_clade"].notna()]
    prefix = Path(prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    for suffix, data in [("all", frame), ("sampled", sampled)]:
        with atomic_writer(str(prefix) + f".{suffix}.tsv") as handle:
            data.sort_values(["group", "leaf_name"]).to_csv(handle, sep="\t", index=False)
    with atomic_writer(str(prefix) + ".nwk") as handle:
        if not sampled.empty:
            tree.prune(sampled["leaf_name"].tolist(), preserve_branch_length=True)
            handle.write(tree.write(parser=9 if topology_only else 0) + "\n")
        # There is no Newick representation of zero tips. An empty file means
        # no contrastive clades; summary.json and the TSV headers disambiguate it.
    return frame.to_dict("records"), sampled.to_dict("records")


def prepare(samples, metadata, traits, tree, outdir, trait="C4", seed=12345):
    from ete4 import Tree
    rows = species_rows(samples)
    annotation = read_species_traits(traits, trait)
    scores = {r["species"]: float(r["busco_percent"]) for r in read_tsv(metadata)}
    eligible = {n for n in rows if annotation.get(n, "") != ""}
    states = sorted({annotation[n] for n in eligible})
    if len(states) != 2:
        raise ValueError(f"contrast analysis requires exactly two observed trait states; found {states}")
    guide = Tree(Path(tree).read_text(), parser=9)
    if set(guide.leaf_names()) != set(rows):
        raise ValueError("NCBI guide differs from selected species")
    # NCBI distances are arbitrary and are never exported as evolutionary lengths.
    for node in guide.traverse():
        node.dist = 1
    guide.prune(sorted(eligible), preserve_branch_length=True)
    out = Path(outdir)
    skim_rows = [{"leaf_name": n, "trait": annotation[n], "busco_percent": scores[n]} for n in sorted(eligible)]
    all_rows, sampled = skim(guide, skim_rows, out / "ncbi_skim", seed, topology_only=True)
    representatives = {r["leaf_name"] for r in sampled}
    if len(representatives) < 4:
        raise ValueError("fewer than four representative species; BUSCO inference requires four")
    write_tsv(out / "samples.tsv", list(next(iter(rows.values()))), [rows[n] for n in sorted(representatives)])
    write_tsv(out / "traits.tsv", ["species", "trait", "busco_percent", "role"],
              [{"species": n, "trait": annotation.get(n, ""), "busco_percent": scores[n],
                "role": "observed" if n in eligible else "missing_trait"}
               for n in rows])
    write_json(out / "selection.json", {"created_at": now(), "species_trait": file_record(traits),
               "trait": trait, "states": states, "dataset_species": len(rows),
               "observed_species": len(eligible), "inference_species": len(representatives),
               "missing_trait_rows": sorted(set(rows) - set(annotation)), "seed": seed,
               "nwkit": "0.27.0", "ncbi_tree": file_record(tree)})


def rooted_tree(path, names, outgroup):
    """Validate the original rooted tree before any traits or exclusions prune it."""
    from ete4 import Tree
    text = Path(path).read_text()
    if text.count(";") != 1:
        raise ValueError("expected one inferred species tree")
    tree = Tree(text, parser=0)
    tips = list(tree.leaf_names())
    if len(tips) != len(set(tips)) or set(tips) != set(names):
        raise ValueError("inferred tree tips differ from species manifest")
    if outgroup not in names:
        raise ValueError("outgroup is absent from the species manifest")
    if len(tree.children) != 2 or not any(c.is_leaf and c.name == outgroup for c in tree.children):
        raise ValueError("inferred tree is not rooted with the selected outgroup")
    if any(n.dist is None or not math.isfinite(n.dist) or n.dist < 0
           for n in tree.traverse() if not n.is_root):
        raise ValueError("inferred tree requires finite nonnegative branch lengths")
    return tree


def summarize(tree, selection_dir, outgroup_file, outdir, seed=12345):
    selection, out = Path(selection_dir), Path(outdir)
    plan = json.loads((selection / "selection.json").read_text())
    all_traits = read_tsv(selection / "traits.tsv")
    first_all = read_tsv(selection / "ncbi_skim.all.tsv")
    first_reps = read_tsv(selection / "ncbi_skim.sampled.tsv")
    names = {r["leaf_name"] for r in first_reps}
    outgroup = Path(outgroup_file).read_text().strip()
    if outgroup not in names:
        raise ValueError("outgroup is absent from the selected representatives")
    inferred = rooted_tree(tree, names, outgroup)
    rows = [{"leaf_name": r["leaf_name"], "trait": r["trait"], "busco_percent": float(r["busco_percent"])}
            for r in first_reps]
    first_representatives = {r["group"]: r["leaf_name"] for r in first_reps}
    tip_for_species = {r["leaf_name"]: first_representatives[r["group"]] for r in first_all}
    return summarize_tree(inferred, rows, all_traits, tip_for_species, out, plan["trait"], seed, outgroup,
                          {"outgroup_file": file_record(outgroup_file), "inferred_tree": file_record(tree),
                           "selection": file_record(selection / "selection.json"),
                           "assignment": "non-representative species inherit their NCBI skim group membership"})


def summarize_tree(inferred, rows, all_traits, tip_for_species, outdir, trait, seed, outgroup, provenance):
    """Shared molecular skim, pair assignment and export for both input routes."""
    from ete4 import Tree
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    if rows:
        second_all, second_reps = skim(inferred, rows, out / "summary_tree", seed)
        summary_tree = Tree((out / "summary_tree.nwk").read_text(), parser=0)
        _, contrast_reps = skim(summary_tree, [
            {k: r[k] for k in ["leaf_name", "trait", "busco_percent"]} for r in second_reps
        ], out / "contrastive", seed, contrastive=True)
    else:
        second_all, second_reps, contrast_reps = [], [], []
        for prefix in ["summary_tree", "contrastive"]:
            fields = ["leaf_name", "trait", "busco_percent", "group"]
            if prefix == "contrastive":
                fields += ["contrastive_clade"]
            for suffix in ["all", "sampled"]:
                write_tsv(out / f"{prefix}.{suffix}.tsv", fields, [])
            with atomic_writer(out / f"{prefix}.nwk") as handle:
                handle.write("")
    # Map original species through their input-tree tips to molecular groups.
    second_group = {r["leaf_name"]: r["group"] for r in second_all}
    final_reps = {r["group"]: r["leaf_name"] for r in second_reps}
    members = {name: second_group[tip] for name, tip in tip_for_species.items()}
    counts = Counter(members.values())
    candidates = {}
    for row in contrast_reps:
        candidates.setdefault(int(row["contrastive_clade"]), []).append(row)
    pairs, paired_groups, unresolved = [], {}, []
    # Canonical ordering makes public pair IDs independent of Newick child order.
    for candidate in sorted(candidates.values(), key=lambda rs: sorted(r["leaf_name"] for r in rs)):
        if len(candidate) != 2 or len({r["trait"] for r in candidate}) != 2:
            unresolved.append(sorted(r["leaf_name"] for r in candidate))
            continue
        left, right = sorted(candidate, key=lambda r: str(r["trait"]))
        pair_id = len(pairs) + 1
        groups = [second_group[r["leaf_name"]] for r in [left, right]]
        for group in groups:
            if group in paired_groups:
                raise ValueError("one clade assigned to more than one contrast pair")
            paired_groups[group] = pair_id
        pairs.append({"contrast_pair_id": pair_id, "state_a": left["trait"], "state_b": right["trait"],
                      "representative_a": left["leaf_name"], "representative_b": right["leaf_name"],
                      "group_a": groups[0], "group_b": groups[1],
                      "n_species_a": counts[groups[0]], "n_species_b": counts[groups[1]]})
    metadata = []
    for row in all_traits:
        name = row["species"]
        group = members.get(name, "")
        metadata.append({"species": name, trait: row["trait"],
                         "role": "outgroup" if name == outgroup else row["role"],
                         "group": group, "representative": final_reps.get(group, ""),
                         "is_representative": int(final_reps.get(group) == name),
                         "n_species_in_group": counts.get(group, ""),
                         "contrast_pair_id": paired_groups.get(group, "")})
    fields = ["contrast_pair_id", "state_a", "state_b", "representative_a", "representative_b",
              "group_a", "group_b", "n_species_a", "n_species_b"]
    write_tsv(out / "contrast_pairs.tsv", fields, pairs)
    write_tsv(out / "species_metadata.tsv", ["species", trait, "role", "group", "representative",
              "is_representative", "n_species_in_group", "contrast_pair_id"], metadata)
    report = {"created_at": now(), "trait": trait, "nwkit": "0.27.0", "outgroup": outgroup,
               "summary_species": len(second_reps), "contrast_pairs": len(pairs), "seed": seed,
               "unresolved_contrastive_clades": unresolved, **provenance}
    write_json(out / "summary.json", report)
    return report


def from_tree(tree, tree_qc, samples, metadata, traits, outdir, trait="C4", seed=12345, exclude_species=()):
    """Assign pairs directly on a completed species tree, optionally pruning species."""
    if not isinstance(exclude_species, (list, tuple)) or any(not isinstance(n, str) for n in exclude_species):
        raise ValueError("exclude_species must be a list of species IDs")
    if len(set(exclude_species)) != len(exclude_species):
        raise ValueError("duplicate excluded species IDs")
    if not isinstance(trait, str) or not trait.strip() or trait in {
            "species", "role", "group", "representative", "is_representative", "n_species_in_group", "contrast_pair_id"}:
        raise ValueError("contrast.trait must name a non-reserved phenotype column")
    if type(seed) is not int or seed <= 0:
        raise ValueError("seed must be a positive integer")
    names = {r["species"] for r in read_tsv(samples)}
    qc = json.loads(Path(tree_qc).read_text())
    if qc.get("species") != len(names):
        raise ValueError("species-tree QC differs from species manifest")
    outgroup = qc.get("outgroup")
    inferred = rooted_tree(tree, names, outgroup)
    annotation = read_species_traits(traits, trait)
    retained = names - set(exclude_species)
    eligible = {n for n in retained if annotation.get(n, "") != ""}
    states = sorted({annotation[n] for n in eligible})
    if len(states) > 2:
        raise ValueError(f"contrast analysis supports at most two observed trait states; found {states}")
    scores = {}
    for row in read_tsv(metadata):
        name = row["species"]
        if name not in eligible:
            continue
        score = float(row["busco_percent"])
        if not math.isfinite(score) or not 0 <= score <= 100:
            raise ValueError(f"invalid BUSCO completeness: {name}")
        if name in scores and scores[name] != score:
            raise ValueError(f"conflicting BUSCO completeness: {name}")
        scores[name] = score
    if eligible - scores.keys():
        raise ValueError("missing BUSCO completeness: " + ", ".join(sorted(eligible - scores.keys())))
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    with atomic_writer(out / "observed_tree.nwk") as handle:
        if eligible:
            inferred.prune(sorted(eligible), preserve_branch_length=True)
            handle.write(inferred.write(parser=0) + "\n")
    rows = [{"leaf_name": n, "trait": annotation[n], "busco_percent": scores[n]} for n in sorted(eligible)]
    all_traits = [{"species": n, "trait": annotation.get(n, ""),
                   "role": "observed" if n in eligible else "missing_trait"} for n in sorted(retained)]
    return summarize_tree(inferred, rows, all_traits, {n: n for n in eligible}, out, trait, seed, outgroup,
                          {"mode": "inferred_tree", "inferred_tree": file_record(tree),
                           "tree_qc": file_record(tree_qc), "samples": file_record(samples),
                           "metadata": file_record(metadata), "species_trait": file_record(traits),
                           "source_species": len(names), "retained_species": len(retained),
                           "observed_species": len(eligible), "states": states,
                           "excluded_species": sorted(names - retained),
                           "missing_trait_species": sorted(retained - eligible),
                           "rooting": "source_root_inherited", "outgroup_in_observed_tree": outgroup in eligible,
                           "reestimated": False,
                           "assignment": "membership on the observed-species subtree of the input molecular tree"})



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    p = sub.add_parser("prepare")
    for name in ["samples", "metadata", "traits", "tree", "outdir"]:
        p.add_argument("--" + name, required=True)
    p.add_argument("--trait", default="C4")
    p.add_argument("--seed", type=int, default=12345)
    p = sub.add_parser("summarize")
    for name in ["tree", "selection-dir", "outgroup-file", "outdir"]:
        p.add_argument("--" + name, required=True)
    p.add_argument("--seed", type=int, default=12345)
    p = sub.add_parser("from_tree")
    for name in ["tree", "tree-qc", "samples", "metadata", "traits", "outdir"]:
        p.add_argument("--" + name, required=True)
    p.add_argument("--trait", default="C4")
    p.add_argument("--seed", type=int, default=12345)
    p.add_argument("--exclude-species", type=json.loads, default=[], help="JSON list of species IDs to prune")
    args = vars(parser.parse_args())
    globals()[args.pop("action")](**args)
