#!/usr/bin/env python3
"""Trait-guided nwkit skims, representative manifests, pair tables and figures."""
import argparse
import json
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


def summarize(tree, selection_dir, outgroup_file, outdir, seed=12345):
    from ete4 import Tree
    selection, out = Path(selection_dir), Path(outdir)
    plan = json.loads((selection / "selection.json").read_text())
    all_traits = read_tsv(selection / "traits.tsv")
    first_all = read_tsv(selection / "ncbi_skim.all.tsv")
    first_reps = read_tsv(selection / "ncbi_skim.sampled.tsv")
    names = {r["leaf_name"] for r in first_reps}
    outgroup = Path(outgroup_file).read_text().strip()
    if outgroup not in names:
        raise ValueError("outgroup is absent from the selected representatives")
    inferred = Tree(Path(tree).read_text(), parser=0)
    if set(inferred.leaf_names()) != names:
        raise ValueError("inferred tree tips differ from representative manifest")
    if len(inferred.children) != 2 or not any(c.is_leaf and c.name == outgroup for c in inferred.children):
        raise ValueError("inferred tree is not rooted with the selected outgroup")
    rows = [{"leaf_name": r["leaf_name"], "trait": r["trait"], "busco_percent": float(r["busco_percent"])}
            for r in first_reps]
    second_all, second_reps = skim(inferred, rows, out / "summary_tree", seed)
    summary_tree = Tree((out / "summary_tree.nwk").read_text(), parser=0)
    contrast_all, contrast_reps = skim(summary_tree, [
        {k: r[k] for k in ["leaf_name", "trait", "busco_percent"]} for r in second_reps
    ], out / "contrastive", seed, contrastive=True)
    # Compose two many-to-one maps; do not rely on IDs from different skims
    # having the same meaning or on an unchecked dict overwriting duplicates.
    second_group = {r["leaf_name"]: r["group"] for r in second_all}
    first_to_second = {r["group"]: second_group[r["leaf_name"]] for r in first_reps}
    final_reps = {r["group"]: r["leaf_name"] for r in second_reps}
    members = {r["leaf_name"]: first_to_second[r["group"]] for r in first_all}
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
        metadata.append({"species": name, plan["trait"]: row["trait"],
                         "role": "outgroup" if name == outgroup else row["role"],
                         "group": group, "representative": final_reps.get(group, ""),
                         "is_representative": int(final_reps.get(group) == name),
                         "n_species_in_group": counts.get(group, ""),
                         "contrast_pair_id": paired_groups.get(group, "")})
    fields = ["contrast_pair_id", "state_a", "state_b", "representative_a", "representative_b",
              "group_a", "group_b", "n_species_a", "n_species_b"]
    write_tsv(out / "contrast_pairs.tsv", fields, pairs)
    write_tsv(out / "species_metadata.tsv", list(metadata[0]), metadata)
    write_json(out / "summary.json", {"created_at": now(), "trait": plan["trait"], "nwkit": "0.27.0",
               "outgroup": outgroup, "outgroup_file": file_record(outgroup_file),
               "inferred_tree": file_record(tree), "selection": file_record(selection / "selection.json"),
               "summary_species": len(second_reps), "contrast_pairs": len(pairs), "seed": seed,
               "unresolved_contrastive_clades": unresolved,
               "assignment": "non-representative species inherit their NCBI skim group membership"})



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
    args = vars(parser.parse_args())
    globals()[args.pop("action")](**args)
