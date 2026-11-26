#!/usr/bin/env python3
"""Trait-guided nwkit skims, representative manifests, pair tables and figures."""
import argparse
import json
import math
from collections import Counter
from pathlib import Path

from phylogeny_outgroup import outgroup_ids, prune_in_tip_order, read_outgroup, validate_root
from types import SimpleNamespace

from common import RANDOM_SEED, atomic_writer, file_record, now, read_tsv, write_json, write_tsv
from phylogeny_root import nwkit_backend, species_rows
from species_traits import read_species_traits
from sample_identity import species_id, traits_for_samples
from expanded_pair_clades import COLUMN, TABLE_FIELDS, annotate
from contrast_metadata import CLADE_COLUMN


def group_traits(tree, rows):
    """Assign maximal uniform trait clades without selecting representatives."""
    nwkit_backend()
    import pandas as pd
    from nwkit.skim import mark_traits_to_nodes, add_group_ids
    args = SimpleNamespace(group_by="trait")
    names = set(tree.leaf_names())
    frame = pd.DataFrame(sorted(rows, key=lambda r: r["leaf_name"]))
    if set(frame["leaf_name"]) != names or frame["leaf_name"].duplicated().any():
        raise ValueError("skim trait rows must match tree tips exactly")
    # No missing traits enter grouping.
    if frame["trait"].eq("").any():
        raise ValueError("skim requires observed traits")
    tree = mark_traits_to_nodes(tree, frame, args)
    return add_group_ids(frame, tree)


def save_skim(tree, frame, prefix, seed, contrastive=False, topology_only=False):
    """Write representative views after their clade memberships are known."""
    from nwkit.skim import add_contrastive_clade_ids, sample_from_groups
    args = SimpleNamespace(group_by="trait", filter_by="busco_percent", filter_mode="descending",
                           retain_per_clade=1, prioritize_non_missing=True, seed=seed)
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
            prune_in_tip_order(tree, sampled["leaf_name"].tolist())
            handle.write(tree.write(parser=9 if topology_only else 0) + "\n")
        # There is no Newick representation of zero tips. An empty file means
        # no contrastive clades; summary.json and the TSV headers disambiguate it.
    return frame.to_dict("records"), sampled.to_dict("records")


def skim(tree, rows, prefix, seed, contrastive=False, topology_only=False):
    """Select NCBI inference representatives or write molecular-tree views."""
    return save_skim(tree, group_traits(tree, rows), prefix, seed, contrastive, topology_only)


def trait_clade_pairs(tree, rows):
    """Find opposite-state sister clades directly on the input molecular tree."""
    frame = group_traits(tree, rows)
    from nwkit.skim import add_contrastive_clade_ids

    candidates = add_contrastive_clade_ids(frame.copy(), tree)
    pairs, unresolved = [], []
    # In a right-facing plot, the last Newick tip appears at the top.
    tip_order = {name: index for index, name in enumerate(reversed(list(tree.leaf_names())))}
    clades = [part for _, part in candidates.dropna(subset=["contrastive_clade"]).groupby("contrastive_clade")]
    for candidate in sorted(clades, key=lambda part: min(tip_order[name] for name in part["leaf_name"])):
        groups = candidate.drop_duplicates("group").sort_values("trait")
        if len(groups) != 2 or groups["trait"].nunique() != 2:
            unresolved.append(sorted(candidate["leaf_name"]))
            continue
        left, right = groups.to_dict("records")
        pairs.append({"contrast_pair_id": len(pairs) + 1,
                      "state_a": left["trait"], "state_b": right["trait"],
                      "group_a": left["group"], "group_b": right["group"]})
    return frame, pairs, unresolved


def prepare(samples, metadata, traits, tree, outdir, trait="C4", seed=RANDOM_SEED):
    from ete4 import Tree
    rows = species_rows(samples)
    annotation = traits_for_samples(rows.values(), read_species_traits(traits, trait))
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
    if len({species_id(rows[n]) for n in representatives}) < 4:
        raise ValueError("fewer than four representative species; BUSCO inference requires four")
    write_tsv(out / "samples.tsv", list(next(iter(rows.values()))), [rows[n] for n in sorted(representatives)])
    write_tsv(out / "traits.tsv", ["species", "species_id", "trait", "busco_percent", "role"],
              [{"species": n, "species_id": species_id(rows[n]), "trait": annotation.get(n, ""), "busco_percent": scores[n],
                "role": "observed" if n in eligible else "missing_trait"}
               for n in rows])
    write_json(out / "selection.json", {"created_at": now(), "species_trait": file_record(traits),
               "trait": trait, "states": states,
               "dataset_species": len({species_id(r) for r in rows.values()}), "dataset_samples": len(rows),
               "observed_species": len({species_id(rows[n]) for n in eligible}), "observed_samples": len(eligible),
               "inference_species": len({species_id(rows[n]) for n in representatives}),
               "inference_samples": len(representatives),
               "missing_trait_rows": sorted(set(rows) - set(annotation)), "seed": seed,
               "nwkit": nwkit_backend().__version__, "ncbi_tree": file_record(tree)})


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
    validate_root(tree, outgroup)
    if any(n.dist is None or not math.isfinite(n.dist) or n.dist < 0
           for n in tree.traverse() if not n.is_root):
        raise ValueError("inferred tree requires finite nonnegative branch lengths")
    return tree


def summarize(tree, selection_dir, outgroup_file, outdir, seed=RANDOM_SEED):
    selection, out = Path(selection_dir), Path(outdir)
    plan = json.loads((selection / "selection.json").read_text())
    all_traits = read_tsv(selection / "traits.tsv")
    first_all = read_tsv(selection / "ncbi_skim.all.tsv")
    first_reps = read_tsv(selection / "ncbi_skim.sampled.tsv")
    names = {r["leaf_name"] for r in first_reps}
    outgroup = read_outgroup(outgroup_file)
    inferred = rooted_tree(tree, names, outgroup)
    rows = [{"leaf_name": r["leaf_name"], "trait": r["trait"], "busco_percent": float(r["busco_percent"])}
            for r in first_reps]
    first_representatives = {r["group"]: r["leaf_name"] for r in first_reps}
    tip_for_species = {r["leaf_name"]: first_representatives[r["group"]] for r in first_all}
    return summarize_tree(inferred, rows, all_traits, tip_for_species, out, plan["trait"], seed, outgroup,
                          {"outgroup_file": file_record(outgroup_file), "inferred_tree": file_record(tree),
                           "selection": file_record(selection / "selection.json"),
                           "trait_labels": plan.get("trait_labels", {}),
                           "assignment": "non-representative species inherit their NCBI skim group membership"})


def summarize_tree(inferred, rows, all_traits, tip_for_species, outdir, trait, seed, outgroup, provenance):
    """Assign molecular clades and blocks before selecting display samples."""
    from ete4 import Tree
    outgroup = outgroup_ids(outgroup)
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    if rows:
        frame, pairs, unresolved = trait_clade_pairs(inferred, rows)
        second_all = frame.to_dict("records")
    else:
        second_all, pairs, unresolved = [], [], []
    # Map original species through their input-tree tips to molecular groups.
    second_group = {r["leaf_name"]: r["group"] for r in second_all}
    members = {name: second_group[tip] for name, tip in tip_for_species.items()}
    sample_counts = Counter(members.values())
    taxa = {r["species"]: r.get("species_id", r["species"]) for r in all_traits}
    counts = {g: len({taxa[n] for n, group in members.items() if group == g}) for g in sample_counts}
    paired_groups = {}
    for pair in pairs:
        for group in [pair["group_a"], pair["group_b"]]:
            if group in paired_groups:
                raise ValueError("one clade assigned to more than one contrast pair")
            paired_groups[group] = pair["contrast_pair_id"]
        for side in ["a", "b"]:
            group = pair[f"group_{side}"]
            pair[f"n_species_{side}"] = counts[group]
            pair[f"n_samples_{side}"] = sample_counts[group]
    metadata = []
    for row in all_traits:
        name = row["species"]
        group = members.get(name, "")
        metadata.append({"species": name, "species_id": taxa[name], trait: row["trait"],
                         "role": "outgroup" if name in outgroup else row["role"],
                         CLADE_COLUMN: group,
                         "n_species_in_trait_clade": counts.get(group, ""),
                         "n_samples_in_trait_clade": sample_counts.get(group, ""),
                         "contrast_pair_id": paired_groups.get(group, "")})
    metadata, expanded_clades, expanded_report = annotate(
        inferred if rows else None, metadata, trait, tip_for_species=tip_for_species)

    # Representatives affect only the exported views and their tip labels.
    # Pair IDs and block membership have already been fixed on the input tree.
    if rows:
        _, second_reps = save_skim(inferred.copy(), frame, out / "summary_tree", seed)
        summary_tree = Tree((out / "summary_tree.nwk").read_text(), parser=0)
        skim(summary_tree, [{k: r[k] for k in ["leaf_name", "trait", "busco_percent"]}
                            for r in second_reps], out / "contrastive", seed, contrastive=True)
    else:
        second_reps = []
        for prefix in ["summary_tree", "contrastive"]:
            fields = ["leaf_name", "trait", "busco_percent", "group"]
            if prefix == "contrastive":
                fields += ["contrastive_clade"]
            for suffix in ["all", "sampled"]:
                write_tsv(out / f"{prefix}.{suffix}.tsv", fields, [])
            with atomic_writer(out / f"{prefix}.nwk") as handle:
                handle.write("")
    final_reps = {r["group"]: r["leaf_name"] for r in second_reps}
    for row in metadata:
        row["representative"] = final_reps.get(row[CLADE_COLUMN], "")
        row["is_representative"] = int(row["representative"] == row["species"])
    for pair in pairs:
        for side in ["a", "b"]:
            pair[f"representative_{side}"] = final_reps[pair[f"group_{side}"]]
    fields = ["contrast_pair_id", "state_a", "state_b", "representative_a", "representative_b",
              "group_a", "group_b", "n_species_a", "n_species_b", "n_samples_a", "n_samples_b"]
    write_tsv(out / "contrast_pairs.tsv", fields, pairs)
    write_tsv(out / "species_metadata.tsv", ["species", "species_id", trait, "role", CLADE_COLUMN, "representative",
              "is_representative", "n_species_in_trait_clade", "n_samples_in_trait_clade", "contrast_pair_id", COLUMN], metadata)
    write_tsv(out / "phylogenetic_blocks.tsv", TABLE_FIELDS, expanded_clades)
    report = {"created_at": now(), "trait": trait, "pair_ids_follow_tree": True,
               "state_order": sorted({r["trait"] for r in all_traits if r["trait"] != ""}),
               "nwkit": nwkit_backend().__version__, "outgroup": outgroup,
               "summary_species": len({taxa[r["leaf_name"]] for r in second_reps}),
               "summary_samples": len(second_reps), "contrast_pairs": len(pairs), "seed": seed,
               "pair_selection": "opposite-state sister clades on the input molecular tree",
               "block_selection": "maximal clades on the input molecular tree",
               "summary_representatives": "display only",
               "unresolved_contrastive_clades": unresolved, **provenance, **expanded_report}
    write_json(out / "summary.json", report)
    return report


def from_tree(tree, tree_qc, samples, metadata, traits, outdir, trait="C4", seed=RANDOM_SEED, exclude_species=(),
              trait_labels=None):
    """Assign pairs directly on a completed species tree, optionally pruning species."""
    if not isinstance(exclude_species, (list, tuple)) or any(not isinstance(n, str) for n in exclude_species):
        raise ValueError("exclude_species must be a list of species IDs")
    if len(set(exclude_species)) != len(exclude_species):
        raise ValueError("duplicate excluded species IDs")
    if not isinstance(trait, str) or not trait.strip() or trait in {
            "species", "species_id", "role", "group", CLADE_COLUMN, "representative", "is_representative",
            "n_species_in_group", "n_samples_in_group", "n_species_in_trait_clade",
            "n_samples_in_trait_clade", "contrast_pair_id", COLUMN, "expanded_pair_clade_id"}:
        raise ValueError("trait must name a non-reserved phenotype column")
    if type(seed) is not int or seed <= 0:
        raise ValueError("seed must be a positive integer")
    sample_rows = read_tsv(samples)
    names = {r["species"] for r in sample_rows}
    qc = json.loads(Path(tree_qc).read_text())
    if qc.get("species") != len(names):
        raise ValueError("species-tree QC differs from species manifest")
    outgroup = outgroup_ids(qc.get("outgroup"))
    inferred = rooted_tree(tree, names, outgroup)
    annotation = traits_for_samples(sample_rows, read_species_traits(traits, trait))
    retained = names - set(exclude_species)
    eligible = {n for n in retained if annotation.get(n, "") != ""}
    states = sorted({annotation[n] for n in eligible})
    if len(states) > 2:
        raise ValueError(f"contrast analysis supports at most two observed trait states; found {states}")
    source_states = sorted({annotation[n] for n in names if annotation.get(n, "") != ""})
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
            prune_in_tip_order(inferred, sorted(eligible))
            handle.write(inferred.write(parser=0) + "\n")
    rows = [{"leaf_name": n, "trait": annotation[n], "busco_percent": scores[n]} for n in sorted(eligible)]
    taxa = {r["species"]: species_id(r) for r in sample_rows}
    all_traits = [{"species": n, "species_id": taxa[n], "trait": annotation.get(n, ""),
                   "role": "observed" if n in eligible else "missing_trait"} for n in sorted(retained)]
    return summarize_tree(inferred, rows, all_traits, {n: n for n in eligible}, out, trait, seed, outgroup,
                          {"mode": "inferred_tree", "inferred_tree": file_record(tree),
                           "tree_qc": file_record(tree_qc), "samples": file_record(samples),
                           "metadata": file_record(metadata), "species_trait": file_record(traits),
                           "source_species": len(set(taxa.values())), "source_samples": len(names),
                           "retained_species": len({taxa[n] for n in retained}), "retained_samples": len(retained),
                           "observed_species": len({taxa[n] for n in eligible}), "observed_samples": len(eligible),
                           "states": states,
                           "trait_labels": trait_labels or {},
                           "state_order": source_states if len(source_states) <= 2 else states,
                           "excluded_species": sorted(names - retained),
                           "missing_trait_species": sorted(retained - eligible),
                           "rooting": "source_root_inherited", "outgroup_in_observed_tree": bool(set(outgroup) & eligible),
                           "retained_outgroup": sorted(set(outgroup) & eligible),
                           "reestimated": False,
                           "assignment": "membership on the observed-species subtree of the input molecular tree"})



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    p = sub.add_parser("prepare")
    for name in ["samples", "metadata", "traits", "tree", "outdir"]:
        p.add_argument("--" + name, required=True)
    p.add_argument("--trait", default="C4")
    p.add_argument("--seed", type=int, default=RANDOM_SEED)
    p = sub.add_parser("summarize")
    for name in ["tree", "selection-dir", "outgroup-file", "outdir"]:
        p.add_argument("--" + name, required=True)
    p.add_argument("--seed", type=int, default=RANDOM_SEED)
    p = sub.add_parser("from_tree")
    for name in ["tree", "tree-qc", "samples", "metadata", "traits", "outdir"]:
        p.add_argument("--" + name, required=True)
    p.add_argument("--trait", default="C4")
    p.add_argument("--seed", type=int, default=RANDOM_SEED)
    p.add_argument("--exclude-species", type=json.loads, default=[], help="JSON list of species IDs to prune")
    args = vars(parser.parse_args())
    globals()[args.pop("action")](**args)
