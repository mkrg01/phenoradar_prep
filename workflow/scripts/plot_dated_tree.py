#!/usr/bin/env python3
"""Plot a fitted time tree with species traits and frozen NCBI taxonomy."""
import argparse
from collections import Counter
from itertools import groupby
import logging
import math
from pathlib import Path

from common import atomic_writer, read_tsv, write_tsv
from infer_phylogeny import read_tree
from sample_identity import species_id
from species_traits import read_species_traits
from taxonomy_check import Taxonomy, load_samples


TRAIT_COLORS = {"0": "#888888", "1": "#006C4D", "": "#EAEAEA"}
TRAIT_LABELS = {"0": "0", "1": "1", "": "Missing"}
TIP_FIELDS = ["plot_row", "species", "scientific_name", "taxid", "order", "family",
              "trait", "trait_state", "trait_label"]


def tip_annotations(tips, samples, taxonomy, traits, trait):
    _, manifest = load_samples(samples)
    values = read_species_traits(traits, trait)
    db = Taxonomy(taxonomy)
    rows = []
    try:
        for i, tip in enumerate(tips, 1):
            sample = manifest.get(tip.name)
            if sample is None:
                raise ValueError(f"dated tree tip is absent from the sample manifest: {tip.name}")
            state = values.get(species_id(sample), "")
            if state not in TRAIT_COLORS:
                raise ValueError(f"expected a normalized trait state for {tip.name}: {state!r}")
            lineage = db.lineage(int(sample["taxid"]))
            rows.append(dict(plot_row=i, species=tip.name, scientific_name=sample["scientific_name"],
                             taxid=sample["taxid"], order=lineage.get("order", (None, ""))[1],
                             family=lineage.get("family", (None, ""))[1], trait=trait,
                             trait_state=state, trait_label=TRAIT_LABELS[state]))
    finally:
        db.db.close()
    return rows


def tree_segments(tree, tips, ages):
    age_rows = read_tsv(ages)
    node_ages = {row["node"]: float(row["age_ma"]) for row in age_rows}
    names = [node.name for node in tree.traverse() if not node.is_leaf]
    if (len(node_ages) != len(age_rows) or len(set(names)) != len(names)
            or set(node_ages) != set(names)
            or any(not math.isfinite(age) or age < 0 for age in node_ages.values())):
        raise ValueError("node ages must match the dated Newick's named internal nodes")
    root_age = node_ages[tree.name]
    if root_age <= 0:
        raise ValueError("root age must be positive")
    x = {node: 0.0 if node.is_leaf else node_ages[node.name] for node in tree.traverse()}
    y = {node: i for i, node in enumerate(tips, 1)}
    internal, terminal = [], []
    for node in tree.traverse("postorder"):
        if node.children:
            low = min(y[child] for child in node.children)
            high = max(y[child] for child in node.children)
            y[node] = (low + high) / 2
            internal.append([(x[node], low), (x[node], high)])
            for child in node.children:
                edges = terminal if child.is_leaf else internal
                edges.append([(x[node], y[child]), (x[child], y[child])])
    return internal, terminal, root_age


def track_runs(values):
    start = 0.5
    for value, items in groupby(values):
        end = start + sum(1 for _ in items)
        yield value, start, end
        start = end


def family_labels(rows, runs, spacing):
    candidates = []
    for value, lo, hi in runs:
        # Run boundaries lie halfway between tips and truncate to slice indices.
        if any(row["trait_state"] == "1" for row in rows[int(lo):int(hi)]):
            candidates.append((value, lo, hi))

    capacity = math.floor(len(rows) / spacing)
    if len(candidates) > capacity:
        logging.getLogger(__name__).warning(
            "%d family runs with trait state 1 exceed the figure's label capacity. See the tips TSV.",
            len(candidates) - capacity)
        candidates.sort(key=lambda run: run[2] - run[1], reverse=True)
        candidates = sorted(candidates[:capacity], key=lambda run: run[1])

    # Preserve the order of the bars and reserve space for the remaining names.
    lower = 0.5 + spacing / 2
    upper = len(rows) + 0.5 - spacing / 2
    for i, (value, lo, hi) in enumerate(candidates):
        midpoint = (lo + hi) / 2
        remaining = len(candidates) - i - 1
        position = min(max(midpoint, lower), upper - remaining * spacing)
        yield value or "Unassigned", midpoint, position
        lower = position + spacing


def draw(internal, terminal, root_age, rows, trait, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection, PolyCollection
    from matplotlib.patches import Patch
    from matplotlib.ticker import MaxNLocator

    count = len(rows)
    show_names = count <= 60
    width = 7
    height = min(9, max(2, 0.85 + 0.135 * count))
    bottom, top = 0.46, 0.25  # Margins in inches.
    tracks = [5, 0.15, 1.1] + ([2] if show_names else [])
    trait_title = "C₄" if trait == "C4" else trait
    with plt.rc_context({"pdf.fonttype": 42, "svg.fonttype": "none", "text.parse_math": False,
                         "font.family": "DejaVu Sans", "font.size": 7.5, "text.color": "#202020"}):
        fig, (tree_ax, trait_ax, family_ax, *name_axes) = plt.subplots(
            1, len(tracks), figsize=(width, height), width_ratios=tracks, sharey=True)
        try:
            fig.subplots_adjust(left=0.032, right=0.976, top=1 - top / height,
                                bottom=bottom / height, wspace=0.04)
            tree_ax.set_ylim(count + 0.5, 0.5)
            tree_ax.set_xlim(root_age * 1.02, 0)
            tree_ax.xaxis.set_major_locator(MaxNLocator(nbins=6))
            tree_ax.set_axisbelow(True)
            tree_ax.grid(axis="x", color="#E6E6E6", linewidth=0.4)
            # Softer terminal edges keep the deeper relationships visible at print size.
            dense = count > 300
            tree_ax.add_collection(LineCollection(terminal, colors="#929292" if dense else "#404040",
                                                   linewidths=0.18 if dense else 0.55))
            tree_ax.add_collection(LineCollection(internal, colors="#404040",
                                                   linewidths=0.28 if dense else 0.55))
            tree_ax.set_xlabel("Age (Ma)", labelpad=3)
            tree_ax.tick_params(axis="x", length=2.5, width=0.5, pad=2.5)
            tree_ax.tick_params(axis="y", left=False, labelleft=False)
            for edge in ["left", "right", "top"]:
                tree_ax.spines[edge].set_visible(False)
            tree_ax.spines["bottom"].set(color="#808080", linewidth=0.5)
            for ax in [trait_ax, family_ax, *name_axes]:
                ax.set_axis_off()
                ax.set_xlim(0, 1)
            trait_ax.set_title(trait_title, loc="right", fontsize=8, weight="bold", pad=4)
            for ax, title in zip([family_ax, *name_axes], ["Family", "Species"]):
                ax.set_title(title, loc="left", fontsize=8, weight="bold", pad=4)

            trait_runs = list(track_runs(row["trait_state"] for row in rows))
            trait_ax.add_collection(PolyCollection(
                [[(0, lo), (1, lo), (1, hi), (0, hi)] for _, lo, hi in trait_runs],
                facecolors=[TRAIT_COLORS[state] for state, _, _ in trait_runs], edgecolors="none"))
            family_runs = list(track_runs(row["family"] for row in rows))
            # A separate bar marks each contiguous run of tips in the same family.
            rows_per_inch = count / (height - bottom - top)
            gap = min(0.24, 0.013 * rows_per_inch)
            family_ax.add_collection(PolyCollection(
                [[(0, lo + gap), (0.07, lo + gap), (0.07, hi - gap), (0, hi - gap)]
                 for _, lo, hi in family_runs],
                facecolors=["#A5A5A5" if value else "#EAEAEA" for value, _, _ in family_runs],
                edgecolors="none"))
            for label, midpoint, position in family_labels(rows, family_runs, 0.12 * rows_per_inch):
                family_ax.plot([0.07, 0.12, 0.18], [midpoint, midpoint, position],
                               color="#777777", linewidth=0.4)
                family_ax.text(0.21, position, label, va="center")
            if show_names:
                for row in rows:
                    name_axes[0].text(0, row["plot_row"], row["scientific_name"],
                                      style="italic", va="center")

            counts = Counter(row["trait_state"] for row in rows)
            handles = [Patch(facecolor=color, edgecolor="#A5A5A5", linewidth=0.35, label=TRAIT_LABELS[state])
                       for state, color in TRAIT_COLORS.items() if counts[state]]
            fig.legend(handles=handles, ncol=len(handles), frameon=False,
                       loc="lower right", bbox_to_anchor=(0.976, 0.12 / height), borderaxespad=0,
                       handlelength=0.9, handletextpad=0.45, columnspacing=1)
            for extension in ["pdf", "svg", "png"]:
                with atomic_writer(outdir / f"species_tree.dated.{extension}", "wb") as handle:
                    fig.savefig(handle, format=extension, dpi=600, facecolor="white")
        finally:
            plt.close(fig)


def plot(tree, ages, samples, taxonomy, traits, trait, outdir):
    phy = read_tree(tree)
    tips = list(phy.leaves())
    rows = tip_annotations(tips, samples, taxonomy, traits, trait)
    internal, terminal, root_age = tree_segments(phy, tips, ages)
    outdir = Path(outdir)
    draw(internal, terminal, root_age, rows, trait, outdir)
    write_tsv(outdir / "species_tree.dated.tips.tsv", TIP_FIELDS, rows)
    print(f"Plotted {len(rows):,} tips with {trait} and NCBI taxonomy")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["tree", "ages", "samples", "taxonomy", "traits", "trait", "outdir"]:
        parser.add_argument("--" + name, required=True)
    plot(**vars(parser.parse_args()))
