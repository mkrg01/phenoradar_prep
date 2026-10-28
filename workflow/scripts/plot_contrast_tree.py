#!/usr/bin/env python3
"""Draw compressed contrast trees without rerunning selection or inference."""
import argparse
import json
import math
from pathlib import Path

from common import read_tsv


def tip_identity(row, samples):
    """Resolve explicit names and runs without splitting computational sample IDs."""
    records = samples.get(row["species"], [row])
    names = {r.get("scientific_name") or r.get("species_id", "").replace("_", " ")
             for r in records}
    names.discard("")
    if len(names) > 1:
        raise ValueError(f"conflicting scientific names for plot tip: {row['species']}")
    name = next(iter(names), row.get("scientific_name") or row.get("species_id", "").replace("_", " "))
    runs = sorted({r["run"] for r in records if r.get("run")})
    # Legacy metadata may contain only an opaque sample ID. Keep it upright.
    return name or row["species"], bool(name), ", ".join(runs)


def paired_clades(phy, rows):
    """Validate that each recorded pair is exactly one two-tip summary clade."""
    pairs = {}
    for tip in phy.leaves():
        pair = rows[tip.name].get("contrast_pair_id", "")
        if pair:
            pairs.setdefault(pair, []).append(tip)
    clades = []
    for pair, tips in pairs.items():
        if len(tips) != 2:
            raise ValueError(f"contrast pair {pair} must have exactly two representative tips")
        root = phy.common_ancestor(tips)
        if set(root.leaf_names()) != {tip.name for tip in tips}:
            raise ValueError(f"contrast pair {pair} is not a summary-tree clade")
        clades.append(root)
    return clades


def plot(tree, metadata, summary, outdir, samples=None):
    """Publication and ID phylograms, sharing topology, counts and pair colors."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.font_manager import FontProperties
    from matplotlib.legend import Legend
    from matplotlib.lines import Line2D
    from matplotlib.textpath import TextToPath
    from ete4 import Tree
    info = json.loads(Path(summary).read_text())
    metadata_rows = read_tsv(metadata)
    rows = {r["species"]: r for r in metadata_rows if r["is_representative"] == "1"}
    if not Path(tree).read_text().strip():
        if rows or info["summary_species"]:
            raise ValueError("empty summary tree has nonempty representatives")
        fig, ax = plt.subplots(figsize=(7.2, 2))
        ax.set_axis_off()
        ax.text(0.5, 0.5, f"No species with observed {info['trait']}\n0 contrast pairs",
                ha="center", va="center", transform=ax.transAxes)
        Path(outdir).mkdir(parents=True, exist_ok=True)
        for stem in ["summary_tree", "summary_tree_id"]:
            for extension in ["pdf", "svg"]:
                fig.savefig(Path(outdir) / f"{stem}.{extension}")
        plt.close(fig)
        return
    phy = Tree(Path(tree).read_text(), parser=0)
    # Ladderization only orders siblings for drawing; topology and lengths stay intact.
    counts = {}
    for node in phy.traverse("postorder"):
        counts[node] = sum(counts[c] for c in node.children) if node.children else 1
        node.children.sort(key=lambda c: counts[c], reverse=True)
    tips = list(phy.leaves())
    if set(phy.leaf_names()) != set(rows):
        raise ValueError("plot labels differ from summary tree")
    clades = paired_clades(phy, rows)
    paired_nodes = {node for root in clades for node in root.traverse()}
    paired_edges = {node for root in clades for node in root.descendants()}
    sample_rows = {}
    for row in read_tsv(samples) if samples else []:
        sample_rows.setdefault(row["species"], []).append(row)
    identities = {node: tip_identity(rows[node.name], sample_rows) for node in tips}
    states = sorted({r[info["trait"]] for r in rows.values()})
    state_order = info.get("state_order", sorted({r[info["trait"]] for r in metadata_rows if r[info["trait"]] != ""}))
    # Preserve binary colors even in legacy summaries containing only one state.
    if set(state_order) <= {"0", "1"}:
        state_order = ["0", "1"]
    if "" in states or len(state_order) > 2 or not set(states) <= set(state_order):
        raise ValueError("plot requires at most two observed trait states matching the recorded state order")
    colors = dict(zip(state_order, ["#595959", "#009E73"]))
    labels = {state: info.get("trait_labels", {}).get(state, state) for state in states}
    distance = {phy: 0.0}
    for node in phy.traverse("preorder"):
        for child in node.children:
            distance[child] = distance[node] + (child.dist or 0)

    # Lay out in physical points so 8 pt labels stay 8 pt in the exported figure.
    italic = FontProperties(family="DejaVu Sans", style="italic", size=8)
    upright = FontProperties(family="DejaVu Sans", size=8)
    text_path = TextToPath()
    def text_width(value, font):
        return text_path.get_text_width_height_descent(str(value), font, ismath=False)[0]
    def suffix(node, detailed):
        row = rows[node.name]
        fields = [f"n = {row['n_species_in_group']}"]
        if detailed:
            run = identities[node][2]
            fields += [f"run = {run}" if run else f"sample = {node.name}",
                       f"group = {row.get('group') or '-'}", f"pair = {row.get('contrast_pair_id') or '-'}"]
        return " (" + "; ".join(fields) + ")"
    def label_width(node, detailed):
        name, scientific, _ = identities[node]
        return text_width(name, italic if scientific else upright) + text_width(suffix(node, detailed), upright)
    margin, gap, row_height, legend_height = 12, 10, 13, 58
    # Keep the same branch-length scale and row spacing in both figures.
    tree_width = max(180, 7.2 * 72 - margin * 2 - gap - max(label_width(n, False) for n in tips))
    name_x = margin + tree_width + gap
    legend_width = max(text_width(info["trait"], upright),
                       sum(text_width(labels[state], upright) + 20 for state in states) + gap,
                       text_width("Contrast pair clade", upright) + 24 if clades else 0)
    height = (len(tips) - 1) * row_height + 44 + margin + legend_height
    maximum = max(distance.values())
    x = {node: margin + value / (maximum or 1) * tree_width for node, value in distance.items()}
    y = {node: 44 + i * row_height for i, node in enumerate(tips)}
    for node in phy.traverse("postorder"):
        if not node.is_leaf:
            y[node] = sum(y[c] for c in node.children) / len(node.children)

    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 8,
                         "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
                         "lines.solid_capstyle": "butt", "savefig.facecolor": "white",
                         "text.parse_math": False}):
        for detailed, stem in [(False, "summary_tree"), (True, "summary_tree_id")]:
            width = max(7.2 * 72, name_x + max(label_width(n, detailed) for n in tips) + margin,
                        margin * 2 + legend_width)
            fig = plt.figure(figsize=(width / 72, height / 72))
            ax = fig.add_axes([0, 0, 1, 1], xlim=(0, width), ylim=(0, height))
            ax.set_axis_off()
            for node in phy.traverse():
                if node.children:
                    ax.plot([x[node]] * 2, [min(y[c] for c in node.children), max(y[c] for c in node.children)],
                            color="#C62828" if node in paired_nodes else "#202020", linewidth=0.7)
                for child in node.children:
                    ax.plot([x[node], x[child]], [y[child]] * 2,
                            color="#C62828" if child in paired_edges else "#202020", linewidth=0.7)
            for node in tips:
                row = rows[node.name]
                ax.plot([x[node] + 3, name_x - 4], [y[node]] * 2, color="#BBBBBB",
                        linewidth=0.4, linestyle=(0, (1.5, 2)), zorder=1)
                ax.scatter(x[node], y[node], color=colors[row[info["trait"]]], s=14,
                           edgecolors="white", linewidths=0.3, zorder=3)
                name, scientific, _ = identities[node]
                font = italic if scientific else upright
                ax.text(name_x, y[node], name, fontproperties=font, va="center", color="#202020")
                ax.text(name_x + text_width(name, font), y[node], suffix(node, detailed),
                        fontproperties=upright, va="center", color="#202020")
            ax.legend(handles=[Line2D([], [], color=colors[s], marker="o", markersize=4,
                                    linestyle="", label=labels[s]) for s in states],
                      title=info["trait"], title_fontsize=8,
                      loc="upper left", bbox_to_anchor=(margin, height - margin), bbox_transform=ax.transData,
                      ncol=len(states), frameon=False, borderaxespad=0,
                      handlelength=0.8, handletextpad=0.5, columnspacing=1.5)
            if clades:
                ax.add_artist(Legend(ax, [Line2D([], [], color="#C62828", linewidth=0.7)],
                                     ["Contrast pair clade"], loc="upper left",
                                     bbox_to_anchor=(margin, height - margin - 32),
                                     bbox_transform=ax.transData, frameon=False, borderaxespad=0,
                                     handlelength=1.2, handletextpad=0.6))
            if maximum > 0:
                target = maximum / 5
                magnitude = 10 ** math.floor(math.log10(target))
                scale = max(value * magnitude for value in [1, 2, 5] if value * magnitude <= target)
                end = margin + scale / maximum * tree_width
                ax.plot([margin, end], [22, 22], color="#202020", linewidth=0.7)
                for position in [margin, end]:
                    ax.plot([position] * 2, [20, 24], color="#202020", linewidth=0.7)
                ax.text(margin, 9, f"{scale:g} substitutions per site", va="center", fontsize=7)
            Path(outdir).mkdir(parents=True, exist_ok=True)
            for extension in ["pdf", "svg"]:
                fig.savefig(Path(outdir) / f"{stem}.{extension}")
            plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["tree", "metadata", "summary", "outdir"]:
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--samples", help="Sample manifest with explicit scientific_name and run columns")
    plot(**vars(parser.parse_args()))
