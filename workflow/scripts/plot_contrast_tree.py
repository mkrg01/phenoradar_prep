#!/usr/bin/env python3
"""Draw compressed contrast trees without rerunning selection or inference."""
import argparse
import json
import math
from pathlib import Path

from common import read_tsv


def plot(tree, metadata, summary, outdir):
    """Vector phylogram with italic names, aligned annotations and a scale bar."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.font_manager import FontProperties
    from matplotlib.lines import Line2D
    from matplotlib.textpath import TextPath
    from ete4 import Tree
    info = json.loads(Path(summary).read_text())
    rows = {r["species"]: r for r in read_tsv(metadata) if r["is_representative"] == "1"}
    phy = Tree(Path(tree).read_text(), parser=0)
    # Ladderization only orders siblings for drawing; topology and lengths stay intact.
    counts = {}
    for node in phy.traverse("postorder"):
        counts[node] = sum(counts[c] for c in node.children) if node.children else 1
        node.children.sort(key=lambda c: counts[c], reverse=True)
    tips = list(phy.leaves())
    if set(phy.leaf_names()) != set(rows):
        raise ValueError("plot labels differ from summary tree")
    states = sorted({r[info["trait"]] for r in rows.values()})
    colors = dict(zip(states, ["#595959", "#009E73"]))
    distance = {phy: 0.0}
    for node in phy.traverse("preorder"):
        for child in node.children:
            distance[child] = distance[node] + (child.dist or 0)

    # Lay out in physical points so 8 pt labels stay 8 pt in the exported figure.
    italic = FontProperties(family="DejaVu Sans", style="italic", size=8)
    upright = FontProperties(family="DejaVu Sans", size=8)
    def text_width(value, font):
        return TextPath((0, 0), str(value), prop=font).get_extents().width
    names = {node: node.name.replace("_", " ") for node in tips}
    name_width = max(text_width(name, italic) for name in names.values())
    columns = [("n", "n_species_in_group"), ("Group", "group"), ("Pair", "contrast_pair_id")]
    widths = [max(text_width(title, upright), *(text_width(r[field], upright)
                  for r in rows.values() if r[field] != "")) for title, field in columns]
    margin, gap, row_height = 12, 12, 13
    annotation_width = name_width + sum(widths) + gap * 4
    width = max(7.2 * 72, margin * 2 + annotation_width + 144)
    tree_width = width - margin * 2 - annotation_width
    name_x = margin + tree_width + gap
    column_x, cursor = [], name_x + name_width
    for column_width in widths:
        cursor += gap + column_width
        column_x.append(cursor)
    height = (len(tips) - 1) * row_height + 76
    header_y = height - 12
    maximum = max(distance.values())
    x = {node: margin + value / (maximum or 1) * tree_width for node, value in distance.items()}
    y = {node: 44 + i * row_height for i, node in enumerate(tips)}
    for node in phy.traverse("postorder"):
        if not node.is_leaf:
            y[node] = sum(y[c] for c in node.children) / len(node.children)

    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 8,
                         "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
                         "lines.solid_capstyle": "butt", "savefig.facecolor": "white"}):
        fig = plt.figure(figsize=(width / 72, height / 72))
        ax = fig.add_axes([0, 0, 1, 1], xlim=(0, width), ylim=(0, height))
        ax.set_axis_off()
        for node in phy.traverse():
            if node.children:
                ax.plot([x[node]] * 2, [min(y[c] for c in node.children), max(y[c] for c in node.children)],
                        color="#202020", linewidth=0.6)
            for child in node.children:
                ax.plot([x[node], x[child]], [y[child]] * 2, color="#202020", linewidth=0.6)
        for node in tips:
            row = rows[node.name]
            ax.plot([x[node] + 3, name_x - 4], [y[node]] * 2, color="#BBBBBB",
                    linewidth=0.4, linestyle=(0, (1.5, 2)), zorder=1)
            ax.scatter(x[node], y[node], color=colors[row[info["trait"]]], s=14,
                       edgecolors="white", linewidths=0.3, zorder=3)
            ax.text(name_x, y[node], names[node], fontproperties=italic, va="center", color="#202020")
            for (title, field), position in zip(columns, column_x):
                ax.text(position, y[node], row[field], ha="right", va="center",
                        color="#202020", fontweight="bold" if title == "Pair" else "normal")
        ax.text(name_x, header_y, "Species", va="center", color="#595959")
        for (title, _), position in zip(columns, column_x):
            ax.text(position, header_y, title, ha="right", va="center", color="#595959")
        ax.plot([name_x, width - margin], [header_y - 9] * 2, color="#CCCCCC", linewidth=0.4)
        labels = {"0": "C₃", "1": "C₄"} if info["trait"] == "C4" else {s: s for s in states}
        ax.legend(handles=[Line2D([], [], color=colors[s], marker="o", markersize=4,
                                linestyle="", label=labels[s]) for s in states],
                  loc="center left", bbox_to_anchor=(margin, header_y), bbox_transform=ax.transData,
                  ncol=len(states), frameon=False, borderaxespad=0,
                  handlelength=0.8, handletextpad=0.5, columnspacing=1.5)
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
            fig.savefig(Path(outdir) / f"summary_tree.{extension}")
        plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["tree", "metadata", "summary", "outdir"]:
        parser.add_argument("--" + name, required=True)
    plot(**vars(parser.parse_args()))
