"""Trait legends preserve labels, colors, and layout across contrast subsets."""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
from matplotlib.figure import Figure
import pytest

from common import write_json, write_tsv
from plot_contrast_tree import paired_clades, plot, tip_identity


@pytest.mark.parametrize("trait,states,state_order,expected_colors", [
    ("C4", ["0", "1"], None, ["#595959", "#009e73"]),
    ("carnivory", ["0", "1"], None, ["#595959", "#009e73"]),
    ("carnivory", ["1"], None, ["#009e73"]),
    ("habit", ["herbaceous", "woody"], ["herbaceous", "woody"], ["#595959", "#009e73"]),
    ("habit", ["woody"], ["herbaceous", "woody"], ["#009e73"]),
    ("A long descriptive trait name " * 4,
     ["absence of the focal phenotype " * 4, "presence of the focal phenotype " * 4],
     None, ["#595959", "#009e73"]),
])
def test_trait_legend_matches_tip_colors_and_fits_above_labels(
        tmp_path, monkeypatch, trait, states, state_order, expected_colors):
    tree = tmp_path / "tree.nwk"
    tree.write_text("((A:1,B:1):1,(C:1,D:1):1);")
    metadata = tmp_path / "metadata.tsv"
    write_tsv(metadata, ["species", "is_representative", trait, "n_species_in_group", "group", "contrast_pair_id"],
              [dict(species=name, is_representative="1", **{trait: states[i % len(states)]},
                    n_species_in_group="1", group=str(i + 1), contrast_pair_id="")
               for i, name in enumerate("ABCD")])
    summary = tmp_path / "summary.json"
    info = {"trait": trait, "summary_species": 4}
    if state_order is not None:
        info["state_order"] = state_order
    write_json(summary, info)
    original = {path: path.read_bytes() for path in [tree, metadata, summary]}
    checked = []
    close = plt.close

    def check_and_close(figure=None):
        try:
            if isinstance(figure, Figure):
                figure.canvas.draw()
                ax = figure.axes[0]
                legend = ax.get_legend()
                assert legend.get_title().get_text() == trait
                assert [text.get_text() for text in legend.get_texts()] == states
                assert [to_hex(handle.get_color()) for handle in legend.legend_handles] == expected_colors
                expected = dict(zip(states, expected_colors))
                assert [to_hex(tip.get_facecolors()[0]) for tip in ax.collections] == [
                    expected[states[i % len(states)]] for i in range(4)]
                bounds = legend.get_window_extent()
                assert bounds.x0 >= 0 and bounds.x1 <= figure.bbox.x1
                assert bounds.y1 <= figure.bbox.y1
                names = [text for text in ax.texts if text.get_text() in "ABCD"]
                assert bounds.y0 > max(text.get_window_extent().y1 for text in names)
                for text in ax.texts:
                    assert text.get_window_extent().x1 <= figure.bbox.x1
                checked.append(True)
        finally:
            close(figure)

    monkeypatch.setattr(plt, "close", check_and_close)
    plot(tree, metadata, summary, tmp_path)
    assert checked == [True, True]
    for stem in ["summary_tree", "summary_tree_id"]:
        assert (tmp_path / f"{stem}.pdf").read_bytes().startswith(b"%PDF")
    svg = (tmp_path / "summary_tree.svg").read_text()
    assert "C₃" not in svg and "C₄" not in svg
    assert {path: path.read_bytes() for path in original} == original


def test_empty_tree_reports_configured_trait(tmp_path):
    tree = tmp_path / "tree.nwk"
    tree.write_text("")
    metadata = tmp_path / "metadata.tsv"
    write_tsv(metadata, ["species", "is_representative", "carnivory"], [])
    summary = tmp_path / "summary.json"
    write_json(summary, {"trait": "carnivory", "summary_species": 0, "state_order": ["0", "1"]})
    plot(tree, metadata, summary, tmp_path)
    for stem in ["summary_tree", "summary_tree_id"]:
        assert "No species with observed carnivory" in (tmp_path / f"{stem}.svg").read_text()
        assert (tmp_path / f"{stem}.pdf").read_bytes().startswith(b"%PDF")
    assert json.loads(summary.read_text())["state_order"] == ["0", "1"]


def test_publication_and_id_labels_keep_accessions_upright_and_color_only_pair_clades(tmp_path, monkeypatch):
    tree = tmp_path / "tree.nwk"
    tree.write_text("((sample_A:1,sample_B:1):1,(sample_C:1,sample_D:1):1);")
    metadata, samples, summary = [tmp_path / name for name in ["metadata.tsv", "samples.tsv", "summary.json"]]
    names = ["Cyperus fuscus", "Cyperus esculentus", "Panicum miliaceum", "Panicum miliaceum"]
    write_tsv(samples, ["species", "scientific_name", "run"],
              [dict(species=f"sample_{tip}", scientific_name=name, run=f"SRR_{tip}_01")
               for tip, name in zip("ABCD", names)])
    write_tsv(metadata, ["species", "is_representative", "C4", "n_species_in_group", "group", "contrast_pair_id"],
              [dict(species=f"sample_{tip}", is_representative="1", C4=str(i % 2),
                    n_species_in_group="7" if i == 0 else "1", group=str(i + 1),
                    contrast_pair_id="9" if i < 2 else "") for i, tip in enumerate("ABCD")])
    write_json(summary, {"trait": "C4", "summary_species": 3})
    original = {path: path.read_bytes() for path in [tree, metadata, samples, summary]}
    figures, geometry = [], []
    close = plt.close

    def check_and_close(figure=None):
        try:
            if not isinstance(figure, Figure):
                return
            figure.canvas.draw()
            ax = figure.axes[0]
            texts = ax.texts
            assert [text.get_text() for text in texts if text.get_fontstyle() == "italic"] == names
            suffixes = [text for text in texts if text.get_text().startswith(" (n = ")]
            assert len(suffixes) == 4 and all(text.get_fontstyle() == "normal" for text in suffixes)
            assert not {"Species", "n", "Group", "Pair"} & {text.get_text() for text in texts}
            if not figures:
                assert [text.get_text() for text in suffixes] == [" (n = 7)", *(" (n = 1)" for _ in range(3))]
                assert all("SRR" not in text.get_text() for text in texts)
            else:
                assert suffixes[0].get_text() == " (n = 7; run = SRR_A_01; group = 1; pair = 9)"
                assert suffixes[2].get_text() == " (n = 1; run = SRR_C_01; group = 3; pair = -)"
                assert suffixes[3].get_text() == " (n = 1; run = SRR_D_01; group = 4; pair = -)"
            red = [line for line in ax.lines if to_hex(line.get_color()) == "#c62828"]
            assert len(red) == 3  # Two descendant edges and their joining vertical line.
            assert all(max(line.get_ydata()) < 70 for line in red)
            for text in texts:
                bounds = text.get_window_extent()
                assert bounds.x0 >= 0 and bounds.x1 <= figure.bbox.x1
                assert bounds.y0 >= 0 and bounds.y1 <= figure.bbox.y1
            geometry.append([(list(line.get_xdata()), list(line.get_ydata())) for line in ax.lines])
            figures.append(True)
        finally:
            close(figure)

    monkeypatch.setattr(plt, "close", check_and_close)
    plot(tree, metadata, summary, tmp_path, samples=samples)
    assert figures == [True, True] and geometry[0] == geometry[1]
    assert {path: path.read_bytes() for path in original} == original


@pytest.mark.parametrize("members,message", [(["A"], "exactly two"), (["A", "B", "C"], "exactly two"),
                                           (["A", "C"], "not a summary-tree clade")])
def test_invalid_pair_does_not_color_unrelated_tips(members, message):
    from ete4 import Tree
    rows = {name: {"contrast_pair_id": "1" if name in members else ""} for name in "ABCD"}
    with pytest.raises(ValueError, match=message):
        paired_clades(Tree("((A,B),(C,D));", parser=9), rows)


def test_legacy_labels_do_not_guess_species_or_split_run_ids():
    assert tip_identity({"species": "opaque_SRR_01"}, {}) == ("opaque_SRR_01", False, "")
    row = {"species": "Cyperus_fuscus_SRR_01", "species_id": "Cyperus_fuscus"}
    assert tip_identity(row, {}) == ("Cyperus fuscus", True, "")
    rows = [{**row, "scientific_name": "Cyperus fuscus", "run": run} for run in ["SRR_01", "SRR_02"]]
    assert tip_identity(row, {row["species"]: rows}) == ("Cyperus fuscus", True, "SRR_01, SRR_02")
