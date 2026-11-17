import pytest

from common import write_tsv
from date_phylogeny import ROOT_MAX_MA, tree_constraints, validate_dated
from infer_phylogeny import read_tree


def test_redundant_ancestor_minimum_does_not_block_treepl(tmp_path):
    tree_file = tmp_path / "tree.nwk"
    tree_file.write_text("(((A:90,B:90)Child:10,C:100)Parent:40,D:140)Root;")
    tree = read_tree(tree_file)
    table = tmp_path / "minima.tsv"
    fields = ["left_species", "right_species", "min_age_ma", "fossil_ids"]
    minima = [
        dict(left_species="A", right_species="B", min_age_ma=85, fossil_ids="10"),
        dict(left_species="A", right_species="C", min_age_ma=50, fossil_ids="11"),
    ]
    write_tsv(table, fields, minima)
    constraints, config, redundant = tree_constraints(tree, table)
    assert redundant == 1
    assert config == ["mrca = Child A B", "min = Child 85", "mrca = Root A D", "max = Root 154"]
    assert {row["node"] for row in constraints} == {"Child", "Parent", "Root"}
    assert validate_dated(tree_file, tree, constraints)[2] == 140

    too_old = tmp_path / "too_old.nwk"
    too_old.write_text("(((A:90,B:90)Child:10,C:100)Parent:55,D:155)Root;")
    with pytest.raises(ValueError, match="beyond native rounding at Root"):
        validate_dated(too_old, tree, constraints)

    minima[0]["min_age_ma"] = ROOT_MAX_MA + 1
    write_tsv(table, fields, minima)
    with pytest.raises(ValueError, match="incompatible AngioCal minimum"):
        tree_constraints(tree, table)


def test_timetree_bounds_replace_angiosperm_root_cap_and_detect_conflict(tmp_path):
    tree_file = tmp_path / "tree.nwk"
    tree_file.write_text("((A:90,B:90)Child:95,(C:80,D:80):105)Root;")
    tree = read_tree(tree_file)
    table = tmp_path / "bounds.tsv"
    rows = [dict(left_species="A", right_species="C", min_age_ma=180, max_age_ma=190),
            dict(left_species="A", right_species="B", min_age_ma=85, max_age_ma=95)]
    write_tsv(table, list(rows[0]), rows)

    constraints, config, _ = tree_constraints(tree, table, "timetree")
    assert "max = Root 190" in config and "min = Child 85" in config and "max = Child 95" in config
    assert validate_dated(tree_file, tree, constraints)[2] == 185

    rows[1].update(min_age_ma=191, max_age_ma=200)
    write_tsv(table, list(rows[0]), rows)
    with pytest.raises(ValueError, match="ancestor A / C has maximum 190 Ma, but descendant A / B has minimum 191 Ma"):
        tree_constraints(tree, table, "timetree")
