"""Public metadata uses canonical names while earlier saved inputs remain readable."""
import json

import pytest

from common import file_record, read_tsv, write_json, write_tsv
from contrast_metadata import canonical_row
from consumer_metadata import create, validate
from test_consumer_metadata import cohort


@pytest.mark.parametrize("suffix", ["", "_all", "_phenotyped", "_representatives"])
def test_legacy_block_columns_accept_equal_aliases_and_reject_conflicts(suffix):
    old, new = "expanded_pair_clade_id" + suffix, "phylogenetic_block_id" + suffix
    assert canonical_row({old: "1", new: "1"}) == {new: "1"}
    with pytest.raises(ValueError, match="conflicting metadata aliases"):
        canonical_row({old: "1", new: "2"})


def test_trait_clade_aliases_apply_only_to_clade_metadata():
    old = {"group": "2", "n_species_in_group": "3", "n_samples_in_group": "4"}
    assert canonical_row(old) == old
    assert canonical_row(old, trait_clades=True) == {
        "trait_clade_id": "2", "n_species_in_trait_clade": "3", "n_samples_in_trait_clade": "4"}
    with pytest.raises(ValueError, match="conflicting metadata aliases"):
        canonical_row({"group": "1", "trait_clade_id": "2"}, trait_clades=True)


def test_legacy_inputs_emit_new_names_and_old_outputs_remain_verifiable(cohort):
    kwargs, paths = cohort
    rows = read_tsv(paths["all"])
    for row in rows:
        row["expanded_pair_clade_id"] = "2" if row["habit"] else ""
    write_tsv(paths["all"], list(rows[0]), rows)
    expected = create(**kwargs, pair_tables={"all": paths["all"]})
    assert "phylogenetic_block_id" in expected[0]
    assert "expanded_pair_clade_id" not in expected[0]
    # Reconstruct an earlier schema-2 output, retaining its own recorded checksum.
    columns = ["expanded_pair_clade_id" if k == "phylogenetic_block_id" else k for k in expected[0]]
    legacy = [{("expanded_pair_clade_id" if k == "phylogenetic_block_id" else k): v
               for k, v in row.items()} for row in expected]
    write_tsv(kwargs["output"], columns, legacy)
    info_path = kwargs["output"].with_suffix(".json")
    info = json.loads(info_path.read_text())
    info.update(schema_version=2, output=dict(file_record(kwargs["output"]), path="metadata/species_metadata.tsv"))
    write_json(info_path, info)
    before = kwargs["output"].read_bytes()
    assert validate(kwargs["output"], "habit") == expected
    assert kwargs["output"].read_bytes() == before
