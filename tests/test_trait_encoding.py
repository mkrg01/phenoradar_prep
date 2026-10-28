"""Text phenotypes have an explicit, auditable binary encoding before analysis."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest
import yaml

from common import read_tsv, write_tsv
from configuration import validate_keys
from trait_encoding import normalize, read_labels

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def trait_inputs(tmp_path):
    samples, traits = tmp_path / "samples.tsv", tmp_path / "traits.tsv"
    rows = [dict(species=f"Plant_{name}_R{i}", species_id=f"Plant_{name}")
            for i, name in enumerate("AABCDEFG")]
    write_tsv(samples, list(rows[0]), rows)
    return samples, traits, tmp_path / "normalized.tsv", tmp_path / "encoding.json"


@pytest.mark.parametrize("trait", ["carnivory", "C4"])
def test_text_encoding_preserves_labels_missing_values_and_sample_identity(trait_inputs, trait):
    samples, traits, output, report = trait_inputs
    write_tsv(traits, ["species", trait],
              [{"species": f"Plant {name}", trait: value}
               for name, value in zip("ABCDEFZ", ["absent", "present", "NA", "NaN", "nan", "", "unused state"])])
    original = traits.read_bytes()
    normalize(samples, output, report, trait, traits, {"absent": 0, "present": 1})
    assert [row[trait] for row in read_tsv(output)] == ["0", "1", "", "", "", "", ""]
    assert read_labels(report, trait) == {"0": "absent", "1": "present"}
    info = json.loads(report.read_text())
    assert info["biological_species"] == 7 and info["missing_trait_species"] == 5
    assert info["observed_source_states"] == ["absent", "present"]
    assert traits.read_bytes() == original


def test_mapping_direction_is_explicit_even_with_one_observed_state(trait_inputs):
    samples, traits, output, report = trait_inputs
    write_tsv(traits, ["species", "habit"], [{"species": "Plant A", "habit": "herbaceous"}])
    normalize(samples, output, report, "habit", traits, {"woody": 0, "herbaceous": 1})
    assert read_tsv(output)[0]["habit"] == "1"
    assert read_labels(report, "habit") == {"0": "woody", "1": "herbaceous"}


@pytest.mark.parametrize("source", [True, False])
def test_numeric_traits_and_absent_source_need_no_mapping(trait_inputs, source):
    samples, traits, output, report = trait_inputs
    write_tsv(traits, ["species", "C4"], [{"species": "Plant A", "C4": "1"}])
    normalize(samples, output, report, "C4", traits if source else None)
    assert read_tsv(output)[0]["C4"] == ("1" if source else "")
    assert read_labels(report, "C4") == {"0": "0", "1": "1"}


@pytest.mark.parametrize("mapping,message", [(None, "trait_encoding is required"),
                                             ({"absent": 0, "present": 1}, "absent from trait_encoding")])
def test_unmapped_labels_fail_before_replacing_outputs(trait_inputs, mapping, message):
    samples, traits, output, report = trait_inputs
    write_tsv(traits, ["species", "habit"], [{"species": "Plant A", "habit": "unexpected"}])
    output.write_text("previous output")
    with pytest.raises(ValueError, match=message):
        normalize(samples, output, report, "habit", traits, mapping)
    assert output.read_text() == "previous output" and not report.exists()


@pytest.mark.parametrize("mapping", [False, [], {}, {"absent": 0}, {"absent": 0, "present": 0},
                                    {"absent": 0, "present": True}, {"absent": 0, "present": 1.0},
                                    {"NA": 0, "present": 1}, {" absent": 0, "present": 1}, {0: 0, 1: 1}])
def test_invalid_mapping_is_rejected_during_configuration(mapping):
    with pytest.raises(ValueError, match="trait_encoding"):
        validate_keys({"trait_encoding": mapping})


def test_encoding_record_must_match_requested_trait(trait_inputs):
    samples, traits, output, report = trait_inputs
    normalize(samples, output, report, "habit")
    with pytest.raises(ValueError, match="differs"):
        read_labels(report, "C4")


def test_normalization_cannot_overwrite_its_inputs(trait_inputs):
    samples, traits, output, report = trait_inputs
    write_tsv(traits, ["species", "C4"], [{"species": "Plant A", "C4": "1"}])
    original = traits.read_bytes()
    with pytest.raises(ValueError, match="distinct"):
        normalize(samples, traits, report, traits=traits)
    assert traits.read_bytes() == original


def test_workflow_normalizes_text_traits_before_consumer_metadata(
        tiny_inputs, workflow_project, command_environment, seed_taxonomy):
    snakemake = os.environ.get("SNAKEMAKE_BIN") or shutil.which("snakemake")
    if not snakemake:
        pytest.skip("Snakemake required")
    source = Path(tiny_inputs["metadata"]).parent
    write_tsv(source / "species_trait.tsv", ["species", "habit"],
              [dict(species=name, habit=value) for name, value in
               [("Alpha plant", "absent"), ("Beta sp-X", "present"), ("Gamma plant", "unselected state")]])
    seed_taxonomy(tiny_inputs["taxonomy_db"])
    config = workflow_project / "override.yaml"
    config.write_text(yaml.safe_dump(dict(run_name="test", input_root=str(source), trait="habit",
                                         trait_encoding={"absent": 0, "present": 1})))
    target = "results/test/metadata/species_metadata.tsv"
    result = subprocess.run([snakemake, "--snakefile", str(ROOT / "workflow/Snakefile"),
                             "--configfile", str(config), "--cores", "1", "--", target],
                            cwd=workflow_project, env=command_environment({"python": sys.executable}),
                            text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    assert result.returncode == 0, result.stdout
    assert [row["habit"] for row in read_tsv(workflow_project / target)] == ["0", "1"]
    encoding = workflow_project / "results/test/metadata/trait_encoding.json"
    assert read_labels(encoding, "habit") == {"0": "absent", "1": "present"}
