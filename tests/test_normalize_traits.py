"""Only binary traits or missing values enter an analysis."""
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
from filter_species import discover
from normalize_traits import normalize
from phase_config import read_yaml

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def trait_inputs(tmp_path):
    samples, traits = tmp_path / "samples.tsv", tmp_path / "traits.tsv"
    rows = [dict(species=f"Plant_{name}_R{i}", species_id=f"Plant_{name}")
            for i, name in enumerate("AABCDEFG")]
    write_tsv(samples, list(rows[0]), rows)
    return samples, traits, tmp_path / "normalized.tsv", tmp_path / "report.json"


@pytest.mark.parametrize("trait", ["carnivory", "C4"])
def test_normalization_preserves_binary_states_and_missing_values(trait_inputs, trait):
    samples, traits, output, report = trait_inputs
    write_tsv(traits, ["species", trait],
              [{"species": f"Plant {name}", trait: value}
               for name, value in zip("ABCDEFZ", ["0", "1", "NA", "NaN", "nan", "", "unused state"])])
    original = traits.read_bytes()
    normalize(samples, output, report, trait, traits)
    assert [row[trait] for row in read_tsv(output)] == ["0", "1", "", "", "", "", ""]
    info = json.loads(report.read_text())
    assert info["biological_species"] == 7 and info["missing_trait_species"] == 5
    assert info["observed_states"] == ["0", "1"]
    assert traits.read_bytes() == original


@pytest.mark.parametrize("source", [True, False])
def test_single_binary_state_and_absent_source(trait_inputs, source):
    samples, traits, output, report = trait_inputs
    write_tsv(traits, ["species", "C4"], [{"species": "Plant A", "C4": "1"}])
    normalize(samples, output, report, traits=traits if source else None)
    assert read_tsv(output)[0]["C4"] == ("1" if source else "")
    assert json.loads(report.read_text())["observed_states"] == (["1"] if source else [])


@pytest.mark.parametrize("invalid_state", ["absent", "2", "0.0"])
def test_nonbinary_states_fail_before_replacing_outputs(trait_inputs, invalid_state):
    samples, traits, output, report = trait_inputs
    write_tsv(traits, ["species", "habit"], [{"species": "Plant A", "habit": invalid_state}])
    output.write_text("previous output")
    with pytest.raises(ValueError, match="must be 0, 1 or missing"):
        normalize(samples, output, report, "habit", traits)
    assert output.read_text() == "previous output" and not report.exists()


def test_legacy_config_only_accepts_null_encoding(tmp_path):
    path = tmp_path / "analysis.yaml"
    path.write_text("trait_encoding: null\n")
    validate_keys({"trait_encoding": None})
    assert read_yaml(path) == {}

    path.write_text("trait_encoding: {absent: 0, present: 1}\n")
    with pytest.raises(ValueError, match="trait_encoding has been removed"):
        read_yaml(path)
    with pytest.raises(ValueError, match="trait_encoding has been removed"):
        validate_keys({"trait_encoding": {"absent": 0, "present": 1}})


def test_normalization_cannot_overwrite_its_inputs(trait_inputs):
    samples, traits, output, report = trait_inputs
    write_tsv(traits, ["species", "C4"], [{"species": "Plant A", "C4": "1"}])
    original = traits.read_bytes()
    with pytest.raises(ValueError, match="distinct"):
        normalize(samples, traits, report, traits=traits)
    assert traits.read_bytes() == original


def test_workflow_normalizes_binary_traits_before_consumer_metadata(
        tiny_inputs, workflow_project, command_environment, seed_taxonomy):
    snakemake = os.environ.get("SNAKEMAKE_BIN") or shutil.which("snakemake")
    if not snakemake:
        pytest.skip("Snakemake required")
    source = Path(tiny_inputs["metadata"]).parent
    write_tsv(source / "species_trait.tsv", ["species", "habit"],
              [dict(species=name, habit=value) for name, value in
               [("Alpha plant", "0"), ("Beta sp-X", "1"), ("Gamma plant", "unselected state")]])
    seed_taxonomy(tiny_inputs["taxonomy_db"])
    config = workflow_project / "override.yaml"
    config.write_text(yaml.safe_dump(dict(run_name="test", input_root=str(source), trait="habit")))
    target = "results/test/metadata/species_metadata.tsv"
    result = subprocess.run([snakemake, "--snakefile", str(ROOT / "workflow/Snakefile"),
                             "--configfile", str(config), "--cores", "1", "--", target],
                            cwd=workflow_project, env=command_environment({"python": sys.executable}),
                            text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    assert result.returncode == 0, result.stdout
    assert [row["habit"] for row in read_tsv(workflow_project / target)] == ["0", "1"]
    analysis = workflow_project / "results/test"
    report = analysis / "metadata/species_trait.json"
    assert json.loads(report.read_text())["observed_states"] == ["0", "1"]
    assert not (analysis / "metadata/trait_encoding.json").exists()
    assert discover(analysis)["sections"]["traits"]["path"] == str(analysis / "metadata/species_trait.tsv")
