"""Configuration shape validation, including partial command-line overrides."""
import json
from pathlib import Path

import pytest
import yaml

from configuration import KEYS, validate_keys
from date_phylogeny import validate_settings as validate_treepl_settings
from versioning import IMAGE_REPOSITORY, read_version, resolve_container_image

ROOT = Path(__file__).resolve().parents[1]


def test_config_and_optional_tool_defaults_cover_supported_keys():
    config = yaml.safe_load((ROOT / "config/config.yaml").read_text())
    validate_keys(config)
    assert "treepl" not in config["phylogeny"]["dating"]
    config["phylogeny"]["dating"]["treepl"] = validate_treepl_settings({})
    for section, keys in KEYS.items():
        values = config
        for part in section.split(".") if section else []:
            values = values[part]
        assert set(values) == set(keys.split()), section
    validate_keys(yaml.safe_load((ROOT / "config/pilot.yaml").read_text()))
    validate_keys({"phylogeny": {"dating": {"calibration_source": "file"}}})
    validate_keys({"run_name": "c4_run1"})


@pytest.mark.parametrize("config,path", [
    ({"unknown": {}}, "unknown"),
    ({"container_image": "auto"}, "container_image"),
    ({"inputs": {"metadata": "elsewhere/metadata.tsv"}}, "inputs"),
    ({"phylogeny": {"busco_full_dir": "elsewhere/busco"}}, "phylogeny.busco_full_dir"),
    ({"phylogeny": {"dating": {"calibrations": "elsewhere/ages.tsv"}}}, "phylogeny.dating.calibrations"),
    ({"phylogeny": {"enable": True}}, "phylogeny.enable"),
    ({"phylogeny": {"dating": {"enable": True}}}, "phylogeny.dating.enable"),
])
def test_unknown_settings_name_the_full_path(config, path):
    with pytest.raises(ValueError, match=f"unknown configuration settings: {path}"):
        validate_keys(config)


@pytest.mark.parametrize("config,section", [
    (None, "configuration"),
    ({"odb": []}, "odb"),
    ({"phylogeny": None}, "phylogeny"),
    ({"phylogeny": {"dating": False}}, "phylogeny.dating"),
])
def test_sections_require_mappings(config, section):
    with pytest.raises(ValueError, match=f"{section} must be a mapping"):
        validate_keys(config)


@pytest.mark.parametrize("value", [True, False])
def test_species_list_is_an_explicit_opt_in(value):
    validate_keys({"selection": {"species_list": value}})


@pytest.mark.parametrize("value", [None, "input/species_list.txt", "false", 0, 1, [], {}])
def test_species_list_rejects_paths_and_nonbooleans(value):
    with pytest.raises(ValueError, match="selection.species_list must be true or false"):
        validate_keys({"selection": {"species_list": value}})


@pytest.mark.parametrize("seed", [1, 12345, 2147483647])
def test_global_seed_accepts_supported_integer_range(seed):
    validate_keys({"seed": seed})


@pytest.mark.parametrize("seed", [None, False, True, 0, -1, 1.5, "12345", [], {}, 2147483648])
def test_global_seed_rejects_unsupported_values(seed):
    with pytest.raises(ValueError, match="seed must be an integer between 1 and 2147483647"):
        validate_keys({"seed": seed})


def test_image_follows_version_file(tmp_path):
    (tmp_path / "VERSION").write_text("0.2.0\n")
    assert resolve_container_image(tmp_path) == f"docker://{IMAGE_REPOSITORY}:v0.2.0"
    (tmp_path / "VERSION").write_text("1.0.1\n")
    assert resolve_container_image(tmp_path) == f"docker://{IMAGE_REPOSITORY}:v1.0.1"


def test_run_record_keeps_resolved_image_separate_from_configuration(tmp_path):
    from record_run import record
    (tmp_path / "VERSION").write_text("0.2.0\n")
    workflow = tmp_path / "workflow"
    workflow.mkdir()
    selection = tmp_path / "selection.json"
    selection.write_text("{}\n")
    output = tmp_path / "run.json"
    image = resolve_container_image(tmp_path)
    record('{"run_name": "test"}', selection, workflow, output, container_image=image)
    report = json.loads(output.read_text())
    assert report["container_image"] == image
    assert report["config"] == {"run_name": "test"}
    validate_keys(report["config"])


@pytest.mark.parametrize("version", ["", "v0.2.0", "0.2", "01.2.3", "1.2.3-rc1", "1.2.3\n2.0.0"])
def test_invalid_release_versions_fail_clearly(tmp_path, version):
    (tmp_path / "VERSION").write_text(version)
    with pytest.raises(ValueError, match="VERSION must contain"):
        read_version(tmp_path)


def test_missing_version_has_archive_instructions(tmp_path):
    with pytest.raises(ValueError, match="complete workflow checkout or release archive"):
        resolve_container_image(tmp_path)


@pytest.mark.parametrize("value", [None, "resources/odb_existing/tlight", "/data/snapshot"])
def test_existing_odb_accepts_null_or_directory(value):
    validate_keys({"odb": {"existing_results": value}})


@pytest.mark.parametrize("value", ["", " ", False, True, 0, 1, [], {}])
def test_existing_odb_rejects_invalid_settings(value):
    with pytest.raises(ValueError, match="odb.existing_results must be null or a snapshot directory"):
        validate_keys({"odb": {"existing_results": value}})


def analysis_config(trees=(), **enabled):
    return {"trait": "C4", "phylogeny": {
        "trees": list(trees),
        **{name: {"enabled": enabled.get(name, False)}
           for name in ["contrast_pairs", "dating", "taxonomy_check"]},
    }}


@pytest.mark.parametrize("trees", [[], ["all"], ["phenotyped"], ["representatives"],
                                   ["all", "phenotyped", "representatives"]])
def test_tree_selection_is_the_only_inference_switch(trees):
    from configuration import validate_analysis
    config = analysis_config(trees)
    validate_keys(config)
    validate_analysis(config, ["all"])
    if trees:
        validate_analysis(config, ["phylogeny", "phylogeny_prepare"])


@pytest.mark.parametrize("trees", [None, False, "all", ["unknown"], ["all", "all"], [["all"]]])
def test_invalid_tree_selection_is_rejected(trees):
    with pytest.raises(ValueError, match="phylogeny.trees must be a list"):
        validate_keys({"phylogeny": {"trees": trees}})


@pytest.mark.parametrize("trait", [None, False, "", "  ", "species", "contrast_pair_id"])
def test_invalid_shared_trait_is_rejected(trait):
    with pytest.raises(ValueError, match="trait"):
        validate_keys({"trait": trait})


@pytest.mark.parametrize("section", ["contrast_pairs", "dating", "taxonomy_check"])
def test_postprocessing_requires_a_selected_tree(section):
    from configuration import validate_analysis
    with pytest.raises(ValueError, match=f"phylogeny.{section}.enabled requires a nonempty"):
        validate_analysis(analysis_config(**{section: True}))


@pytest.mark.parametrize("section", ["dating", "taxonomy_check"])
def test_unsupported_representative_postprocessing_is_explicit(section):
    from configuration import validate_analysis
    with pytest.raises(ValueError, match=f"phylogeny.{section} does not yet support representatives"):
        validate_analysis(analysis_config(["representatives"], **{section: True}))


@pytest.mark.parametrize("target,section", [("contrast_pairs", "contrast_pairs"), ("timetree", "dating"),
    ("phylogeny_calibrations", "dating"), ("taxonomy_check", "taxonomy_check")])
def test_explicit_targets_cannot_enable_disabled_postprocessing(target, section):
    from configuration import validate_analysis
    with pytest.raises(ValueError, match=f"phylogeny.{section}.enabled: true"):
        validate_analysis(analysis_config(["all"]), [target])
    validate_analysis(analysis_config(["all"], **{section: True}), [target])


@pytest.mark.parametrize("trees,target,expected", [
    ([], "phylogeny", "requires a nonempty phylogeny.trees"),
    (["all"], "contrast_pairs", "phylogeny.contrast_pairs.enabled: true"),
])
def test_snakefile_rejects_conflicting_targets_before_execution(workflow_project, trees, target, expected):
    import os
    import shutil
    import subprocess
    snakemake = os.environ.get("SNAKEMAKE_BIN") or shutil.which("snakemake")
    if not snakemake:
        pytest.skip("Snakemake required")
    override = workflow_project / "override.yaml"
    override.write_text(yaml.safe_dump(analysis_config(trees)))
    result = subprocess.run([snakemake, "--snakefile", str(ROOT / "workflow/Snakefile"),
                             "--configfile", str(override), "--dry-run", "--", target],
                            cwd=workflow_project, text=True, capture_output=True,
                            env={**os.environ, "XDG_CACHE_HOME": str(workflow_project / "cache")})
    output = result.stdout + result.stderr
    assert result.returncode != 0 and expected in output, output
    assert not (workflow_project / "results").exists()


@pytest.mark.parametrize("trees,target", [
    ([], "phylogeny/all/species_tree.nwk"),
    (["representatives"], "phylogeny/all/species_tree.nwk"),
    (["all"], "phylogeny/representatives/species_tree.nwk"),
    (["all"], "phylogeny/all/contrast/contrast_pairs.tsv"),
    (["representatives"], "phylogeny/representatives/contrast/contrast_pairs.tsv"),
    (["all"], "phylogeny/all/dating/species_tree.dated.nwk"),
    (["all"], "phylogeny/all/taxonomy_check"),
])
def test_output_targets_cannot_schedule_unselected_trees_or_disabled_steps(workflow_project, trees, target):
    import os
    import shutil
    import subprocess
    snakemake = os.environ.get("SNAKEMAKE_BIN") or shutil.which("snakemake")
    if not snakemake:
        pytest.skip("Snakemake required")
    override = workflow_project / "override.yaml"
    override.write_text(yaml.safe_dump({"run_name": "test", **analysis_config(trees)}))
    result = subprocess.run([snakemake, "--snakefile", str(ROOT / "workflow/Snakefile"),
                             "--configfile", str(override), "--dry-run", "--", f"results/test/{target}"],
                            cwd=workflow_project, text=True, capture_output=True,
                            env={**os.environ, "XDG_CACHE_HOME": str(workflow_project / "cache")})
    output = result.stdout + result.stderr
    assert result.returncode != 0 and "MissingRuleException" in output, output
    assert not (workflow_project / "results").exists()
