"""Configuration shape validation, including partial command-line overrides."""
from pathlib import Path

import pytest
import yaml

from configuration import KEYS, validate_keys
from versioning import IMAGE_REPOSITORY, read_version, resolve_container_image

ROOT = Path(__file__).resolve().parents[1]


def test_config_keys_match_documented_defaults():
    config = yaml.safe_load((ROOT / "config/config.yaml").read_text())
    validate_keys(config)
    for section, keys in KEYS.items():
        values = config
        for part in section.split(".") if section else []:
            values = values[part]
        assert set(values) == set(keys.split()), section
    validate_keys(yaml.safe_load((ROOT / "config/pilot.yaml").read_text()))
    validate_keys({"phylogeny": {"dating": {"calibration_source": "file"}}})
    validate_keys({"odb": {"existing_results": "/data/snapshot"}})
    validate_keys({"run_name": "c4_run1"})


@pytest.mark.parametrize("config,path", [
    ({"analysis": "pilot"}, "analysis"),
    ({"tools": {}}, "tools"),
    ({"taxonomy": {"source": None}}, "taxonomy"),
    ({"taxonomy_audit": {"enabled": True}}, "taxonomy_audit"),
    ({"odb": {"mem_mb": 8000}}, "odb.mem_mb"),
    ({"odb": {"min_free_gb": 750}}, "odb.min_free_gb"),
    ({"odb": {"reference_min_free_gb": 200}}, "odb.reference_min_free_gb"),
    ({"odb": {"allow_nonlocal": True}}, "odb.allow_nonlocal"),
    ({"odb": {"keep_work": True}}, "odb.keep_work"),
    ({"odb": {"chunk_size": 100}}, "odb.chunk_size"),
    ({"odb": {"batch_size": 64}}, "odb.batch_size"),
    ({"phylogeny": {"sequence_mode": "protein"}}, "phylogeny.sequence_mode"),
    ({"phylogeny": {"sequence_dir": "/data/proteins"}}, "phylogeny.sequence_dir"),
    ({"phylogeny": {"sequence_suffix": ".fa"}}, "phylogeny.sequence_suffix"),
    ({"phylogeny": {"busco_full_suffix": ".tsv"}}, "phylogeny.busco_full_suffix"),
    ({"phylogeny": {"min_occupancy": 0.5}}, "phylogeny.min_occupancy"),
    ({"phylogeny": {"seed": 12345}}, "phylogeny.seed"),
    ({"phylogeny": {"dating": {"treepl": {}}}}, "phylogeny.dating.treepl"),
    ({"phylogeny": {"dating": {"lsd2": {}}}}, "phylogeny.dating.lsd2"),
    ({"phylogeny": {"dating": {"lsd2": {"variance": 1}}}}, "phylogeny.dating.lsd2"),
    ({"phylogeny": {"dating": {"lsd2": {"variance_parameter": None}}}}, "phylogeny.dating.lsd2"),
    ({"phylogeny": {"dating": {"lsd2": {"numsites": None}}}}, "phylogeny.dating.lsd2"),
    ({"phylogeny": {"dating": {"timetree": {"min_clade_taxa": 8}}}},
     "phylogeny.dating.timetree"),
    ({"phylogeny": {"dating": {"timetree": {"representatives": None}}}},
     "phylogeny.dating.timetree"),
    ({"phylogeny": {"dating": {"timetree": {"offline": False}}}},
     "phylogeny.dating.timetree"),
    ({"phylogeny": {"dating": {"timetree": {}}}}, "phylogeny.dating.timetree"),
    ({"phylogeny": {"dating": {"timetree": {"max_representatives": 64}}}}, "phylogeny.dating.timetree"),
    ({"phylogeny": {"dating": {"timetree": {"max_queries": 32}}}}, "phylogeny.dating.timetree"),
    ({"phylogeny": {"dating": {"timetree": {"min_studies": 5}}}}, "phylogeny.dating.timetree"),
    ({"kegg": {"thread": 2}}, "kegg.thread"),
    ({"phenoradar": {"alignments": True}}, "phenoradar"),
    ({"phenoradar": {}}, "phenoradar"),
    ({"species_filter": {}}, "species_filter"),
])
def test_unknown_settings_name_the_full_path(config, path):
    with pytest.raises(ValueError, match=f"unknown configuration settings: {path}"):
        validate_keys(config)


@pytest.mark.parametrize("path", [
    "odb.threads", "odb.mem_gb", "alignment.threads", "alignment.mem_gb",
    "kegg.threads", "kegg.mem_gb", "phylogeny.align_threads", "phylogeny.tree_threads",
    "phylogeny.astral_threads", "phylogeny.preparation_mem_gb", "phylogeny.alignment_mem_gb",
    "phylogeny.trimming_mem_gb", "phylogeny.tree_mem_gb", "phylogeny.astral_mem_gb",
    "phylogeny.dating.mem_gb", "taxonomy_check.mem_gb",
])
def test_resource_overrides_are_not_workflow_config(path):
    config = 8
    for key in reversed(path.split(".")):
        config = {key: config}
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


@pytest.mark.parametrize("seed", [1, 12345, 2147483647])
def test_global_seed_accepts_supported_integer_range(seed):
    validate_keys({"seed": seed})


@pytest.mark.parametrize("seed", [None, False, True, 0, -1, 1.5, "12345", [], {}, 2147483648])
def test_global_seed_rejects_unsupported_values(seed):
    with pytest.raises(ValueError, match="seed must be an integer between 1 and 2147483647"):
        validate_keys({"seed": seed})


@pytest.mark.parametrize("image", [
    None, False, 12, 0.1, [], {}, "", "  ", "latest", "v0.1.0", "0.1", "01.2.3",
    "1.2.3-rc1", "1.2.3+build", " 0.1.0", "0.1.0\n", "1.2.3\n2.0.0",
    "/data/custom image.sif", "images/release.sif",
    "docker://ghcr.io/mkrg01/phenoradar_prep:v0.1.0",
    "docker://ghcr.io/mkrg01/phenoradar_prep@sha256:" + "a" * 64,
    "docker://example/image:v1",
])
def test_invalid_container_images(image, tmp_path):
    with pytest.raises(ValueError, match="container_image must be auto or"):
        validate_keys({"container_image": image})
    for enabled in (True, False):
        with pytest.raises(ValueError, match="container_image must be auto or"):
            resolve_container_image(image, tmp_path, enabled=enabled)


def test_auto_image_follows_version_without_git(tmp_path):
    (tmp_path / "VERSION").write_text("0.2.0\n")
    assert resolve_container_image("auto", tmp_path, enabled=True) == f"docker://{IMAGE_REPOSITORY}:v0.2.0"
    (tmp_path / "VERSION").write_text("1.0.1\n")
    assert resolve_container_image("auto", tmp_path, enabled=True).endswith(":v1.0.1")
    assert not (tmp_path / ".git").exists()


@pytest.mark.parametrize("version", ["0.1.0", "1.2.3", "10.20.30"])
def test_explicit_version_does_not_need_version_file(tmp_path, version):
    validate_keys(yaml.safe_load(f'container_image: "{version}"'))
    assert resolve_container_image(version, tmp_path, enabled=True) == f"docker://{IMAGE_REPOSITORY}:v{version}"


def test_explicit_version_overrides_checkout_version(tmp_path):
    (tmp_path / "VERSION").write_text("0.2.0\n")
    assert resolve_container_image("0.1.0", tmp_path, enabled=True) == f"docker://{IMAGE_REPOSITORY}:v0.1.0"


@pytest.mark.parametrize("image", ["auto", "0.1.0"])
def test_native_execution_does_not_need_a_container_or_version_file(tmp_path, image):
    assert resolve_container_image(image, tmp_path, enabled=False) is None


@pytest.mark.parametrize("version", ["", "v0.2.0", "0.2", "01.2.3", "1.2.3-rc1", "1.2.3\n2.0.0"])
def test_invalid_release_versions_fail_clearly(tmp_path, version):
    (tmp_path / "VERSION").write_text(version)
    with pytest.raises(ValueError, match="VERSION must contain"):
        read_version(tmp_path)


def test_missing_version_has_archive_instructions(tmp_path):
    with pytest.raises(ValueError, match="complete workflow checkout or release archive"):
        resolve_container_image("auto", tmp_path, enabled=True)
