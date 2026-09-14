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
    validate_keys({"phylogeny": {"dating": {"lsd2": {"variance": 0}}}})
    validate_keys({"odb": {"existing_results": "/data/snapshot"}})
    validate_keys({"run_name": "c4_run1"})


@pytest.mark.parametrize("config,path", [
    ({"analysis": "pilot"}, "analysis"),
    ({"tools": {}}, "tools"),
    ({"odb": {"mem_mb": 8000}}, "odb.mem_mb"),
    ({"phylogeny": {"min_occupancy": 0.5}}, "phylogeny.min_occupancy"),
    ({"phylogeny": {"dating": {"treepl": {}}}}, "phylogeny.dating.treepl"),
    ({"phylogeny": {"dating": {"timetree": {"min_clade_taxa": 8}}}},
     "phylogeny.dating.timetree.min_clade_taxa"),
    ({"kegg": {"thread": 2}}, "kegg.thread"),
    ({"phenoradar": {"alginments": True}}, "phenoradar.alginments"),
    ({"species_filter": {}}, "species_filter"),
])
def test_unknown_settings_name_the_full_path(config, path):
    with pytest.raises(ValueError, match=f"unknown configuration settings: {path}"):
        validate_keys(config)


@pytest.mark.parametrize("config,section", [
    (None, "configuration"),
    ({"odb": []}, "odb"),
    ({"phylogeny": None}, "phylogeny"),
    ({"phylogeny": {"dating": False}}, "phylogeny.dating"),
    ({"phylogeny": {"dating": {"lsd2": 1}}}, "phylogeny.dating.lsd2"),
])
def test_sections_require_mappings(config, section):
    with pytest.raises(ValueError, match=f"{section} must be a mapping"):
        validate_keys(config)


@pytest.mark.parametrize("image", [False, 12, [], {}, "", "  "])
def test_invalid_container_images(image):
    with pytest.raises(ValueError, match="container_image must be null or"):
        validate_keys({"container_image": image})


def test_auto_image_follows_version_without_git(tmp_path):
    (tmp_path / "VERSION").write_text("0.2.0\n")
    assert resolve_container_image("auto", tmp_path, enabled=True) == f"docker://{IMAGE_REPOSITORY}:v0.2.0"
    (tmp_path / "VERSION").write_text("1.0.1\n")
    assert resolve_container_image("auto", tmp_path, enabled=True).endswith(":v1.0.1")
    assert not (tmp_path / ".git").exists()


@pytest.mark.parametrize("image", [None, "/data/custom image.sif", "docker://example/image:v1"])
def test_explicit_image_does_not_need_version_file(tmp_path, image):
    assert resolve_container_image(image, tmp_path, enabled=True) == image


def test_native_auto_does_not_need_a_container_or_version_file(tmp_path):
    assert resolve_container_image("auto", tmp_path, enabled=False) is None


@pytest.mark.parametrize("version", ["", "v0.2.0", "0.2", "01.2.3", "1.2.3-rc1", "1.2.3\n2.0.0"])
def test_invalid_release_versions_fail_clearly(tmp_path, version):
    (tmp_path / "VERSION").write_text(version)
    with pytest.raises(ValueError, match="VERSION must contain"):
        read_version(tmp_path)


def test_missing_version_has_archive_instructions(tmp_path):
    with pytest.raises(ValueError, match="complete workflow checkout or release archive"):
        resolve_container_image("auto", tmp_path, enabled=True)
