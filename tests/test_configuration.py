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
    validate_keys({"run_name": "c4_run1"})


@pytest.mark.parametrize("config,path", [
    ({"unknown": {}}, "unknown"),
    ({"inputs": {"metdata": "input/metadata.tsv"}}, "inputs.metdata"),
    ({"odb": {"existing_results": "/data/snapshot"}}, "odb.existing_results"),
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
