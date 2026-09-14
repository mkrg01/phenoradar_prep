"""Configuration shape validation, including partial command-line overrides."""
from pathlib import Path

import pytest
import yaml

from configuration import KEYS, validate_keys

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
