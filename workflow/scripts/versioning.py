"""One release version for source archives, image selection, and CI."""
from pathlib import Path
import os
import re


IMAGE_REPOSITORY = "ghcr.io/mkrg01/phenoradar_prep"


def parse_version(value):
    """Accept stable major.minor.patch versions, without a leading v."""
    if not isinstance(value, str) or not re.fullmatch(r"(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)", value):
        raise ValueError("VERSION must contain a stable major.minor.patch version, for example 0.2.0")
    return tuple(int(part) for part in value.split("."))


def read_version(root):
    path = Path(root) / "VERSION"
    try:
        value = path.read_text().strip()
    except FileNotFoundError as error:
        raise ValueError(f"Missing {path}; use a complete workflow checkout or release archive") from error
    parse_version(value)
    return value


def resolve_container_image(root):
    """Use the GHCR release recorded in the checkout's VERSION file."""
    frozen = os.environ.get('PHENORADAR_FROZEN_IMAGE')
    if frozen:
        # Set only by the saved-analysis launcher after validating its manifest.
        if not re.fullmatch(r'docker://' + re.escape(IMAGE_REPOSITORY) +
                            r'(?::v[0-9]+\.[0-9]+\.[0-9]+|@sha256:[0-9a-f]{64})', frozen):
            raise ValueError('invalid frozen analysis container image')
        return frozen
    return f"docker://{IMAGE_REPOSITORY}:v{read_version(root)}"
