"""One release version for source archives, image selection, and CI."""
from pathlib import Path
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


def resolve_container_image(image, root, *, enabled):
    """Resolve auto only for container execution; preserve explicit overrides."""
    if image != "auto":
        return image
    if not enabled:
        return None
    return f"docker://{IMAGE_REPOSITORY}:v{read_version(root)}"
