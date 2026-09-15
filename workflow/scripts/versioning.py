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


def validate_container_image(image):
    """Accept auto or a release version, never an image URI or local path."""
    if image == "auto":
        return
    try:
        parse_version(image)
    except ValueError:
        raise ValueError(
            "container_image must be auto or a major.minor.patch version without v "
            "(e.g. '0.1.0'); image URIs and SIF paths are no longer supported"
        ) from None


def resolve_container_image(image, root, *, enabled):
    """Resolve the selected release to GHCR only for container execution."""
    validate_container_image(image)
    if not enabled:
        return None
    version = read_version(root) if image == "auto" else image
    return f"docker://{IMAGE_REPOSITORY}:v{version}"
