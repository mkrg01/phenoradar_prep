"""Fetch the release image before workflows that combine Conda and Apptainer."""
from pathlib import Path
import sys

ROOT = Path(workflow.basedir).parent
sys.path.insert(0, str(ROOT / "workflow/scripts"))
from versioning import resolve_container_image

containerized: resolve_container_image(ROOT)

# No conda directive: image retrieval must precede any Conda environment lookup.
rule prepare_container:
    threads: 1
    resources: mem_mb=1000
    shell: "conda info --json > /dev/null"
