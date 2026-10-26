"""Fetch the release image before workflows that combine Conda and Apptainer."""
from pathlib import Path
import sys

ROOT = Path(workflow.basedir).parent
sys.path.insert(0, str(ROOT / "workflow/scripts"))
# No container/Conda directive: automatic image pulls happen outside rule
# resources. Pull explicitly so compression observes the allocated budget.
rule prepare_container:
    threads: 1
    resources: mem_mb=16000
    params:
        python=sys.executable, root=str(ROOT), code=str(ROOT / 'workflow/scripts/container_pull.py'),
        cache=lambda wc: str(workflow.persistence.container_img_path)
    shell:
        '{params.python:q} {params.code:q} --root {params.root:q} '
        '--cpus {threads} --memory-mb {resources.mem_mb} --cache-dir {params.cache:q}'
