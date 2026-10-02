#!/usr/bin/env bash
# Metadata preparation never submits assembly or analysis jobs.
set -euo pipefail
root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd "$root"
exec "${WORKFLOW_PYTHON:-${DATASET_PYTHON:-python}}" "$root/workflow/scripts/metadata_catalog.py" "$@"
