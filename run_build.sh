#!/usr/bin/env bash
# Run from any directory. Heavy work is submitted only by the explicit submit command.
set -euo pipefail
root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd "$root"
exec "${WORKFLOW_PYTHON:-${DATASET_PYTHON:-python}}" "$root/workflow/scripts/dataset.py" "$@"
