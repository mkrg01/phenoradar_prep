#!/usr/bin/env bash
# Compatibility command name: dataset construction is now the build phase.
set -euo pipefail
root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
printf '%s\n' 'run_dataset.sh is now run_build.sh; downstream analyses use run_analysis.sh (see docs/datasets.md).' >&2
exec "$root/run_build.sh" "$@"
