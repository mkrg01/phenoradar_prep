#!/usr/bin/env bash
#SBATCH --job-name=metadata
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=3-00:00:00
#SBATCH --output=slurm-metadata-%j.out

# Run directly or submit from the repository root with sbatch.
set -euo pipefail
root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
if [[ ! -f "$root/workflow/scripts/metadata_catalog.py" ]]; then
    # Slurm executes a copy of the submitted script from its spool directory.
    if [[ -f "$PWD/workflow/scripts/metadata_catalog.py" ]]; then
        root=$PWD
    elif [[ -n "${SLURM_SUBMIT_DIR:-}" && -f "$SLURM_SUBMIT_DIR/workflow/scripts/metadata_catalog.py" ]]; then
        root=$SLURM_SUBMIT_DIR
    else
        echo "Metadata repository not found; submit from its root or use sbatch --chdir <repository>." >&2
        exit 1
    fi
fi
cd "$root"
exec "${WORKFLOW_PYTHON:-${DATASET_PYTHON:-python}}" "$root/workflow/scripts/metadata_catalog.py" "$@"
