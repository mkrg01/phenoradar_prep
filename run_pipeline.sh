#!/usr/bin/env bash
#SBATCH --job-name=phenoradar_prep
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=192G
#SBATCH --time=21-00:00:00
#SBATCH --partition=debug
#SBATCH --output=pipeline-%j.out
#SBATCH --error=pipeline-%j.err

# Adjust SBATCH settings for your cluster; see docs/running.md and docs/containers.md.
# Make snakemake and singularity available on PATH and edit config/config.yaml:
# sbatch run_pipeline.sh
# Direct execution uses --cores and --resources mem_gb=...; SBATCH lines are ignored.
set -euo pipefail

die() {
    printf '%s\n' "$1" >&2
    exit 2
}

allocation_args=()
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    # Slurm spools the script elsewhere; use the submission directory or --chdir.
    root=$PWD
    [[ -f "$root/workflow/Snakefile" ]] || die 'Submit from the repository root, or set sbatch --chdir to that directory.'
    [[ "${SLURM_JOB_NUM_NODES:-${SLURM_NNODES:-1}}" == 1 && "${SLURM_NTASKS:-1}" == 1 ]] ||
        die 'Use one node and one task; set parallelism with sbatch --cpus-per-task.'
    batch_cores=${SLURM_CPUS_PER_TASK:-1}
    [[ "$batch_cores" =~ ^[1-9][0-9]*$ ]] || die 'SLURM_CPUS_PER_TASK must be a positive integer.'
    if [[ -n "${SLURM_MEM_PER_NODE:-}" ]]; then
        batch_memory_mib=$SLURM_MEM_PER_NODE
    elif [[ "${SLURM_MEM_PER_CPU:-}" =~ ^[1-9][0-9]*$ ]]; then
        batch_memory_mib=$((SLURM_MEM_PER_CPU * batch_cores))
    else
        die 'Request a finite memory allocation with sbatch --mem or --mem-per-cpu.'
    fi
    [[ "$batch_memory_mib" =~ ^[1-9][0-9]*$ ]] || die 'The memory allocation must be positive; --mem=0 is not supported.'
    # Convert Slurm MiB to Snakemake MB, reserving 4 GB for workflow overhead.
    workflow_memory_mb=$((batch_memory_mib * 1048576 / 1000000 - 4000))
    ((workflow_memory_mb > 0)) || die 'Allocate more than 4 GB to leave memory for processing steps.'
    printf 'Slurm job: %s\nWorking directory: %s\n' "$SLURM_JOB_ID" "$root"
    printf 'Workflow budget: %s CPUs, %s.%03d GB (4 GB reserved for overhead)\n' \
        "$batch_cores" "$((workflow_memory_mb / 1000))" "$((workflow_memory_mb % 1000))"
    allocation_args=(--cores "$batch_cores" --resources "mem_mb=$workflow_memory_mb")
else
    root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
fi
cd "$root"

snakemake_bin=${SNAKEMAKE_BIN:-snakemake}
command -v "$snakemake_bin" >/dev/null 2>&1 ||
    die "Snakemake is unavailable: $snakemake_bin. Install Snakemake and add it to PATH, or set SNAKEMAKE_BIN; see README.md requirements."
export XDG_CACHE_HOME="$root/.cache"
mkdir -p "$XDG_CACHE_HOME"

# Convert mem_gb only within --resources. Leave targets after -- in "$@".
workflow_args=()
reading_resources=false
while (($#)); do
    arg=$1
    shift
    case "$arg" in
        --) break ;;
        --resources|--res)
            arg=--resources
            reading_resources=true ;;
        --resources=*|--res=*)
            workflow_args+=(--resources)
            arg=${arg#*=}
            reading_resources=true ;;
        -*) reading_resources=false ;;
    esac
    if "$reading_resources" && [[ "$arg" == mem_gb=* ]]; then
        memory_gb=${arg#*=}
        [[ "$memory_gb" =~ ^[1-9][0-9]*$ && ${#memory_gb} -le 15 ]] ||
            die 'mem_gb must be a positive integer in GB (for example, mem_gb=192).'
        arg="mem_mb=$((memory_gb * 1000))"
    fi
    workflow_args+=("$arg")
done

# Quote the bind path for Snakemake's shell command. User deployment flags can
# override these defaults; the local executor and allocation limits apply last.
printf -v container_root '%q' "$root"
exec "$snakemake_bin" --printshellcmds --rerun-incomplete \
    --snakefile workflow/Snakefile \
    --software-deployment-method conda apptainer \
    --apptainer-args "--cleanenv --bind $container_root" "${workflow_args[@]}" \
    --executor local "${allocation_args[@]}" -- "$@"
