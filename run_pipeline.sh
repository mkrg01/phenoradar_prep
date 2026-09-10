#!/usr/bin/env bash
#SBATCH --job-name=phenoradar_prep
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=192G
#SBATCH --time=21-00:00:00
#SBATCH --partition=debug
#SBATCH --signal=B:INT@120
#SBATCH --output=logs/pipeline-%j.log
#SBATCH --open-mode=append
#SBATCH --export=ALL

# Slurm: activate the workflow environment, then submit from the repository root:
# mkdir -p logs
# sbatch --partition=YOUR_PARTITION run_pipeline.sh --configfile config/mydata.yaml
# Direct: ./run_pipeline.sh --software-deployment-method conda --cores 16 --resources mem_gb=192
# All steps run locally. The SBATCH comments only apply when using sbatch.
set -Eeuo pipefail

deployment_args=()
allocation_args=()
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    # Slurm copies this script to its spool directory. Use the job's working
    # directory (submission directory or sbatch --chdir), not the script path.
    root=$PWD
    if [[ ! -f "$root/workflow/Snakefile" ]]; then
        printf '%s\n' 'Submit from the repository root, or set sbatch --chdir to that directory.' >&2
        exit 2
    fi
    if [[ "${SLURM_JOB_NUM_NODES:-${SLURM_NNODES:-1}}" != 1 || "${SLURM_NTASKS:-1}" != 1 ]]; then
        printf '%s\n' 'Use one node and one task; set parallelism with sbatch --cpus-per-task.' >&2
        exit 2
    fi
    batch_cores=${SLURM_CPUS_PER_TASK:-1}
    if [[ ! "$batch_cores" =~ ^[1-9][0-9]*$ ]]; then
        printf '%s\n' 'SLURM_CPUS_PER_TASK must be a positive integer.' >&2
        exit 2
    fi
    if [[ -n "${SLURM_MEM_PER_NODE:-}" ]]; then
        batch_memory_mib=$SLURM_MEM_PER_NODE
    elif [[ "${SLURM_MEM_PER_CPU:-}" =~ ^[1-9][0-9]*$ ]]; then
        batch_memory_mib=$((SLURM_MEM_PER_CPU * batch_cores))
    else
        printf '%s\n' 'Request a finite memory allocation with sbatch --mem or --mem-per-cpu.' >&2
        exit 2
    fi
    if [[ ! "$batch_memory_mib" =~ ^[1-9][0-9]*$ ]]; then
        printf '%s\n' 'The memory allocation must be positive; --mem=0 is not supported.' >&2
        exit 2
    fi
    # Slurm reports MiB; Snakemake mem_mb uses decimal MB. Leave 4 GB for
    # the workflow process and environment management outside rule budgets.
    batch_memory_mb=$((batch_memory_mib * 1048576 / 1000000))
    workflow_memory_mb=$((batch_memory_mb - 4 * 1000))
    if ((workflow_memory_mb <= 0)); then
        printf '%s\n' 'Allocate more than 4 GB to leave memory for processing steps.' >&2
        exit 2
    fi
    printf 'Slurm job: %s\nWorking directory: %s\n' "$SLURM_JOB_ID" "$root"
    printf 'Workflow budget: %s CPUs, %s.%03d GB (4 GB reserved for overhead)\n' \
        "$batch_cores" "$((workflow_memory_mb / 1000))" "$((workflow_memory_mb % 1000))"
    deployment_args=(--software-deployment-method conda)
    allocation_args=(--cores "$batch_cores" --resources "mem_mb=$workflow_memory_mb")
else
    root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
fi
cd "$root"
if ! command -v "${SNAKEMAKE_BIN:-snakemake}" >/dev/null 2>&1; then
    printf '%s\n' 'Snakemake is unavailable. Activate the workflow Conda environment first.' >&2
    exit 2
fi
export XDG_CACHE_HOME="${PHENORADAR_CACHE_DIR:-$root/.cache}"
mkdir -p "$XDG_CACHE_HOME"

# Translate the user-facing GB budget to Snakemake's standard memory resource.
# Only convert values belonging to --resources, leaving config overrides intact.
append_resource() {
    local resource=$1 memory_gb
    if [[ "$resource" == mem_gb=* ]]; then
        memory_gb=${resource#mem_gb=}
        if [[ ! "$memory_gb" =~ ^[1-9][0-9]*$ || ${#memory_gb} -gt 15 ]]; then
            printf '%s\n' 'mem_gb must be a positive integer in GB (for example, mem_gb=192).' >&2
            exit 2
        fi
        resource="mem_mb=$((memory_gb * 1000))"
    fi
    workflow_args+=("$resource")
}

# Keep targets after -- while applying the executor and allocation limits last.
workflow_args=()
targets=()
reading_resources=false
while (($#)); do
    if [[ "$1" == -- ]]; then
        shift
        targets=("$@")
        break
    fi
    case "$1" in
        --resources|--res)
            workflow_args+=(--resources)
            reading_resources=true ;;
        --resources=*|--res=*)
            workflow_args+=(--resources)
            append_resource "${1#*=}"
            reading_resources=true ;;
        -*)
            workflow_args+=("$1")
            reading_resources=false ;;
        *)
            if "$reading_resources"; then
                append_resource "$1"
            else
                workflow_args+=("$1")
            fi ;;
    esac
    shift
done

# exec forwards signals to Snakemake and preserves its exit status in both modes.
exec "${SNAKEMAKE_BIN:-snakemake}" --printshellcmds --rerun-incomplete \
    --snakefile workflow/Snakefile \
    "${deployment_args[@]}" "${workflow_args[@]}" \
    --executor local "${allocation_args[@]}" -- "${targets[@]}"
