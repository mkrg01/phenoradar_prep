#!/usr/bin/env bash
#SBATCH --job-name=phenoradar_prep
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=256G
#SBATCH --time=21-00:00:00
#SBATCH --partition=debug
#SBATCH --output=pipeline-%j.out
#SBATCH --error=pipeline-%j.err

# Adjust SBATCH settings for your cluster; see docs/running.md and docs/containers.md.
# Normally use run_build.sh / run_analysis.sh with their phase configurations.
# Low-level execution requires --configfile <prepared phase>/pipeline.yaml.
# Direct execution uses --cores and --resources mem_gb=...; SBATCH lines are ignored.
set -euo pipefail
launcher_args=("$@")

die() {
    printf '%s\n' "$1" >&2
    exit 2
}

# A PATH-installed Snakemake may use a Python outside PATH. Use that interpreter
# for the common helpers as well, unless the caller explicitly selected one.
workflow_python=${WORKFLOW_PYTHON:-}
if [[ -z "$workflow_python" ]]; then
    snakemake_executable=$(command -v "${SNAKEMAKE_BIN:-snakemake}" || true)
    if [[ -f "$snakemake_executable" ]]; then
        IFS= read -r first_line < "$snakemake_executable" || true
        interpreter=${first_line#\#!}
        if [[ "$first_line" == '#!'* && "$interpreter" != *[[:space:]]* && -x "$interpreter" && "${interpreter##*/}" == python* ]]; then
            workflow_python=$interpreter
        fi
    fi
    workflow_python=${workflow_python:-$(command -v python || command -v python3 || true)}
fi

executor=auto
cluster_mode=false
explicit_local=false
inspection_mode=false
for launcher_option in "$@"; do
    [[ "$launcher_option" == -- ]] && break
    case "$launcher_option" in
        --slurm) executor=slurm; cluster_mode=true ;;
        --age) executor=age; cluster_mode=true ;;
        --cluster) cluster_mode=true ;;
        --cores|-c|--cores=*|--prepare-container) explicit_local=true ;;
        --dry-run|-n|--list-rules|--list-target-rules|--dag|--rulegraph|--filegraph|--summary|--detailed-summary|--lint|--list-changes|--list-untracked) inspection_mode=true ;;
    esac
done
allocation_args=()
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    # Slurm spools the script elsewhere; use the submission directory or --chdir.
    root=$PWD
    [[ -f "$root/workflow/Snakefile" ]] || die 'Submit from the repository root, or set sbatch --chdir to that directory.'
    [[ "${SLURM_JOB_NUM_NODES:-${SLURM_NNODES:-1}}" == 1 && "${SLURM_NTASKS:-1}" == 1 ]] ||
        die 'Use one node and one task; set parallelism with sbatch --cpus-per-task.'
    if ! "$cluster_mode"; then
    executor=local
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
    fi
elif [[ "${JOB_ID:-}" =~ ^[0-9]+$ && -n "${NSLOTS:-}" ]]; then
    root=$PWD
    [[ -f "$root/workflow/Snakefile" ]] || die 'Submit from the repository root with qsub -cwd.'
    if ! "$cluster_mode"; then
        executor=local
        allocation=$("$workflow_python" "$root/workflow/scripts/scheduler.py" allocation)
        read -r batch_cores workflow_memory_mb <<< "$allocation"
        allocation_args=(--cores "$batch_cores" --resources "mem_mb=$workflow_memory_mb")
    fi
else
    root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
fi
cd "$root"

if [[ "$executor" == auto ]]; then
    if "$explicit_local" && ! "$cluster_mode"; then
        executor=local
    elif ! "$cluster_mode" && [[ ! -f "$root/workflow/scripts/scheduler.py" ]]; then
        executor=local
    else
        executor=$("$workflow_python" "$root/workflow/scripts/scheduler.py" detect --optional)
        if "$cluster_mode" && [[ "$executor" == local ]]; then
            die 'No scheduler detected; activate the Slurm or AGE client environment.'
        fi
    fi
fi
if [[ "$executor" != local ]]; then
    export PHENORADAR_SCHEDULER=$executor
fi

if [[ "$executor" == slurm ]]; then
    export PHENORADAR_REAL_SBATCH=${PHENORADAR_REAL_SBATCH:-$(command -v sbatch)}
    export PHENORADAR_SUBMIT_PYTHON=${PHENORADAR_SUBMIT_PYTHON:-$workflow_python}
    export PHENORADAR_SUBMIT_STATE="$root/.snakemake/slurm_submissions"
    export PATH="$root/workflow/slurm:$PATH"
fi

snakemake_bin=${SNAKEMAKE_BIN:-snakemake}
command -v "$snakemake_bin" >/dev/null 2>&1 ||
    die "Snakemake is unavailable: $snakemake_bin. Install Snakemake and add it to PATH, or set SNAKEMAKE_BIN; see README.md requirements."
export XDG_CACHE_HOME="$root/.cache"
mkdir -p "$XDG_CACHE_HOME"

# Convert mem_gb only within --resources. Leave targets after -- in "$@".
workflow_args=()
prepare_container=false
reading_resources=false
pipeline_config=workflow/pipeline_defaults.yaml
profile_supplied=false
expect_config=false
while (($#)); do
    arg=$1
    shift
    case "$arg" in
        --) break ;;
        --slurm|--age|--cluster) reading_resources=false; continue ;;
        --prepare-container) prepare_container=true; reading_resources=false; continue ;;
        --resources|--res)
            arg=--resources
            reading_resources=true ;;
        --resources=*|--res=*)
            workflow_args+=(--resources)
            arg=${arg#*=}
            reading_resources=true ;;
        -*) reading_resources=false ;;
    esac
    if "$expect_config"; then pipeline_config=$arg; expect_config=false; fi
    [[ "$arg" == --configfile ]] && expect_config=true
    [[ "$arg" == --configfile=* ]] && pipeline_config=${arg#*=}
    [[ "$arg" == --profile || "$arg" == --profile=* ]] && profile_supplied=true
    if "$reading_resources" && [[ "$arg" == mem_gb=* ]]; then
        memory_gb=${arg#*=}
        [[ "$memory_gb" =~ ^[1-9][0-9]*$ && ${#memory_gb} -le 15 ]] ||
            die 'mem_gb must be a positive integer in GB (for example, mem_gb=192).'
        arg="mem_mb=$((memory_gb * 1000))"
    fi
    workflow_args+=("$arg")
done

# Bootstrap without Conda environment declarations: Snakemake 9.8 can query
# Conda inside the image before its normal image-pull step.
container_setup_args=()
if "$prepare_container"; then
    [[ "$executor" == local ]] || die "--prepare-container is local; omit --slurm."
    (($# == 0)) || die '--prepare-container does not accept analysis targets.'
    container_setup_args=(--snakefile workflow/container.smk --software-deployment-method apptainer)
    set -- prepare_container
fi

export PHENORADAR_STORAGE_CONFIG=${PHENORADAR_STORAGE_CONFIG:-$root/$pipeline_config}
[[ "$pipeline_config" == /* ]] && export PHENORADAR_STORAGE_CONFIG=$pipeline_config
export PHENORADAR_PROJECT_ROOT=${PHENORADAR_PROJECT_ROOT:-$root}
export PHENORADAR_USAGE_ROOT=${PHENORADAR_USAGE_ROOT:-$root/.snakemake/storage}
if [[ "$executor" == local && "${PHENORADAR_STORAGE_ACTIVE:-}" != 1 && -f "$root/workflow/scripts/storage_runtime.py" ]] && ! "$prepare_container" && ! "$inspection_mode"; then
    exec "$workflow_python" "$root/workflow/scripts/storage_runtime.py" run \
        --pipeline "$PHENORADAR_STORAGE_CONFIG" --root "$root" --usage "$PHENORADAR_USAGE_ROOT" \
        --label pipeline -- bash "$root/run_pipeline.sh" "${launcher_args[@]}"
fi
if [[ "$executor" == age ]]; then
    if ! "$profile_supplied"; then
        profile=$("$workflow_python" "$root/workflow/scripts/scheduler.py" profile --root "$root" --pipeline "$PHENORADAR_STORAGE_CONFIG")
        workflow_args+=(--profile "$profile")
    fi
    executor=cluster-generic
fi

# Mount the repository, including dataset snapshots. In local mode the allocation
# limits win; Slurm mode uses per-rule resources and an independent controller.
printf -v container_root '%q' "$root"
exec "$snakemake_bin" --printshellcmds --rerun-incomplete \
    --snakefile workflow/Snakefile \
    --software-deployment-method conda apptainer \
    --apptainer-args "--cleanenv --bind $container_root" "${workflow_args[@]}" \
    "${container_setup_args[@]}" --executor "$executor" "${allocation_args[@]}" -- "$@"
