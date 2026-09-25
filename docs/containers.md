# Container deployment

[Documentation](index.md) · [Running the workflow](running.md)

## Set up Singularity

Use Linux x86-64 with host Snakemake and `singularity` on compute nodes' `PATH`.
Apptainer must provide its `singularity` compatibility command. Pipeline tools
run inside the container by default.

Use a published release checkout: the GHCR image tag follows [VERSION](../VERSION).
Before the first mapping/analysis job or DAG dry-run, fetch and check it:

```bash
./run_pipeline.sh --prepare-container --cores 1 --resources mem_gb=4
```

Alternatively, submit `sbatch --cpus-per-task=1 --mem=8G run_pipeline.sh --prepare-container`
and wait for completion. This needs no dataset inputs and reuses an existing
image. Repeat after changing releases; unpublished images cannot be pulled.

The repository is mounted automatically, including `input/`, `builds/`,
`analyses/`, `imports/`, and `resources/`. Keep imported products inside it.

## GeneGalleon for incremental datasets

[build.yaml](../config/build.yaml) pins GeneGalleon's version, exact source
revision, and OCI image digest separately. Missing source/SIF files are fetched
when upstream work is needed, or ahead of time with `./run_build.sh fetch-software`.
They are cached under `genegalleon.cache_dir` and checksummed in each build.

`repository` and `image` allow local overrides. An HTTPS SIF URL requires
`image_sha256`; OCI pulls record the resulting SIF hash. Update version, revision,
and image digest together. Incompatible existing products cause reuse conflicts;
see [builds](datasets.md#pinned-genegalleon-source-and-sif).

## Native execution

To use host Conda environments instead of the container:

```bash
./run_pipeline.sh --software-deployment-method conda \
  --cores 16 --resources mem_gb=192 \
  --configfile analyses/analysis001/pipeline.yaml -- all
```

## Build and publish with GitHub Actions

Commit environment and installer changes; CI generates the Dockerfile and saves
it as the `container-recipe` artifact. Container builds use that generated file.
To publish a matching image, update [VERSION](../VERSION) and follow
[releases](releases.md). No Dockerfile commit is needed.

For local image builds only, run `python workflow/scripts/generate_container.py`
first using the pinned [workflow environment](../environment.yaml).
