# Container deployment

[Documentation](index.md) · [Running the workflow](running.md)

## Set up Singularity

Use Linux x86-64 with `snakemake` and `singularity` on compute nodes' `PATH`.
Apptainer must provide its `singularity` compatibility command.
The launcher uses the container by default; host Conda activation is unnecessary.

Use a published release checkout. The GHCR image matches [VERSION](../VERSION):
`0.2.1` selects `docker://ghcr.io/mkrg01/phenoradar_prep:v0.2.1`.
Before the first analysis or dry-run, fetch the image and check Conda inside it:

```bash
./run_pipeline.sh --prepare-container --cores 1 --resources mem_gb=4
```

For Slurm, use `sbatch --cpus-per-task=1 --mem=8G run_pipeline.sh --prepare-container`
and wait for successful completion before submitting the analysis. This dedicated
step avoids Snakemake 9.8 querying Conda in an image before it has been downloaded.
It needs no dataset inputs, uses the same image cache as analysis jobs, and reuses
an existing image. Repeat after switching to a new release; unpublished versions
cannot be pulled. A cached image still permits offline runs.

The launcher mounts the repository automatically. Keep dataset files in `input/`
and references in `resources/`; no additional bind settings are needed.

## GeneGalleon for incremental datasets

The upstream GeneGalleon dependency is pinned separately in
[config/build.yaml](../config/build.yaml). Build preparation fetches its
source and SIF only when upstream work is missing. Use
`./run_build.sh fetch-software --config config/build.local.yaml` to prepare
these dependencies in advance. See [pinned GeneGalleon source and SIF](datasets.md#pinned-genegalleon-source-and-sif)
for cache reuse, local overrides, and updating the pins.

## Native execution

To create/use workflow environments with host Conda:

```bash
./run_pipeline.sh --software-deployment-method conda \
  --cores 16 --resources mem_gb=192
```

To use tools already on `PATH`, invoke Snakemake directly without deployment flags.

## Build and publish with GitHub Actions

After changing environment definitions or installers, regenerate and commit the
Dockerfile with `python workflow/scripts/generate_container.py`.
[environment.yaml](../environment.yaml) pins Snakemake 9.8.0 for generation.

See [development](development.md) for CI/container checks and
[releases](releases.md) for publication. Changed environment recipes require
a matching image.
