# Container deployment

[Documentation](index.md) · [Running the workflow](running.md)

## Set up Singularity

Use Linux x86-64 with `snakemake` and `singularity` on compute nodes' `PATH`.
Apptainer must provide its `singularity` compatibility command.
The launcher uses the container by default; host Conda activation is unnecessary.

Use a published release checkout. The GHCR image matches [VERSION](../VERSION):
`0.2.1` selects `docker://ghcr.io/mkrg01/phenoradar_prep:v0.2.1`.
Snakemake downloads/caches it on first use. Switch checkouts to use another release;
unpublished versions cannot be pulled.

## External files

External inputs/references, including symlink targets, need bind mounts:

```bash
./run_pipeline.sh --cores 16 --resources mem_gb=192 \
  --singularity-args "--cleanenv --bind '$PWD' --bind /data"
```

Custom arguments replace defaults: keep `--cleanenv` and the repository bind.
`--apptainer-args` is an alias.

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
