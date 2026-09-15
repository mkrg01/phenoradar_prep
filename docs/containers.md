# Container deployment

[Documentation](index.md) · [Running the workflow](running.md)

## Set up Singularity

Use Linux x86-64 with `snakemake` and `singularity` on `PATH`, including Slurm
compute nodes. Apptainer must provide its `singularity` compatibility command.
The launcher uses the container by default; host Conda activation is unnecessary.

Use a published release checkout. `container_image: auto` reads [VERSION](../VERSION).
To select a specific GHCR release, set its version without the `v` prefix:

```yaml
container_image: "0.1.0"
```

This selects `docker://ghcr.io/mkrg01/phenoradar_prep:v0.1.0`.
Snakemake downloads and caches the image on first use. Choose a published version
compatible with your checkout.

## External files

External input/reference directories, including symlink targets, need bind mounts:

```bash
./run_pipeline.sh --configfile config/mydata.yaml \
  --cores 4 --resources mem_gb=16 \
  --singularity-args "--cleanenv --bind '$PWD' --bind /data" -- prepare
```

Custom arguments replace the defaults: keep `--cleanenv` and the repository bind.
`--apptainer-args` is an alias. Keep the Release's `image.json` and
[reference snapshots](references.md) with analysis records.

## Native execution

To create and use the workflow environments with host Conda:

```bash
./run_pipeline.sh --software-deployment-method conda \
  --configfile config/mydata.yaml --cores 4 --resources mem_gb=16 -- prepare
```

To use tools already on `PATH`, invoke Snakemake directly without deployment flags.

## Build and publish with GitHub Actions

After changing environment definitions or install scripts, run
`python workflow/scripts/generate_container.py` with Snakemake 9.8.0 and commit
the Dockerfile. [environment.yaml](../environment.yaml) provides a pinned setup.

The [Tests workflow](../.github/workflows/container.yml) runs the Python test suite
and checks the Dockerfile on pull requests and pushes to any branch; manual runs
also build and test the image.
See [releasing a version](releases.md) to publish. Development changes to
environment recipes need a matching image; unpublished versions cannot be pulled.
