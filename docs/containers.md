# Container deployment

[Documentation](index.md) · [Running the workflow](running.md)

## Set up Singularity

Use Linux x86-64 with `snakemake` and `singularity` on `PATH`, including Slurm
compute nodes. Apptainer must provide its `singularity` compatibility command.
The launcher uses the container by default; host Conda activation is unnecessary.

Use a published release checkout. `container_image: auto` selects
`docker://ghcr.io/mkrg01/phenoradar_prep:v<VERSION>` from [VERSION](../VERSION).
Snakemake downloads and caches the image on first use. An existing checkout
keeps its version until you update it.

## Custom images and external files

Override `container_image` with a release URI, a digest URI from the Release's
`image.json`, or an absolute SIF path. To prepare a SIF on a networked host:

```bash
singularity pull phenoradar_prep.sif \
  docker://ghcr.io/mkrg01/phenoradar_prep@sha256:RELEASE_DIGEST
```

Replace `RELEASE_DIGEST`, then configure a path accessible from compute nodes:

```yaml
container_image: /absolute/path/to/phenoradar_prep.sif
```

External input/reference directories, including symlink targets, need bind mounts:

```bash
./run_pipeline.sh --configfile config/mydata.yaml \
  --cores 4 --resources mem_gb=16 \
  --singularity-args "--cleanenv --bind '$PWD' --bind /data" -- prepare
```

Custom arguments replace the defaults: keep `--cleanenv` and the repository bind.
`--apptainer-args` is an alias. Keep image digests/SIF checksums and
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

The [Container workflow](../.github/workflows/container.yml) checks pull requests
and pushes to `main`; manual runs also build and test the image.
See [releasing a version](releases.md) to publish. Development changes to
environment recipes need a matching image; unpublished versions cannot be pulled.
