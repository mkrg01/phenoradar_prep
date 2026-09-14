# Container deployment

[Documentation](index.md)

Singularity is the default execution mode for `run_pipeline.sh`, both directly
and under Slurm. Snakemake runs on the host; processing tools run in a Linux x86-64
image containing all nine workflow Conda environments, including LSD2, MonoPhy,
and ASTRAL-IV int128.

## Set up Singularity

Make the `singularity` command available on the execution host, including Slurm
compute nodes. SingularityCE and Apptainer are supported; the tested Snakemake
9.8.0 calls `singularity`, so an Apptainer installation must provide its
`singularity` compatibility command on `PATH`. Check it with:

```bash
singularity --version
```

Use a published release checkout. The default `container_image: auto` reads
[VERSION](../VERSION) and selects `docker://ghcr.io/mkrg01/phenoradar_prep:v<VERSION>`.
No image setting is required. This also works with release ZIP/tar archives.
Snakemake pulls and caches the image on first use.

To override this choice, set `container_image` to another URI or an absolute SIF
path. Each new [release](releases.md) includes an `image.json` asset with the
exact image digest. Replace `RELEASE_DIGEST` below with that digest to fix the
image contents:

```yaml
container_image: docker://ghcr.io/mkrg01/phenoradar_prep@sha256:RELEASE_DIGEST
```

To prepare a local SIF before submitting a job, download it on a host with
network access:

```bash
singularity pull phenoradar_prep.sif \
  docker://ghcr.io/mkrg01/phenoradar_prep@sha256:RELEASE_DIGEST
```

Then set `container_image` to that SIF's absolute path on storage accessible to
the execution host:

```yaml
container_image: /absolute/path/to/phenoradar_prep.sif
```

## Run the workflow

After making `snakemake` available on `PATH` and configuring your inputs:

```bash
./run_pipeline.sh --configfile config/mydata.yaml \
  --cores 4 --resources mem_gb=16 -- prepare
```

For Slurm, use `sbatch run_pipeline.sh --configfile config/mydata.yaml` after
adjusting the [allocation settings](running.md#slurm). Both modes enable
`--software-deployment-method conda apptainer` automatically. Snakemake calls
this deployment method `apptainer` for both SingularityCE and Apptainer; the
`conda` component activates the environments already built into the image.
The workflow supplies the image's Conda base path so Snakemake does not require
a host Conda installation in this mode.
An explicit `container_image: null` stops container execution with setup
instructions; use `auto` to enable automatic selection.

The launcher passes `--cleanenv` and binds the repository at its existing path.
For external input, reference, and work directories, including symlink targets,
add bind mounts that keep paths identical inside and outside the container:

```bash
./run_pipeline.sh --configfile config/mydata.yaml \
  --cores 4 --resources mem_gb=16 \
  --singularity-args "--cleanenv --bind '$PWD' --bind /data" -- prepare
```

`--singularity-args` and `--apptainer-args` are aliases. Supplying either replaces
the launcher's default arguments, so include `--cleanenv` and the repository bind
alongside your additional mounts. See [running the workflow](running.md) for
targets, pilots, resource budgets, and resuming.

## Native execution

For native Conda execution, install Conda on the host and explicitly select Conda
to create and use the
environments in `workflow/envs/` on the host:

```bash
./run_pipeline.sh --software-deployment-method conda \
  --configfile config/mydata.yaml --cores 4 --resources mem_gb=16 -- prepare
```

This override also works under Slurm. To use tools already installed on `PATH`,
invoke Snakemake directly without deployment flags. Environment definitions pin
the main packages but are not complete dependency lockfiles.

## Build and publish with GitHub Actions

The [Container workflow](../.github/workflows/container.yml) runs as follows:

| Trigger | Action |
| --- | --- |
| Pull request or push to `main` | Check Dockerfile consistency and run related tests |
| Manual **Run workflow** | Build the image and test all nine environments and real tools |

The separate [Release workflow](../.github/workflows/release.yml) publishes when
`VERSION` changes on `main`. It derives the Git tag, image tag, and GitHub Release
from that file, without committing to `main`. See [releasing a version](releases.md)
for the maintainer steps and recovery after a failed publication. The Container
workflow only checks/builds images; it does not publish them.

After changing environment definitions or install scripts, run
`python workflow/scripts/generate_container.py` in the Snakemake 9.8.0 workflow
environment and commit the generated Dockerfile with those changes. The helper
uses Snakemake's `--containerize` feature; CI checks that the result is current.

## Reproducibility

Archive the workflow commit, resolved configuration, image digest (or SIF
checksum), and reference snapshots with an analysis. `run.json` records the
selected image and its environment/source manifest. Retain the built image:
rebuilding can resolve different transitive Conda dependencies. Reference data
is managed separately, and API-backed stages still need network access;
see [references](references.md).

Automatic selection follows the local `VERSION`, not `latest`. Updating an image
does not change an existing checkout. For development changes to environment
recipes, build a matching image or finish a new release before using `auto`;
unpublished versions cannot be downloaded yet.
