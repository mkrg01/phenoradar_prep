# Container deployment

[Documentation](index.md)

GitHub Actions builds and publishes a Linux x86-64 image containing all nine
workflow Conda environments, including LSD2, MonoPhy and ASTRAL-IV int128.
The image uses the same verified upstream installers as native Conda deployment;
no custom Conda packages are required.

## Build and publish with GitHub Actions

The [Container workflow](../.github/workflows/container.yml) runs as follows:

| Trigger | Action |
| --- | --- |
| Pull request or push to `main` | Check Dockerfile consistency and run related tests |
| Manual **Run workflow** | Build the image and test all nine environments and real tools |
| Push a `v*` tag | Run the checks, build and test the image, then publish |

Releases publish to `ghcr.io/mkrg01/phenoradar_prep` with the version tag and
`latest`, using `GITHUB_TOKEN`. Manual runs do not publish images.

After changing environment definitions or install scripts, run
`python workflow/scripts/generate_container.py` in the Snakemake 9.8.0 workflow
environment and commit the generated Dockerfile with those changes. The helper
uses Snakemake's `--containerize` feature; CI checks that the result is current.

## Run with Apptainer

After a release is published, set `container_image` in your dataset configuration.
Replace `RELEASE_DIGEST` below with its actual image digest to fix the version:

```yaml
container_image: docker://ghcr.io/mkrg01/phenoradar_prep@sha256:RELEASE_DIGEST
```

Snakemake runs on the host using [environment.yaml](../environment.yaml), with
Apptainer available on the compute node. Run from the matching workflow checkout:

```bash
./run_pipeline.sh --software-deployment-method conda apptainer \
  --configfile config/mydata.yaml --cores 4 --resources mem_gb=16 \
  --apptainer-args "--cleanenv --bind $PWD" -- prepare
```

The same flags work with `sbatch run_pipeline.sh`. Add external input, reference
and work directories, including symlink targets, to `--bind`, keeping their
paths identical inside and outside the container. An existing SIF file can also
be selected by its absolute path in `container_image`.

For native Conda execution, leave `container_image: null` and use only
`--software-deployment-method conda`.

## Reproducibility

Archive the workflow commit, resolved configuration, image digest (or SIF
checksum), and reference snapshots with an analysis. `run.json` records the
selected image and its environment/source manifest. Retain the built image:
rebuilding can resolve different transitive Conda dependencies. Reference data
is managed separately, and API-backed stages still need network access;
see [references](references.md).
