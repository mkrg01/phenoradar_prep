# Releasing a version

[Documentation](index.md) · [Development](development.md)

## Publish a release

1. Set [VERSION](../VERSION) to an unused, higher `major.minor.patch` value without
   a leading `v`.
2. Commit and push or merge into `main`.
3. Wait for **Actions → Release** to finish before using the image or starting
   another release.

The [release workflow](../.github/workflows/release.yml) runs tests and container
smoke checks, then publishes the image, Git tag, and GitHub Release. CI generates
the Dockerfile from the environment definitions; no manual regeneration is needed.

Tags use `v<VERSION>` and images use `ghcr.io/mkrg01/phenoradar_prep:v<VERSION>`.
The release's `image.json` records the exact image digest and source commit.
`VERSION` is the only version number to edit.

## Repository setup

Allow `GITHUB_TOKEN` to create tags/releases (`contents: write`) and publish
packages (`packages: write`). The GHCR package must be public. If the first run
fails that check, change the package visibility to public and rerun it.

For a fork, update `IMAGE_REPOSITORY` in
[versioning.py](../workflow/scripts/versioning.py).

## Failure recovery

Use **Re-run failed jobs** on the original run to retry the same commit.
Matching images and tags are reused; completed releases are skipped.

If source changes are needed after an image or tag was published, choose a new
version. If nothing was published, fix the source and use **Actions → Release →
Run workflow** on `main` to retry the unused version.
