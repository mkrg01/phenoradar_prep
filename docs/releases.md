# Releasing a version

[Documentation](index.md) · [Container setup](containers.md)

## Maintainer steps

1. If environment definitions or install scripts changed, regenerate and commit
   the [Dockerfile](containers.md#build-and-publish-with-github-actions).
2. Update [VERSION](../VERSION) to an unused, higher `major.minor.patch` value
   without a leading `v`.
3. Commit and push or merge into `main`, then wait for **Actions → Release** to
   finish before using the new image. Finish one release before starting another.

The [Release workflow](../.github/workflows/release.yml) tests and publishes the
matching container, Git tag, and GitHub Release. For `VERSION` equal to `0.2.0`,
the tag is `v0.2.0` and the image is `ghcr.io/mkrg01/phenoradar_prep:v0.2.0`.
The Release's `image.json` records its exact digest and source commit.

`VERSION` is the only number to edit. CI never commits to `main`, so no
synchronization pull is needed. Ordinary commits leave published versions intact.

## Repository setup

Allow `GITHUB_TOKEN` to create tags and Releases (`contents: write`) and publish
packages (`packages: write`). The GHCR package must be public. If the first run
stops at the public-access check, set **Package settings → Change visibility →
Public**, then re-run it.

For a fork, update `IMAGE_REPOSITORY` in
[versioning.py](../workflow/scripts/versioning.py) to its GHCR repository.

## Failure recovery

Use **Re-run failed jobs** on the original run to retry the same commit. Existing
matching images and tags are reused; completed Releases are skipped.

If source changes are needed after an image or tag was published, choose a new
version. If nothing was published, fix the source and use **Actions → Release →
Run workflow** on `main` to retry the unused version.
