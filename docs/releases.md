# Releasing a version

[Documentation](index.md) · [Container setup](containers.md)

`VERSION` is the only release number to edit. It contains a stable
`major.minor.patch` value without a leading `v`. The runtime and release workflow
derive image URIs and tags from it. CI never changes this file or pushes commits
to `main`.

## Maintainer steps

1. Finish the changes to release. If environment recipes changed, regenerate and
   commit the Dockerfile using [the container guide](containers.md#build-and-publish-with-github-actions).
2. Update `VERSION` to an unused version newer than existing release tags.
3. Commit the changes and push or merge them into `main`.
4. Follow **Actions → Release** and wait for successful completion before using
   the new version's automatic image selection.

For example, changing `VERSION` from `0.1.0` to `0.2.0` publishes the source commit
as `v0.2.0`, the image as `ghcr.io/mkrg01/phenoradar_prep:v0.2.0`, and a GitHub
Release named `v0.2.0`. Users of that release can keep `container_image: auto`.
No separate tag command, version edits elsewhere, or synchronization pull is
needed. Fetch tags locally only if you want to see the new remote tag.

An ordinary commit with the same version does not publish anything. Existing
version tags/images are never moved to another commit. The initial `VERSION`
matches the existing `v0.1.0` tag; adding this automation does not republish it.
Choose a newer number for the first release using this process.

## What Actions does

The [Release workflow](../.github/workflows/release.yml) runs on pushes to `main`
that change `VERSION`. **Run workflow** on `main` can also start or retry it.

1. Validate the version and inspect existing tags/Releases.
2. Run the shared configuration, launcher, release, and Dockerfile checks.
3. Build the exact push commit for Linux x86-64 and test all nine environments and
   real tools without network access.
4. Push a staging image (`build-<commit>`) and the version image to GHCR. If a
   previous attempt already uploaded this version, verify and reuse its exact
   image rather than replacing it with a rebuild.
5. Verify anonymous access to the image, then push an annotated Git tag pointing
   to the tested source commit. Only this tag is pushed; branches are untouched.
6. Create a draft Release with generated notes and an `image.json` asset containing
   the source commit, version, and digest-pinned image URI.
7. Update `latest` for the newest version, then publish the Release.

Release runs are serialized. Finish one release before starting another; GitHub
may replace an older pending run when several are queued. A resumed older
release never moves `latest` backwards. All publication happens in this workflow;
it does not depend on a CI-created tag triggering another workflow.

## Repository setup

The workflow uses `GITHUB_TOKEN` with `contents: write` and `packages: write` in
the publication job. Repository/organization rules must allow the workflow to
create `v*` tags and publish packages. It does not need permission to push to a
protected `main` branch.

GHCR packages must be public so users can run without registry credentials.
On the first publication, GitHub may create a private package. If the anonymous
access check stops the workflow, open the package's **Package settings → Change
visibility → Public**, then re-run the failed workflow. This is a one-time
package setting, not a step for each release. Do not supply a personal token in
dataset configuration files.

For a fork that will publish its own images, update `IMAGE_REPOSITORY` in
`workflow/scripts/versioning.py` to the fork's GHCR repository first. The release
planner checks that this matches the GitHub repository receiving the push.

## Failure recovery

Test/build failures occur before publication. Later failures can leave a staging
or version image, a tag, or a draft Release. No public GitHub Release is created
until all preceding steps succeed.

Use **Re-run failed jobs** on the original run to retry the same commit. A matching
existing image is reused, and a matching tag/draft Release is completed. A fully
published release is skipped. Starting **Run workflow** after `main` has advanced
does not replace a tag belonging to an earlier commit.

If source changes are needed after an image/tag has been published, use a new
version. If failure occurred entirely before publication, fix the source and
use **Run workflow** on `main` to retry the still-unused version. Registry access
errors and version/commit mismatches stop publication instead of being treated
as permission to overwrite existing artifacts.

## Validation

Local release tests use temporary Git repositories and simulated registry/GitHub
responses. They cover tag-only pushes, immutable image reuse, public-access
failures, and retry after tag or asset-upload failures. Actual GitHub permissions,
GHCR transfers, and image execution are checked by Actions during publication.
