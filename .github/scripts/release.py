#!/usr/bin/env python3
"""Publish one tested VERSION without creating or updating branch commits."""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "workflow/scripts"))
from versioning import IMAGE_REPOSITORY, parse_version, read_version


def command(args, root=ROOT, *, check=True):
    result = subprocess.run(args, cwd=root, text=True, capture_output=True)
    if check and result.returncode:
        raise RuntimeError(f"{args[0]} failed: {result.stderr.strip() or result.stdout.strip()}")
    return result


def git(*args, root=ROOT):
    return command(["git", *args], root).stdout.strip()


def release_info(repository, tag, root=ROOT):
    result = command(["gh", "api", f"repos/{repository}/releases/tags/{tag}"], root, check=False)
    if result.returncode:
        if "HTTP 404" in result.stderr:
            return None
        raise RuntimeError(f"Cannot check the GitHub release: {result.stderr.strip()}")
    return json.loads(result.stdout)


def stable_tags(root=ROOT):
    versions = {}
    for tag in git("tag", "--list", "v*", root=root).splitlines():
        try:
            versions[tag] = parse_version(tag[1:])
        except ValueError:
            continue
    return versions


def plan(repository, root=ROOT):
    if not re.fullmatch(r"[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+", repository):
        raise ValueError("Expected a GitHub owner/repository")
    if f"ghcr.io/{repository.lower()}" != IMAGE_REPOSITORY:
        raise ValueError("For a fork, update IMAGE_REPOSITORY in workflow/scripts/versioning.py first")
    version = read_version(root)
    if git("status", "--porcelain", root=root):
        raise ValueError("Commit all release changes, including VERSION, before publishing")
    tag = f"v{version}"
    revision = git("rev-parse", "HEAD", root=root)
    tags = stable_tags(root)
    result = {"version": version, "tag": tag, "revision": revision,
              "image_repository": IMAGE_REPOSITORY, "release": "true"}
    if tag in tags:
        tagged_revision = git("rev-list", "-n", "1", tag, root=root)
        if tagged_revision != revision:
            if parse_version(version) < max(tags.values()):
                raise ValueError("VERSION is older than an existing release; choose a newer version")
            result.update(release="false", reason="This version is already tagged; ordinary development does not republish it")
        else:
            info = release_info(repository, tag, root)
            if info and not info["draft"]:
                result.update(release="false", reason="This release is already published")
    elif tags and parse_version(version) <= max(tags.values()):
        raise ValueError("VERSION must be newer than existing release tags")
    return result


def image_digest(reference, root=ROOT):
    result = command(["docker", "manifest", "inspect", "--verbose", reference], root, check=False)
    if result.returncode:
        if any(message in result.stderr.lower() for message in ("manifest unknown", "no such manifest")):
            return None
        raise RuntimeError(f"Cannot inspect {reference}: {result.stderr.strip()}")
    manifest = json.loads(result.stdout)
    digest = manifest.get("Descriptor", {}).get("digest") if isinstance(manifest, dict) else None
    if not isinstance(digest, str) or not re.fullmatch(r"sha256:[0-9a-f]{64}", digest):
        raise ValueError("Expected a single-platform image manifest with a SHA-256 digest")
    return digest


def check_image_labels(reference, version, revision, root=ROOT):
    labels = json.loads(command([
        "docker", "image", "inspect", "--format", "{{json .Config.Labels}}", reference
    ], root).stdout) or {}
    if (labels.get("org.opencontainers.image.version") != version or
            labels.get("org.opencontainers.image.revision") != revision):
        raise ValueError("Image version/source commit does not match; refusing to overwrite a release image")


def publish(repository, local_image, root=ROOT):
    git("fetch", "origin", "--tags", root=root)
    current = plan(repository, root)
    if current["release"] == "false":
        print(current["reason"])
        return
    version, tag, revision = (current[key] for key in ("version", "tag", "revision"))
    reference = f"{IMAGE_REPOSITORY}:{tag}"
    check_image_labels(local_image, version, revision, root)
    # Staging also creates the GHCR package on its first publication, so a missing
    # version manifest can be distinguished from registry authorization errors.
    staging = f"{IMAGE_REPOSITORY}:build-{revision}"
    command(["docker", "tag", local_image, staging], root)
    command(["docker", "push", staging], root)
    digest = image_digest(reference, root)
    if digest:
        # A prior attempt may have pushed the image before a later step failed.
        # Reuse that exact image even if a rebuild resolves newer dependencies.
        immutable = f"{IMAGE_REPOSITORY}@{digest}"
        command(["docker", "pull", immutable], root)
        check_image_labels(immutable, version, revision, root)
        command(["docker", "tag", immutable, local_image], root)
        command(["docker", "run", "--rm", "--network", "none", "-v", f"{root}:/workspace:ro",
                 "-w", "/workspace", local_image, "python", "tests/container_smoke.py"], root)
    else:
        command(["docker", "tag", local_image, reference], root)
        command(["docker", "push", reference], root)
        digest = image_digest(reference, root)
        if digest is None:
            raise RuntimeError("Published image could not be found in the registry")
        immutable = f"{IMAGE_REPOSITORY}@{digest}"

    with tempfile.TemporaryDirectory(prefix="phenoradar-docker-auth-") as anonymous:
        result = command(["docker", "--config", anonymous, "manifest", "inspect", immutable], root, check=False)
    if result.returncode:
        raise RuntimeError("Image is not publicly readable. Set the GHCR package visibility to public, "
                           "then re-run this workflow. Git tag/Release publication was not attempted.")

    tags = stable_tags(root)
    if tag not in tags:
        git("-c", "user.name=github-actions[bot]", "-c",
            "user.email=41898282+github-actions[bot]@users.noreply.github.com",
            "tag", "-a", tag, revision, "-m", f"Release {tag}\n\nImage: {immutable}", root=root)
    else:
        annotation = git("for-each-ref", "--format=%(contents)", f"refs/tags/{tag}", root=root)
        if f"Image: {immutable}" not in annotation:
            raise ValueError("Existing release tag does not record this image digest; refusing to change it")
    # Only the tag is pushed. Never push HEAD or update main from CI.
    git("push", "origin", f"refs/tags/{tag}", root=root)

    info = release_info(repository, tag, root)
    if info and not info["draft"]:
        raise ValueError("Release became public during this run; refusing to modify it")
    if info is None:
        command(["gh", "release", "create", tag, "--repo", repository, "--verify-tag",
                 "--draft", "--generate-notes", "--title", tag], root)
    with tempfile.TemporaryDirectory(prefix="phenoradar-release-") as temporary:
        manifest = Path(temporary) / "image.json"
        manifest.write_text(json.dumps({"version": version, "revision": revision,
                                       "image": f"docker://{immutable}"}, indent=2) + "\n")
        command(["gh", "release", "upload", tag, str(manifest), "--repo", repository, "--clobber"], root)

    git("fetch", "origin", "--tags", root=root)
    latest = parse_version(version) == max(stable_tags(root).values())
    if latest:
        command(["docker", "tag", local_image, f"{IMAGE_REPOSITORY}:latest"], root)
        command(["docker", "push", f"{IMAGE_REPOSITORY}:latest"], root)
    command(["gh", "release", "edit", tag, "--repo", repository, "--draft=false",
             f"--latest={str(latest).lower()}"], root)
    print(f"Published {tag}: docker://{immutable}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=["plan", "publish"])
    parser.add_argument("--repository", default=os.environ.get("GITHUB_REPOSITORY"))
    parser.add_argument("--output", type=Path)
    parser.add_argument("--image", default="phenoradar_prep:ci")
    args = parser.parse_args()
    if not args.repository:
        parser.error("--repository or GITHUB_REPOSITORY is required")
    try:
        if args.action == "publish":
            publish(args.repository, args.image)
        else:
            result = plan(args.repository)
            print(json.dumps(result, indent=2))
            if args.output:
                with args.output.open("a") as handle:
                    for key in ("version", "tag", "revision", "image_repository", "release"):
                        handle.write(f"{key}={result[key]}\n")
    except (ValueError, RuntimeError) as error:
        parser.exit(1, f"{error}\n")


if __name__ == "__main__":
    main()
