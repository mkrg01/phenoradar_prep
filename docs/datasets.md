# Build a database and run analyses

[Documentation](index.md) · [Configuration](configuration.md)

A **build** produces reusable sample data: assembly, BUSCO, quantification,
OrthoDB mapping, and OG expression. An **analysis** selects samples from a
completed build and combines their results. Changing a trait, BUSCO threshold,
or species selection needs a new analysis, while keeping the same database.

## Build a reusable database

Complete [installation](running.md#installation-and-normal-execution), prepare
[metadata](inputs.md), and edit [config/build.yaml](../config/build.yaml).
Set `reuse_from: null` for your first build, or choose a
[completed database](#reusing-completed-databases) to reuse.

```bash
./run_build.sh plan
./run_build.sh submit --name leaf
./run_build.sh status --build results/leaf
```

Review the plan for `reuse`, `pending`, or `conflict` before submitting. Conflicts
must be resolved; ODB reuse is confirmed after translation provides protein hashes.
`submit` returns after scheduling Slurm jobs. A successful build publishes
`results/leaf/database/` and a completion record.

Every included sample must complete all stages. Build execution records BUSCO
scores without filtering samples. For curated datasets, the separate `record`
command applies the dataset's BUSCO threshold and adds below-threshold new runs
to exclusions for future metadata updates. Analysis thresholds remain separate.

To inspect an intermediate stage:

```bash
./run_build.sh submit --name pilot --until busco
# After inspection, finish the same build:
./run_build.sh submit --build results/pilot --until database
```

Endpoints are `assembly`, `busco`, `quant`, and `database` (default), each including
missing prerequisites. For a subset pilot, add `--species-list pilot.txt` with
one biological species ID or exact sample ID per line. Submit again without the
list to finish all samples; an incomplete pilot cannot publish a database.

## Run an analysis

Edit [config/analysis.yaml](../config/analysis.yaml): choose sample selection,
`trait`, and optional branches. The supplied config uses `trait: carnivory` and
enables representative trees and contrast pairs. For expression only, set
`phylogeny.trees: []` and disable `phylogeny.contrast_pairs.enabled`.
Set `inputs.species_trait: null` if traits are unused.

```bash
./run_analysis.sh plan --build results/leaf
./run_analysis.sh submit --build results/leaf --name carnivory
./run_analysis.sh status --analysis results/leaf/downstream/carnivory
```

`trait` selects the column in your trait table; `--name` names the result directory.
Analysis inherits the build's lineage, genetic code, and ODB mapping clade.
It reuses verified sample products and reports missing or modified products as
errors. Several analyses can share a database; optional alignment, KO, and tree
results are computed separately for each analysis.

The default target `all` runs enabled branches and collects
[PhenoRadar inputs](phenoradar_inputs.md). To run one branch, use an
[analysis target](running.md#targets).

## Saved settings and previews

The first `submit` freezes inputs and settings for that name. Repeating it uses
the saved conditions, even if source files have changed. Use a new name for
scientific changes; use [resource overrides](running.md#resource-budgets) for
CPU, memory, time, or concurrency changes.

With `name_mode: timestamp`, `name` becomes the dataset prefix. Each `prepare`
or new `submit` creates `results/<name>_YYYYMMDDTHHMMSSZ/` using UTC preparation
time; same-second builds receive `_02`, `_03`, and so on. `--name` sets an exact
name. Resume with `submit --build results/<printed-build-name>`, including after
`--dry-run`; omitting `--build` creates a new snapshot.

`results/<name>_latest` points to the completed timestamp-mode build with the
newest preparation time. Use a specific build path to pin an analysis.

Set `metadata_provenance` to the adopted metadata's `provenance.json`. The
recorded metadata checksum must match the input table, including with
`--metadata`. Build and database manifests preserve preparation and metadata
history; database manifests also record completion. Fetch dates require
acquisition evidence, and reused products retain their original history.

| Command | Effect |
| --- | --- |
| `plan` | Preview current inputs and settings without saving a run |
| `prepare` | Save inputs and settings without submitting |
| `submit --dry-run` | Save inputs/settings and preview job scripts without submitting |
| `submit` | Save a new run or resume an existing one, then submit jobs |

`--dry-run` fixes the conditions for that name; it does not check the full
Snakemake DAG. Generated `pipeline.yaml` files should not be edited.

## Failure and recovery behavior

Inspect `status`, `squeue`, and the run's `jobs/logs/`, then resubmit:

```bash
./run_build.sh submit --build results/leaf
./run_analysis.sh submit --analysis results/leaf/downstream/carnivory
```

Completed work is reused and failed or missing stages are retried. Active or
unresolved submissions block retries. Failed array tasks block dependent jobs;
let remaining jobs finish or cancel them before resubmitting. A timed-out
controller can leave workers running, so inspect those too.

A retry may restart a failed stage: rnaSPAdes assembly starts from scratch, and
partial read downloads are not guaranteed reusable. Successful samples are
unaffected. For repeated timeouts or memory errors, adjust
[resources](running.md#resource-budgets). See [storage cleanup](outputs.md#storage-cleanup)
for removing intermediates from completed jobs.

## Manually excluding unusable accessions

Edit [config/excluded_accessions.tsv](../config/excluded_accessions.tsv):

```tsv
accession	reason
SRR123456	download_failed_repeatedly
LOCAL_BAD1	unusable_reads
```

`accession` matches the metadata `run` exactly; `reason` is optional. Select this
file with `excluded_accessions` in build settings, or use `null` to disable it.
Exclusions apply before computation or reuse. Prepare a new build after edits;
previous builds remain unchanged. The build saves the original metadata,
exclusions, effective metadata, and `excluded_runs.tsv` for review.

## Updating samples

Edit metadata and submit a new build name with `reuse_from: auto`, or point it to
a previous build or completed database. Added samples run missing work; removed samples leave the
new outputs. Changing a run creates a new sample with its own assembly and
expression. Historical builds remain available.

## Reusing completed databases

Set `reuse_from` in `config/build.yaml` to a completed database or a list:

```yaml
reuse_from:
  - results/leaf_set_a/database/
  - results/leaf_set_b/database/
```

Only samples in the new metadata, after exclusions, are considered. The pipeline
checks identities, checksums, lineage, genetic code, ODB reference, and build
conditions. Conflicting products stop preparation regardless of list order.
Use compatible sources, or `null` for a fresh build.

Build commands cache successful checksum verification in
`.cache/verification.sqlite`. Cached checks require the same path, expected
SHA256, device, inode, size, mtime, and ctime. Changed files are checked again;
missing or corrupt files still fail validation. This cache is optional: it can
be removed, and unavailable caches fall back to full verification. Database
manifests and scientific reuse conditions are unchanged. The first verification
after transferring products still reads their contents.

Verified products are staged into the new build, which then no longer depends
on the source paths. Missing samples run normally. Large products may share
hard links: **do not edit generated files in place**. Reuse origins are recorded
in `reuse.json` and the completed database's provenance.

An unfinished build can resume under its own name or supply its verified
assembly, BUSCO, and quantification stages to a new build. Set `reuse_from` to
its build directory (or `build.json`), or use `reuse_from: auto` to discover
previous builds under `results/`. Automatic reuse checks each sample and stage;
incompatible settings are cache misses and conflicting or damaged matching
products stop preparation. It never imports unregistered partial outputs.

Automatic reuse copies sample-stage receipts; mapping, translation, and expression
caches are reused through explicitly selected completed databases. `null`
disables all reuse. Discovery is frozen during preparation, and the new build
does not depend on source product paths afterward.

For curated representative datasets, use the
[metadata update and adoption commands](../datasets/README.md). An excluded run
can be replaced by another eligible run with the same original taxid, including
the same BioProject. Adoption records can be made before a database completes.
QC recording automatically adds new runs with completed BUSCO results below the
dataset threshold to its accession exclusions. Update metadata and prepare a new
build to process replacements; failed or unfinished BUSCO stages stay eligible
for retry.

## Copying a completed build to another project

Copy the entire **`database/` directory** into the destination project, for example
`results/baseline/database/`, then run:

```bash
./run_analysis.sh submit --build results/baseline --name carnivory
```

Both the build directory and its `database/` path are accepted. Original reads
and working caches are unnecessary; the destination supplies analysis settings,
traits, software, and reference resources.

## Pinned GeneGalleon source and SIF

GeneGalleon software is fetched automatically when assembly, BUSCO, or
quantification needs it. To download in advance:

```bash
./run_build.sh fetch-software
```

Build settings pin the source revision and container. `genegalleon.image_uri`
must use `docker://<registry>/<image>@sha256:<64 lowercase hexadecimal characters>`.
Verified downloads are cached under `resources/software/genegalleon/`;
`repository` and `image` are optional local overrides. AMALGKIT rRNA and
contamination filtering are disabled by this workflow.
