# Build a database and run analyses

[Documentation](index.md) · [Configuration](configuration.md)

A **build** produces reusable sample data: assembly, BUSCO, quantification,
OrthoDB mapping, and OG expression. An **analysis** selects samples from a
completed database. Change traits, sample selection, or the analysis BUSCO
threshold by creating a new analysis from the same database.

## Build a reusable database

Complete [installation](running.md#installation-and-normal-execution), prepare
[metadata](inputs.md), and edit [config/build.yaml](../config/build.yaml).
Set `reuse_from: null` for your first build, or choose a
[completed database](#reusing-completed-databases) to reuse.

Use the build directory name printed by `submit` for `{build_name}` below.

```bash
./run_build.sh plan
./run_build.sh submit
./run_build.sh status --build results/{build_name}
```

| Build command | Purpose |
| --- | --- |
| `plan` | Preview reuse and processing counts; save details to `reports/plan/plan.tsv` |
| `submit` | Save settings, verify reuse candidates on Slurm, then schedule required work |
| `status --build PATH` | Show recorded sample progress and database publication state |

`plan` creates reports without creating a build. Reuse counts are candidates
until verified. `submit` returns after scheduling preparation; monitor the build
with `status`. A successful build publishes `results/{build_name}/database/`.

Each sample runs assembly, longest-CDS extraction, BUSCO, quantification, and
merge in one Slurm allocation. All included samples must finish before publication.
BUSCO scores are recorded without filtering samples. Curated datasets apply
their threshold when [recording results](../datasets/README.md#recording-results).

Use `--until quant` to finish sample processing without publishing a database.
For a pilot, add `--species-list pilot.txt` with one species or sample ID per line.
Resume with the same `--build` path and without the list to finish all samples.

## Check a plan or progress

For builds, `plan` checks identities and settings; `status` reads progress records.
Neither hashes large product files by default. Add `--verify` to check product
integrity immediately:

```bash
./run_build.sh plan --verify
./run_build.sh status --build results/{build_name} --verify
```

`submit` verifies reusable products on Slurm even when `plan` skipped checksums.

| Option | Command | Result |
| --- | --- | --- |
| `--outdir DIR` | `plan` | Write `plan.tsv` and `needs_attention.tsv` to `DIR` |
| `--split` | `plan` | Also write `reuse.tsv`, `pending.tsv`, `conflicts.tsv`, and `excluded.tsv` |
| `--outdir DIR` | `status` | Write `status.tsv` and `failed.tsv` to `DIR` |
| `--json` | `status` | Print the detailed report as JSON |

## Run an analysis

Edit [config/analysis.yaml](../config/analysis.yaml) to choose samples, the trait
column, and optional branches. The supplied config uses `trait: carnivory` and
enables representative trees and contrast pairs. For expression only, set
`phylogeny.trees: []` and disable `phylogeny.contrast_pairs.enabled`.
Set `inputs.species_trait: null` if traits are unused.
Replace `{analysis_name}` with a name for your analysis.

```bash
./run_analysis.sh plan --build results/{build_name}
./run_analysis.sh submit --build results/{build_name} --name {analysis_name}
./run_analysis.sh status --analysis results/{build_name}/downstream/{analysis_name}
```

Analyses inherit the database's lineage, genetic code, and ODB clade and require
intact sample products. Several analyses can share a database; optional alignment,
KO, and tree results are computed separately for each analysis.

The default target `all` runs enabled branches and collects
[PhenoRadar inputs](phenoradar_inputs.md). To run one branch, use an
[analysis target](running.md#targets).

## Saved settings and previews

The first submission saves inputs and settings. Resume with `--build` or
`--analysis` to use those saved conditions. Create a new build or analysis for
scientific changes; use [resource overrides](running.md#resource-budgets) to
change CPU, memory, time, or concurrency.

With `name_mode: timestamp` (the supplied config), `name` is the dataset prefix.
Each build submission without `--name` creates `results/{name}_YYYYMMDDTHHMMSSZ/`
using UTC preparation time, with numbered suffixes for same-second builds.
Use `--build` to retry, including after `--dry-run`.

`--name` sets an exact name. With `name_mode: fixed`, or when `name_mode` is
absent, `name` is the exact build directory name and repeating `submit` resumes
that build.

`results/<name>_latest` points to the completed timestamp-mode build with the
newest preparation time. Use a specific build path to pin an analysis.

Set `metadata_provenance` to the adopted metadata's `provenance.json`; its
checksum must match the input table. The build preserves metadata history and
the origins of reused products.

`./run_build.sh submit --dry-run` saves a build and previews its preparation job
without submitting. Sample jobs are determined after verification, so this is
not a full workflow preview. Keep generated `pipeline.yaml` unchanged. The legacy
build `prepare` command remains available for compatibility.

## Failure and recovery behavior

Check the run's `status` and `jobs/logs/`, then retry with its saved path:

```bash
./run_build.sh submit --build results/{build_name}
./run_analysis.sh submit --analysis results/{build_name}/downstream/{analysis_name}
```

Completed samples are skipped; incomplete samples resume from verified
checkpoints. A failed sample does not block other samples, but all must succeed
before the database is published. Active or unresolved Slurm submissions block
duplicate submissions. Add `--resources retry.yaml` to adjust job resources.

During build preparation, missing or modified sample products are rebuilt when
their reads are available. Metadata, settings, or shared manifest conflicts stop
the build. Check `needs_attention.tsv` for issues found by `plan`, or `status`
for preparation failures. If a damaged reference is pinned by `reference_id`,
create a new build without that pin.

After a code or software update, create a new build and reuse compatible products
from the previous one. Existing builds remain bound to their saved implementation.
For interrupted jobs, see [storage cleanup](outputs.md#storage-cleanup). To abandon
incomplete samples, use
[`record --exclude-failed`](../datasets/README.md#abandoning-incomplete-samples).

## Manually excluding unusable accessions

Edit [config/excluded_accessions.tsv](../config/excluded_accessions.tsv):

```tsv
accession	reason
SRR123456	download_failed_repeatedly
LOCAL_BAD1	unusable_reads
```

`accession` matches the metadata `run` exactly; `reason` is optional. Set
`excluded_accessions` to this file, or `null` to disable exclusions. Exclusions
apply before computation and reuse. After edits, create a new build; review its
excluded samples in `excluded_runs.tsv`.

## Updating samples

Edit metadata and submit a new build with `reuse_from: auto`, or select a
previous build or database. Added samples run missing work; removed samples
are omitted from the new outputs. Changing a run creates a new sample.
Previous builds remain available.

## Reusing completed databases

Set `reuse_from` in `config/build.yaml` to a completed database or a list:

```yaml
reuse_from:
  - results/leaf_set_a/database/
  - results/leaf_set_b/database/
```

Reuse considers samples in the new metadata after exclusions and requires
matching identities, lineage, genetic code, ODB reference, and build conditions.
Intact products can be reused without the original FASTQs; see
[private input guidelines](inputs.md#file-formats).

Product checksums are verified on Slurm. Successful checks are cached for
unchanged files; the first check of a copied database can take longer. Keep
source databases available until preparation finishes importing their products.
The new build then uses its own cache.

Set `reuse_from: auto` to discover previous builds under `results/`, or select an
unfinished build by its directory or `build.json`. Only registered assembly,
BUSCO, and quantification checkpoints are reused; incompatible stages run again.
Translation, mapping, and expression caches are reused from explicitly selected
completed databases. Set `reuse_from: null` to process samples from scratch.

Reuse origins are saved in `reuse.json` and database provenance. Products may
share hard links: **do not edit generated files in place**.

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

Slurm preparation fetches GeneGalleon when sample processing needs it. To download
in advance:

```bash
./run_build.sh fetch-software
```

Build settings pin the source revision and container; downloads are cached under
`resources/software/genegalleon/`. `genegalleon.image_uri`
must use `docker://<registry>/<image>@sha256:<64 lowercase hexadecimal characters>`.
`repository` and `image` are optional local overrides. AMALGKIT rRNA and
contamination filtering are disabled by this workflow.
