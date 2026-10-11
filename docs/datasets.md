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
| `plan` | Verify reusable products, show processing counts, and save `reports/plan/plan.tsv` |
| `submit` | Save settings, verify reuse candidates on Slurm, then schedule required work |
| `status --build PATH` | Show recorded sample progress and database publication state |

Monitor submitted work with `status`. The completed database is saved in
`results/{build_name}/database/`. Builds record BUSCO scores without filtering samples.

## Plans and progress

`plan` checks sample identities, processing conditions, and product checksums.
Successful checks are cached for unchanged files. `status` reports progress;
use `status --verify` when you also need to check product integrity.

| Option | Command | Result |
| --- | --- | --- |
| `--outdir DIR` | `plan` | Write `plan.tsv` and `needs_attention.tsv` to `DIR` |
| `--split` | `plan` | Also write `reuse.tsv`, `pending.tsv`, `conflicts.tsv`, and `excluded.tsv` |
| `--outdir DIR` | `status` | Write `status.tsv` and `failed.tsv` to `DIR` |
| `--json` | `status` | Print the detailed report as JSON |

## Run an analysis

Edit [config/analysis.yaml](../config/analysis.yaml) to choose samples, the trait
column, and optional analyses. For expression only, set `phylogeny.trees: []`,
disable `phylogeny.contrast_pairs.enabled`, and set `inputs.species_trait: null`
if traits are unused.
Replace `{analysis_name}` with a name for your analysis.

```bash
./run_analysis.sh plan --build results/{build_name}
./run_analysis.sh submit --build results/{build_name} --name {analysis_name}
./run_analysis.sh status --analysis results/{build_name}/downstream/{analysis_name}
```

Analyses inherit the database's lineage, genetic code, and ODB clade. The default
target `all` runs enabled analyses and collects [PhenoRadar inputs](phenoradar_inputs.md).
Several analyses can share one database. Use an [analysis target](running.md#targets)
to run a single branch.

## Saved settings and previews

The first submission saves inputs and settings. Resume with `--build` or
`--analysis`; source config edits do not change a saved run. Create a new run
for changed scientific settings, code, or software. Use
[resource overrides](running.md#resource-budgets) to change job resources.

Build names include UTC timestamps by default; retry with the path printed by
`submit`. `--name` sets an exact name. With `name_mode: fixed`, repeating `submit`
resumes the configured name. See [build settings](configuration.md#build-settings).

`submit --dry-run` saves a build and previews its preparation job without
submitting. Sample jobs are created during preparation. Resume with `--build`
to use that saved build.

## Failure and recovery behavior

Check the run's `status` and `jobs/logs/`, then retry with its saved path:

```bash
./run_build.sh submit --build results/{build_name}
./run_analysis.sh submit --analysis results/{build_name}/downstream/{analysis_name}
```

Completed samples are skipped; incomplete samples resume from verified
checkpoints. A failed sample does not block other samples, but all must succeed
before publication. Active or unresolved Slurm jobs block duplicate submissions.
Add `--resources retry.yaml` to adjust job resources.

Missing or modified sample products are rebuilt when their reads are available.
Metadata, settings, or shared manifest conflicts stop the build; check
`needs_attention.tsv` or `status` for the reason. A damaged reference pinned by
`reference_id` requires a new build without that pin.

For interrupted jobs, see [storage cleanup](outputs.md#storage-cleanup).
To abandon incomplete samples in a curated dataset, use
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

For a pilot, add `--species-list pilot.txt` with one species or sample ID per line.
Resume with the same `--build` path and without the list to finish all samples.
Use `--until quant` to finish sample processing without publishing a database.

## Reusing completed databases

Set `reuse_from` in `config/build.yaml` to a completed database or a list:

```yaml
reuse_from:
  - results/leaf_set_a/database/
  - results/leaf_set_b/database/
```

Reuse requires matching sample identities and processing conditions. Intact
products can be reused without the original FASTQs.

Keep source databases available until Slurm preparation finishes verifying and
importing their products. The new build then uses its own cache.

Set `reuse_from: auto` to discover previous builds under `results/`, or select an
unfinished build by its directory or `build.json`. Only registered assembly,
BUSCO, and quantification checkpoints are reused; incompatible stages run again.
Translation, mapping, and expression caches are reused from explicitly selected
completed databases. Set `reuse_from: null` to process samples from scratch.

Products may share hard links: **do not edit generated files in place**.

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

Build settings pin the source revision and container. Downloads are cached under
`resources/software/genegalleon/`; `repository` and `image` select local copies.
