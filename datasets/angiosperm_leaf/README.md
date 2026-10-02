# Angiosperm leaf

One representative leaf Illumina RNA-seq run per original NCBI taxid.
[Selection settings](selection.yaml) define the search and BUSCO threshold;
[build settings](build.yaml) pin the software and enable `reuse_from: auto`.
See the [shared selection policy](../README.md#selection-policy).

Uses the standard AMALGKIT `plantae` [curation rules](select_rules.tsv).
Historical exclusions are supplied; accepted runs are registered from reviewed
builds.

## Setup

Complete [installation](../../docs/running.md#installation-and-normal-execution)
and run from the repository root:

```bash
conda activate phenoradar_prep
./run_build.sh fetch-software --config datasets/angiosperm_leaf/build.yaml
```

AMALGKIT runs inside the pinned image. For metadata-only work, a smaller host
[environment](../environment.yaml) is available.

## Updating metadata

Edit curation rules and exclusions as needed, then submit an update:

```bash
sbatch run_metadata.sh update
```

The job retrieves NCBI metadata, curates it, and publishes representative
metadata. It prints a UTC timestamped work path under
`work/datasets/angiosperm_leaf/`. To review first, add `--dry-run`, then use the
printed `accept` command after inspection.

To select replacements after recording BUSCO failures or changing exclusions:

```bash
sbatch run_metadata.sh update \
  --metadata work/datasets/angiosperm_leaf/20261002T000000Z/curate/metadata/metadata.tsv
```

This reuses curated metadata without querying NCBI. Substitute the work path
printed by the original job. For curation rule changes or processing retries, see
[restart options](../README.md#choosing-where-to-restart).

## Building the database

After the metadata job finishes, submit a new build:

```bash
./run_build.sh plan --config datasets/angiosperm_leaf/build.yaml
./run_build.sh submit --config datasets/angiosperm_leaf/build.yaml
```

Builds use `results/angiosperm_leaf_<UTC timestamp>/` and reuse compatible stages
automatically. The command prints the exact path; `--name` sets a name manually.
Current products are in `results/angiosperm_leaf_latest/database/`. Use an exact
build path to pin an analysis. See [versioning and provenance](../../docs/datasets.md#saved-settings-and-previews).

## Recording QC and retrying

Inspect QC decisions, then record the reviewed results:

```bash
./run_metadata.sh record --build results/angiosperm_leaf_20261002T000000Z --dry-run
./run_metadata.sh record --build results/angiosperm_leaf_20261002T000000Z
```

`record` updates acceptance and exclusion tables, including from partial builds.
Use `--runs SRR123 SRR456` for a reviewed subset. New runs below the BUSCO
threshold enter the exclusion table; the next metadata update selects alternatives.
For manual exclusions, preserve the table's columns and leave unused cells blank.

For execution failures, retry the same build:

```bash
./run_build.sh submit --build results/angiosperm_leaf_20261002T000000Z
```

Add `--resources <retry.yaml>` for resource changes. Successful stages are reused
and the build keeps its original sample selection. Repeating `submit --config`
without `--build` prepares another build; use the printed path when retrying,
including after `--dry-run`.
