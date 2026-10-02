# Angiosperm leaf

One representative leaf Illumina RNA-seq run per original NCBI taxid.
[Selection settings](selection.yaml) define the search and BUSCO threshold;
[build settings](build.yaml) pin the software and enable `reuse_from: auto`.
See the [shared selection policy](../README.md#selection-policy).

The [effective rules](rules/select_rules.tsv) include fronds and exclude
single-cell, spatial, cell-type/LCM, and special RNA assays.
[Migration notes](rules/initial_migration.json) record their origins.
Historical exclusions are supplied; accepted runs are registered from reviewed
builds.

## Setup

Complete [installation](../../docs/running.md#installation-and-normal-execution)
and run from the repository root:

```bash
conda activate phenoradar_prep
export WORKFLOW_PYTHON="$(command -v python)"
./run_build.sh fetch-software --config datasets/angiosperm_leaf/build.yaml
```

AMALGKIT runs inside the pinned image. For metadata-only work, a smaller host
[environment](../environment.yaml) is available.

## Updating metadata

Edit curation rules and exclusions as needed, then choose a new attempt:

```bash
attempt=work/datasets/angiosperm_leaf/refresh_01
sbatch run_metadata.sh update --work "$attempt"
```

The job retrieves NCBI metadata, curates it, and writes representative metadata.
To preview, add `--dry-run`; publish the reviewed candidate with
`./run_metadata.sh accept --candidate "$attempt/candidate"`.

To select replacements after recording BUSCO failures or changing exclusions:

```bash
sbatch run_metadata.sh update --work work/datasets/angiosperm_leaf/reselect_01 \
  --metadata "$attempt/curate/metadata/metadata.tsv"
```

This reuses the saved candidates without querying NCBI. For curation rule changes
or processing retries, see [restart options](../README.md#choosing-where-to-restart).

## Building the database

After the metadata job finishes, submit a new build. The configuration keeps
`angiosperm_leaf` as its dataset prefix and automatically appends a UTC
preparation timestamp:

```bash
./run_build.sh plan --config datasets/angiosperm_leaf/build.yaml
./run_build.sh submit --config datasets/angiosperm_leaf/build.yaml
```

The command prints its exact build path, for example
`results/angiosperm_leaf_20261002T000000Z/`. Completed products are in
that build's `database/`. Compatible stages from earlier builds are reused
automatically. An optional `--name <exact-name>` sets a name explicitly.

`results/angiosperm_leaf_latest` points to the completed build with the
newest preparation time. Use its `database/` to access current products; use an
exact timestamped path to pin an analysis. The database manifest records build
preparation/completion dates and metadata selection/acceptance/fetch dates;
`database/provenance/metadata_provenance.json` preserves the adopted history.

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
