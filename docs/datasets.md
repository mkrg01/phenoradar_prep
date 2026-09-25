# Reusable builds and independent analyses

[Documentation](index.md) · [Configuration](configuration.md)

`run_build.sh` prepares reusable species products through ODB mapping.
`run_analysis.sh` selects species from a completed build and produces expression,
alignment, phylogeny, and PhenoRadar outputs.

## Prepare a manually curated dataset

Maintain `input/metadata.tsv` with **all desired species**, one run per species,
from NCBI or local FASTQs. See [input formats](inputs.md) for required columns and
local-read paths. Keep known unusable runs in the
[manual exclusion list](#manually-excluding-unusable-accessions).

Edit `config/build.yaml` and `config/analysis.yaml` directly. These are the default
configs; `prepare` freezes settings and inputs for each named run. Install the
[workflow environment](../environment.yaml) and prepare the
[workflow image](running.md#installation-and-normal-execution) before submitting mapping or analysis jobs.

## Pinned GeneGalleon source and SIF

Build settings pin GeneGalleon's version, source revision, and image digest.
Missing software is fetched only when assembly, BUSCO, or quantification is needed;
`plan` and reuse-only builds do not fetch it. To download in advance:

```bash
./run_build.sh fetch-software
```

Software is cached under `genegalleon.cache_dir`; `repository` and `image` allow
local overrides. Incompatible product conditions are reported as conflicts;
use a separate `store` for a deliberate rebuild with different conditions.

## Build through mapping

```bash
./run_build.sh plan
./run_build.sh prepare --name build001
./run_build.sh submit --build builds/build001 --until busco --dry-run
./run_build.sh submit --build builds/build001 --until busco
./run_build.sh status --build builds/build001

# After inspection, finish the build:
./run_build.sh submit --build builds/build001 --until mapping
```

Endpoints are `assembly`, `busco`, `quant`, and `mapping` (default). Each includes
missing prerequisites; assembly includes longest-CDS generation. Every species
remaining after run exclusions must finish CDS, full BUSCO, quantification, and
mapping. BUSCO acceptance thresholds apply later, in analysis.

The plan reports `reuse`, `pending`, or `conflict`. ODB reuse is confirmed after
translation provides protein hashes. Conflicts stop preparation/submission.
For a pilot, add `--species-list pilot.txt` to `submit`; later submit without it
to finish all species. An incomplete pilot does not start mapping or publish a
completed build.

After successful mapping and validation, build publishes
**`builds/<id>/products/`** and `completed.json`. New CDS/BUSCO/quant files are not
written to the top-level `input/`; see the [output layout](outputs.md#directory-layout).
If mapping was run separately, `run_build.sh complete --build builds/<id>`
validates and publishes it. Normal submission does this automatically.

## Run an analysis

In `analysis.yaml`, choose the BUSCO threshold, species selection, traits, and
optional branches. Set `build` there or pass `--build`:

```bash
./run_analysis.sh plan --build builds/build001
./run_analysis.sh prepare --build builds/build001 --name analysis001
./run_analysis.sh submit --analysis analyses/analysis001 --dry-run
./run_analysis.sh submit --analysis analyses/analysis001
./run_analysis.sh status --analysis analyses/analysis001
```

Analysis reuses verified proteins and mappings; it never runs assembly, BUSCO,
CDS translation, or ODB-mapper. Missing or modified build products cause an error.
Multiple analyses can share one build, inheriting its lineage, genetic code, and
ODB node. Optional analysis outputs are not automatically shared across analysis IDs.

`exclude_species` removes exact species IDs **before computation**.
`inputs.species_trait` supplies traits (`null` when unused); other auxiliary inputs
are described in [configuration](configuration.md). Changed inputs or scientific
settings require a new analysis name. Names beginning with `build_` are reserved.

The default target `all` collects `results/<analysis>/phenoradar_inputs/` on success.
Use `submit --target phylogeny` or another [analysis target](running.md#targets)
for partial execution; collect again after additional branches finish.

## Slurm and retries

`submit` validates inputs, submits jobs, and returns without waiting for completion.
CDS translations are cached by CDS content and genetic code; `register --products`
also imports the bundle's verified translations.

Assembly/BUSCO/quant use species arrays; mapping and analysis use Snakemake
controllers with separate rule jobs. See [execution and resources](running.md)
for concurrency, time limits, and resource overrides.

Before retrying, inspect `status`, Slurm jobs, and `jobs/logs/` under the prepared
build/analysis. Queued/running or unresolved submissions block overlapping retries.
A controller timeout can leave worker jobs running; inspect those jobs too.

## Failure and recovery behavior

Failures do not automatically exclude species. Successful, validated species
stages are registered; failed or missing stages remain pending. Resubmit the
**same build** after inspecting the queue to retry them. Completed work is reused.
Failure details are recorded in `jobs/status/` and logs; a hard kill may leave an
attempt marked `running` even though the job has stopped.

Array batches/phases depend on predecessor success (`afterok`), so a failure can
block later jobs. Let remaining jobs finish or cancel them before resubmitting.
Retries are explicit, not an automatic retry loop.

Reads and temporary work are retained, but stage retry is not always checkpoint
resume. The pinned GeneGalleon uses AMALGKIT `getfastq --redo no` to reuse valid
read state; partially downloaded bytes are not guaranteed reusable. Its rnaSPAdes
path restarts a failed species' assembly from scratch. Successful species are
unaffected. For repeated timeouts or memory errors, adjust
[resources](running.md#resource-budgets) before retrying.

## Manually excluding unusable accessions

Edit [config/excluded_accessions.tsv](../config/excluded_accessions.tsv):

```tsv
accession	reason
SRR123456	download_failed_repeatedly
LOCAL_BAD1	unusable_reads
```

`accession` matches the metadata **`run`** column exactly; `reason` is optional.
Additional review columns are retained. Project/experiment/species IDs do not
exclude their associated runs. The list is selected by `build.yaml`'s
`excluded_accessions`; `null` disables it.

Exclusions apply before cached-product lookup or FASTQ requirements in a new
build. Excluding a species' only run removes it from that build and its analyses.
Prepare a **new build ID** after edits; previous builds and jobs are unchanged.
Eligible completed products remain reusable. The list does not retroactively
invalidate an existing CDS reference from another run of the same species.

Each build preserves `source_metadata.tsv`, the exclusion list, effective
`metadata.tsv`, and `excluded_runs.tsv` with reasons. Keep decisions for accessions
absent from current metadata too, so future metadata preparation can avoid them.

## Updating species and importing existing work

Edit metadata and prepare a **new build ID**. Added species run missing work;
removed species leave the new outputs. Historical builds and caches remain.
Changing only a run can reuse CDS/BUSCO/mapping and quantify the new run.
Mappings are stored per species; updating membership links only the selected tables.
There is no combined mapping database to rebuild.

**Fresh builds need no registration.** To import legacy CDS/BUSCO/quant files,
use the [import layout](inputs.md#importing-existing-products):

```bash
./run_build.sh register --input-dir imports/legacy --metadata input/metadata.tsv
./run_build.sh register --odb-results imports/old_odb --odb-only
```

Registration validates existing work without recomputation. Missing full BUSCO
tables or quantification remain pending; a BUSCO summary alone is insufficient.
Keep source files referenced by the species store. With multiple CDS references,
select `reference_id` in metadata explicitly.

ODB registration copies or hard-links a validated snapshot into the automatic
cache, preserving its source. Repeat `--odb-results` for multiple snapshots.
No `odb.existing_results` setting is needed. See [ODB imports](references.md#reusing-existing-odb-results).

## Copying a completed build to another project

Copy **`builds/<id>/products/` as a unit**, including its manifest, metadata,
CDS/BUSCO/quant, proteins, ODB mappings, and provenance. Put it inside the destination
project, for example `imports/baseline/`:

```bash
./run_analysis.sh prepare --build imports/baseline --name analysis001
./run_analysis.sh submit --analysis analyses/analysis001
```

Analysis needs no registration or original reads/workspaces. The destination
supplies analysis settings, traits, software, and shared reference resources.
Both `builds/<id>/` and its `products/` directory are valid build arguments.
Products are immutable; publication may use hard links, so never edit them in place.

To seed a **new build with additional species**, register the copied bundle once:

```bash
./run_build.sh register --products imports/baseline
./run_build.sh plan
./run_build.sh prepare --name expansion001
./run_build.sh submit --build builds/expansion001 --until mapping
```

Keep the copied bundle: species registration references it. Use compatible
lineage, genetic code, and ODB node settings. Destination run exclusions apply.
Registered translations and species mapping tables are reused by new builds.
Copying products does not resume an interrupted build.

Legacy bundles containing an ODB SQLite database can be imported without using
that database:

```bash
./run_build.sh register --input-dir imports/legacy --metadata imports/legacy/metadata.tsv \
  --odb-results imports/legacy/odb
```

Prepare a new build afterward. Legacy annotations are split and validated once;
subsequent builds reuse the species tables.

## Migrating old configurations

`config/config.yaml` and `config/dataset.yaml` have been replaced by the two phase
configs. To convert saved settings without changing existing data:

```bash
python workflow/scripts/migrate_phase_config.py \
  --legacy-dataset config/dataset.local.yaml --legacy-config saved-config.yaml \
  --build-output config/build.migrated.yaml --analysis-output config/analysis.migrated.yaml
```

A saved `datasets/<old-id>/dataset.json` also works as `--legacy-dataset`.
Review the converted settings, apply them to `build.yaml`/`analysis.yaml`, import
old products if needed, and prepare a new build. Older schema-1 builds are not
portable bundles; resume their jobs with the original checkout. `run_dataset.sh`
remains an alias for the build interface; the old combined `--until all` is removed.
