# Reusable builds and independent analyses

[Documentation](index.md) · [Configuration](configuration.md)

`run_build.sh` prepares a reusable sample database through ODB mapping and OG expression.
`run_analysis.sh` selects samples from a completed build and produces expression,
alignment, phylogeny, and PhenoRadar outputs under `results/<build>/downstream/<name>/`.

The boundary separates work performed for each sample from work on a chosen sample set.
Assembly, BUSCO, quantification, mapping, and OG expression can be reused when
samples are added or removed. Selecting samples, combining tables, aligning genes,
and inferring trees belong to downstream conditions. A changed BUSCO threshold or
trait therefore needs a new downstream name, while retaining the same database.

Build stages metadata and local FASTQs in a separate GeneGalleon workspace for
each sample and schedules its assembly, BUSCO, and quantification stages.
GeneGalleon retrieves public reads through AMALGKIT `getfastq`, or processes the
staged local reads, and produces the assembly, longest CDS, BUSCO results, and
quantification. Build enables `getfastq` for the assembly and quantification
stages; valid reads retained in that workspace can be reused. Snakemake then
translates CDS, maps proteins to OrthoDB, and computes per-sample OG expression.

## Prepare a manually curated dataset

Maintain `input/metadata.tsv` with **all desired samples**, one run per row,
from NCBI or local FASTQs. Multiple samples may share a biological species. See [input formats](inputs.md) for required columns and
local-read paths. Keep known unusable runs in the
[manual exclusion list](#manually-excluding-unusable-accessions).

Edit `config/build.yaml` and `config/analysis.yaml` directly. These are the default
configs; the first `submit` freezes settings and inputs for each named run. Install the
[workflow environment](../environment.yaml) and prepare the
[workflow image](running.md#installation-and-normal-execution) before submitting mapping or analysis jobs.

## Pinned GeneGalleon source and SIF

Build settings pin GeneGalleon's version, source revision, and image digest.
`genegalleon.image_uri` is required and must use
`docker://<registry>/<image>@sha256:<64 lowercase hexadecimal characters>`.
Tag-only references and direct HTTPS SIF URLs are rejected; `image_sha256`
is no longer a supported setting.
Missing software is fetched only when assembly, BUSCO, or quantification is needed;
`plan` and reuse-only builds do not fetch it. To download in advance:

```bash
./run_build.sh fetch-software
```

Software is cached under `resources/software/genegalleon/` by default; no
`cache_dir` setting is needed in `config/build.yaml`. Matching cached source and SIF
files are verified by their recorded SHA256 checksums and reused without downloading
again. `repository` and `image` need not be set for automatic retrieval and reuse;
they remain optional local overrides. Incompatible product conditions are reported as conflicts;
use `reuse_from: null` and a new build name for a deliberate rebuild with different conditions.

## Build a reusable database

Set `name` in `config/build.yaml`, for example `angiosperm_leaf_20260925`.
Use a new name for changed inputs or settings; `submit --name NAME` overrides the config.
Names are literal: update the date yourself, adding `_v2` for same-day revisions.

```bash
mkdir -p results
./run_build.sh plan > results/build_plan.tsv
./run_build.sh submit
./run_build.sh status --build results/angiosperm_leaf_20260925

# To stop a new build at BUSCO for inspection:
./run_build.sh submit --name busco_check --until busco
# After inspection, finish that build:
./run_build.sh submit --name busco_check --until database
```

Endpoints are `assembly`, `busco`, `quant`, and `database` (default).
`mapping` remains a compatibility alias for the full database endpoint. Each includes
missing prerequisites; assembly includes longest-CDS generation. To publish a
completed database, every sample remaining after run exclusions must finish CDS,
full BUSCO, quantification, mapping, and per-sample OG expression. BUSCO acceptance
thresholds apply later, in downstream selection.

Review `results/build_plan.tsv` before submission (for example, with
`less -S results/build_plan.tsv`). The plan reports `reuse`, `pending`, or `conflict`.
ODB reuse is confirmed after
translation provides protein hashes. Conflicts stop preparation/submission.
For a pilot, add `--species-list pilot.txt` to `submit`; list biological species IDs
or exact sample IDs, one per line. Later submit without it to finish all samples.
An incomplete pilot does not start mapping or publish a completed build.

After successful mapping, OG expression, and validation, build publishes
**`results/<id>/database/`** and `completed.json`. New CDS/BUSCO/quant files are not
written to the top-level `input/`; see the [output layout](outputs.md#directory-layout).
If the low-level `database` target was run separately, `run_build.sh complete --build results/<id>`
validates and publishes it. Normal submission does this automatically.

## Run an analysis

In `analysis.yaml`, choose the BUSCO threshold, species selection, traits, and
optional branches. Set `build` there or pass `--build`.
Use a name that identifies the trait and date, such as `c4_photosynthesis_20260929`.
For this example, set `trait: C4` and supply a `C4` column in `inputs.species_trait`.
The `trait` setting selects the phenotype; `--name` names the output directory.
Add a suffix such as `_v2` when changing conditions on the same day.

```bash
./run_analysis.sh plan --build results/angiosperm_leaf_20260925
./run_analysis.sh submit --build results/angiosperm_leaf_20260925 --name c4_photosynthesis_20260929
./run_analysis.sh status --analysis results/angiosperm_leaf_20260925/downstream/c4_photosynthesis_20260929
```

Downstream reuses verified proteins, mappings, and per-sample OG expression;
it never runs assembly, BUSCO, CDS translation, or ODB-mapper. Missing or modified build products cause an error.
Multiple analyses can share one build, inheriting its lineage, genetic code, and
ODB node. Optional analysis outputs are not automatically shared across analysis IDs.

`exclude_species` accepts biological `species_id` values (all samples of that species) or exact analysis sample IDs **before computation**.
`inputs.species_trait` supplies traits (`null` when unused); other auxiliary inputs
are described in [configuration](configuration.md). Changed inputs or scientific
settings require a new downstream name within the collection. Different collections
may use the same downstream name.

The default target `all` collects `results/<build>/downstream/<analysis>/phenoradar_inputs/` on success.
Use `submit --target phylogeny` or another [analysis target](running.md#targets)
for partial execution; collect again after additional branches finish.

## Slurm and retries

The first `submit` saves inputs and settings, validates them, submits jobs, and
returns without waiting for completion. Repeating the same command reuses saved
conditions. Source edits require a new build or downstream name; resource changes
can be applied with `--resources`.

Add `--dry-run` to save the conditions and preview job scripts without submitting.
This also fixes the conditions for that name; remove `--dry-run` to submit them.
To save inputs without generating job scripts, the optional `prepare` command is
still available. `plan` previews current inputs and settings without saving a run.

Existing runs can also be submitted directly, without their source configs:

```bash
./run_build.sh submit --build results/angiosperm_leaf_20260925
./run_analysis.sh submit --analysis results/angiosperm_leaf_20260925/downstream/c4_photosynthesis_20260929
```

CDS translations are cached by CDS content and genetic code and reused automatically.

Assembly/BUSCO/quant use sample arrays; mapping and analysis use Snakemake
controllers with separate rule jobs. See [execution and resources](running.md)
for concurrency, time limits, and resource overrides.

Before retrying, inspect `status`, Slurm jobs, and `jobs/logs/` under the prepared
build/analysis. Queued/running or unresolved submissions block overlapping retries.
A controller timeout can leave worker jobs running; inspect those jobs too.

## Failure and recovery behavior

Failures do not automatically exclude samples. Successful, validated sample
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
path restarts a failed sample's assembly from scratch. Successful samples are
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

## Updating samples

Edit metadata and prepare a **new build ID**. Added samples run missing work;
removed samples leave the new outputs. Historical builds and caches remain.
Changing a run creates a new sample ID and requires its own assembly, BUSCO,
quantification, mapping, and OG expression. Existing samples reuse their own completed products.
Mappings are stored per sample; updating membership links only the selected tables.
There is no combined mapping database to rebuild.

## Reusing completed databases

Set `reuse_from` in `config/build.yaml` to the completed database to reuse:

```yaml
name: angiosperm_leaf_20260930
reuse_from: results/angiosperm_leaf_20260928/database/
```

Several databases can supply samples for one build:

```yaml
reuse_from:
  - results/leaf_set_a/database/
  - results/leaf_set_b/database/
```

Use `null` (or `[]`) for a fresh build. Only samples selected by the new metadata
and accession exclusions are considered. The pipeline verifies sample identities,
checksums, genetic code, BUSCO lineage, ODB node/reference, and recorded build
conditions before reuse. Duplicate paths are ignored. When databases contain the
same sample, its scientific products and conditions must agree; conflicting
products stop preparation regardless of list order. If matching sources differ
only in whether OG expression has been saved, the saved expression is reused.
Use a compatible set of sources or `null` for a deliberate rebuild.

CDS, BUSCO, quantification, proteins, mappings, and saved OG expression are reused.
A schema-4 sample bundle without saved expression is also accepted; its expression
is computed once. Missing samples run the normal pipeline. `plan` only inspects
sources and reports reuse, pending work, or conflicts.

Preparation stages verified products under the new build's `work/cache/` and
records their origin in `reuse.json` (also saved in `database/provenance/` on
completion). Prepared builds no longer need the source database paths to remain
available. Large immutable products may share hard links; never edit generated
files in place. Final products are published under the new build's `database/`.
There are no product writes to `migrations/` and no separately configured product
store or ODB cache. Old directories are not removed automatically.

An unfinished build resumes from its own `work/` and stage receipts under the same
name. To reuse it in a different build, finish and publish its database first.
No databases or unfinished builds are discovered implicitly by scanning `results/`.

## Copying a completed build to another project

Copy **`results/<id>/database/` as a unit**, including its manifest, metadata,
CDS/BUSCO/quant, proteins, ODB mappings, OG expression, and provenance. Put it inside the destination
project, for example `results/baseline/database/`:

```bash
./run_analysis.sh submit --build results/baseline --name c4_photosynthesis_20260929
```

Analysis reads the completed database directly and needs no original reads,
workspaces, or working caches. The destination supplies analysis settings, traits,
software, and shared reference resources.
Both `results/<id>/` and its `database/` directory are valid build arguments.
A renamed local database uses the destination collection name for downstream
outputs. Bundles elsewhere inside the project are also accepted; downstream goes
to `results/<manifest-build-id>/downstream/`. If that ID conflicts with a local
build, copy the bundle under a unique `results/<name>/database/` first.
Products are immutable; publication may use hard links, so never edit them in place.
