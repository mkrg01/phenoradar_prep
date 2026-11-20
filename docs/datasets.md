# Build a database and run analyses

[Documentation](index.md) · [Configuration](configuration.md)

[Build](#build-a-reusable-database) · [Analyze](#run-an-analysis) ·
[Reuse](#reusing-completed-databases) · [Reprocess samples](#reprocess-selected-samples) ·
[Update OrthoDB](#update-orthodb-for-all-samples)

A **build** produces a reusable sample database: assembly, BUSCO, quantification,
OrthoDB mapping, and OG expression. An **analysis** uses that database with its
own sample selection, traits, and optional analyses.

## Build a reusable database

Complete [installation](running.md#installation-and-normal-execution), prepare
[metadata](inputs.md), and edit [config/build.yaml](../config/build.yaml).
Set `reuse_from: null` for your first build, or [reuse existing results](#reusing-completed-databases).

Use the build directory name printed by `submit` for `{build_name}` below.

```bash
./run_build.sh plan
./run_build.sh submit
./run_build.sh status --build results/{build_name}
```

`submit` saves inputs and settings, verifies reusable products, then
submits the required jobs through Slurm or AGE. Add `--prepare-locally` to
import already verified inputs on the submit host.

New builds first check every candidate TaxID against one frozen taxonomy snapshot,
including builds that reuse all sample products. `taxonomy.missing_taxid: exclude`
omits unresolved IDs before import or sample dispatch; merged IDs remain included
with their original identity. Use `error` to stop instead, or `null` for legacy
selection. `taxonomy.database` can select an existing SQLite snapshot.

Review `taxonomy_selection.tsv` and `taxonomy_excluded_runs.tsv` in the build.
`taxonomy_selection.json` pins the decisions and reference checksum for resumes;
mapping uses that same snapshot. Source metadata and manual exclusions are unchanged.
Published `excluded_runs.tsv` reflects the final decisions; taxonomy records are
copied into `provenance/`.
An unreadable reference or an empty selected set stops the build. Plans show cached
decisions, or report taxonomy eligibility as pending without downloading a reference.

The completed database is `results/{build_name}/database/`. Builds record BUSCO
scores for every sample; filtering belongs to the analysis.

### Sample processing

Each run is assembled and quantified separately. With the supplied GeneGalleon
source and settings:

- AMALGKIT retrieves reads and normally applies fastp with a minimum retained
  read length of 25 nt. A fallback can recover original public FASTQs directly
  without fastp. The worker log identifies this route; while getfastq outputs
  are retained, `getfastq_completion.json` records `read_source: public-original`.
- rnaSPAdes assembles the reads. Inputs above 30 Gbp are subsampled to about
  30 Gbp with seed `11`.
- cdskit searches all six frames, preferring complete ORFs to partial candidates.
  Identical CDSs are collapsed, incomplete codons are padded, and one longest
  CDS is kept per rnaSPAdes gene group.
- BUSCO assesses the representative CDS set in transcriptome mode using
  `busco.lineage`. The analysis threshold defaults to 50% complete BUSCOs,
  counting single-copy and duplicated hits.

The [Kew/OrthoDB screen](../datasets/README.md#sample-taxonomy-screening) is a
separate step after sample processing. Reused samples keep their recorded
generating conditions; see the [sample processing table](#sample-processing-table).

## Plans and progress

Build `plan` shows sample and mapping reuse and saves `reports/plan/plan.tsv`.
It checks compatibility and checksums, reusing checks only for unchanged files.
`status` shows progress and log paths for failed jobs;
`status --verify` also checks product integrity.

| Option | Result |
| --- | --- |
| `plan --outdir DIR` | Save `plan.tsv` and `needs_attention.tsv` in `DIR` |
| `plan --split` | Also group reports by reuse, pending, conflict, and exclusion |
| `status --outdir DIR` | Save `status.tsv` and `failed.tsv` in `DIR` |
| `status --json` | Print the detailed report as JSON |

## Run an analysis

Edit [config/analysis.yaml](../config/analysis.yaml) to choose samples, the trait
column, and optional analyses. For OG expression only, set `alignment.enabled`
and `kegg.enabled` to `false`, use `phylogeny.trees: []`, and disable its
dating, taxonomy check, and contrast pair flags. Set `inputs.species_trait: null`
if traits are unused.
Replace `{analysis_name}` with a name for your analysis.

```bash
./run_analysis.sh plan --build results/{build_name}
./run_analysis.sh submit --build results/{build_name} --name {analysis_name}
./run_analysis.sh status --analysis results/{build_name}/downstream/{analysis_name}
```

Analysis `plan` verifies database products and prints a summary. It saves
`plan.tsv` (all samples) and `excluded.tsv` in `reports/analysis/plan/`.
Use `--outdir DIR` to change the directory, relative to the repository root.
Reports are overwritten on each run; `--json` prints detailed JSON instead of
the summary while still saving TSVs. Selected counts precede tree-specific
trait filtering and representative selection.

Use `status` to check the controller, workers, and logs. Use `status --verify`
to check frozen inputs and completed outputs, including large files. Use
`submit --analysis PATH --dry-run` to preview pending work from a saved analysis.
A run is complete after output validation and publication of `completed.json`.

Analyses inherit the database's BUSCO lineage, genetic code, and OrthoDB reference.
The default target `all` runs enabled analyses and validates native outputs.
[Use those outputs directly in PhenoRadar](phenoradar_inputs.md).
Several analyses can share one database;
use an [analysis target](running.md#targets) to run a single branch.

## Saved settings and previews

Inputs and settings are saved at first submission. Config edits do not change
a saved run: resume with `--build` or `--analysis`, or create a new run for
changed software or scientific settings. Use
[resource overrides](running.md#resource-budgets) for job resources.

| Naming or preview option | Behavior |
| --- | --- |
| Default naming | New builds append a UTC timestamp; resume using the printed path |
| `--name NAME` | Use an exact name |
| `name_mode: fixed` | Repeated `submit` resumes the configured build |
| Build `--dry-run` | Save the build and preview jobs without submitting |
| Analysis `--dry-run` | Preview without saving or changing the analysis |

## Failure and recovery behavior

Check the run's `status` and `jobs/logs/`, then retry with its saved path:

```bash
./run_build.sh submit --build results/{build_name}
./run_analysis.sh submit --analysis results/{build_name}/downstream/{analysis_name}
```

Completed samples are skipped; incomplete samples resume from verified
checkpoints. Add `--resources retry.yaml` to adjust resources. Active or unresolved
Slurm jobs block duplicate submissions.

Missing or modified sample products are rebuilt when reads are available.
Conflicting metadata, settings, or shared manifests stop the build; check
`needs_attention.tsv` or `status`. For a damaged reference pinned by `reference_id`,
create a new build without that pin. Publication requires every included sample
to succeed; other samples can continue after one fails.

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

Edit metadata and submit a new build with `reuse_from` pointing to previous
results. Added samples run missing work; removed samples are omitted. Changing
a run creates a new sample. Previous builds remain available.

For a pilot, add `--species-list pilot.txt` with one species or sample ID per line.
Resume with the same `--build` path and without the list to finish all samples.
Use `--until quant` to finish sample processing without publishing a database.

## Reusing completed databases

Set `reuse_from` in [config/build.yaml](../config/build.yaml):

| Value | Reusable results |
| --- | --- |
| `results/previous_build/database` | Sample processing, translated proteins, mapping, and OG expression |
| `results/previous_build/build.json` | Registered sample processing checkpoints |
| `auto` | Sample processing checkpoints discovered under `results/` |
| `null` | Process all samples from scratch |

A list of paths is also accepted. Reuse requires matching sample identities and
scientific settings. Assembly, BUSCO, and quantification form one sample unit:

| Sample state | Behavior in a new build |
| --- | --- |
| Complete | Reuse across GeneGalleon releases; retain the generating version |
| Incomplete, same software | Resume verified checkpoints |
| Incomplete, updated software | Restart from assembly |

Completed samples and checkpoints are matched by run accession, scientific
name, TaxID, and sample IDs. A run selected under an updated name or TaxID
starts from assembly unless a compatible result already exists for that
identity. See the [metadata selection policy](../datasets/README.md#selection-policy).

Keep sources available until preparation finishes importing products; the new
build then uses its own cache. Complete products can be reused without FASTQs.
Published products may share hard links: **do not edit them in place**.

### Reprocess selected samples

For a fix affecting existing results, list exact sample IDs or run accessions:

```yaml
reuse_from: results/previous_build/database
rebuild_samples: [SRR123456, Species_name_SRR789012]
```

Submit a new build. Selected samples restart from assembly; other compatible
samples are reused. Unknown IDs are rejected. Remove any explicit metadata
`reference_id` for a sample being rebuilt.

For subsequent builds, clear `rebuild_samples` and point `reuse_from` to the fixed
database. `auto` stops on conflicting products from different builds.

### Update OrthoDB for all samples

Set `odb.version` in the build config, then submit a **new build**. Each sample's
recorded OrthoDB version is checked automatically:

| Existing result | Action |
| --- | --- |
| Same release and clade | Reuse mapping and OG expression |
| Different release or clade | Redo mapping and OG expression |
| No mapping yet | Run mapping and OG expression |

Compatible assembly, BUSCO, quantification, and translated proteins are reused.
Example for OrthoDB v12:

```yaml
reuse_from: results/previous_build/database
rebuild_samples: []
odb:
  version: v12
  ncbi_tax_id: 3193
  chunk_size: 50
```

```bash
./run_build.sh plan
./run_build.sh submit --name orthodb_refresh
```

If earlier builds used different releases, only samples with a different release
or clade are remapped. The finished database uses the configured release for all
samples. BUSCO lineage is independent; the mapper must support the selected release.

Resume this build with `submit --build results/orthodb_refresh`; completed chunks
remain reusable. Run downstream analyses against the new database to use its
updated OGs. See [reference storage and checks](references.md#preparing-and-verifying-the-reference).

### Sample processing table

New builds include `database/sample_processing.tsv`: one row per sample with
generating versions, methods, references, and BUSCO completeness.
To export it from an older database without rerunning analyses:

```bash
./run_build.sh export-processing --build results/{build_name} \
  --outdir reports/sample_processing
```

Reused samples keep their generating versions; missing historical records are
`unknown`. See the [column definitions](outputs.md#sample-processing-table).

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

Build preparation fetches GeneGalleon when sample processing needs it. To download
in advance:

```bash
./run_build.sh fetch-software
```

`version` and `revision` pin the mounted workflow source; `image_uri` pins the
runtime by OCI SHA256 digest. Image labels and file checksums are recorded in
receipts. Downloads are cached under `resources/software/genegalleon/`;
`repository` and `image` select local copies.

The workflow source and runtime image are pinned separately. Upstream updates
take effect only after these pins are deliberately changed for a new build.
Report tool versions from the image used for the sample; assembly tools and
downstream analysis tools can use different versions. The image's GeneGalleon
label identifies its build, while `version` identifies the mounted workflow
source. Reused samples retain their original versions and settings.
