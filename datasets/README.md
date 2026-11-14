# Curated datasets

Retrieve and curate RNA-seq metadata with AMALGKIT in a pinned GeneGalleon
container, then select representative runs. Metadata updates and database
builds use separate commands.

| Dataset | Scope |
| --- | --- |
| [angiosperm_leaf](angiosperm_leaf/README.md) | One leaf RNA-seq run per original NCBI taxid |

## Selection policy

New candidates must pass AMALGKIT curation for the configured `sample_group`.
For each original `taxid` (or resolved species TaxID when configured):

1. Apply exact accession exclusions. Other runs from the same taxid and
   BioProject remain eligible.
2. Keep the run in `accepted_samples.tsv`, even when a larger run appears.
   Exclusions take precedence.
3. Otherwise, choose the greatest positive `total_bases`, breaking ties by
   ascending accession.
4. Apply any optional `overrides.tsv` choice from the eligible candidates.

Set `metadata_selection.taxonomy_database` in `build.yaml` to a fixed NCBI taxonomy SQLite
snapshot to select one run per current species. Set `taxonomy.database` in the
build configuration to the same path. Overrides use resolved species TaxIDs.
Selection preserves the submitted name and TaxID and reports unresolved species.

A confirmed change to another species releases the previous accepted identity.
The run competes with eligible candidates in its current species, and the old
species can select a replacement. Selection records the SRA and BioSample
confirmation and retires the old acceptance when the candidate is published.

For each NCBI refresh, prepare a new snapshot at a new path, then point both
configurations to it before selecting or building:

```bash
python workflow/scripts/prepare_taxonomy.py --destination resources/taxonomy/2026-11-13/taxa.sqlite
```

Existing snapshots are reused when configured. A saved build checks the
snapshot it used so that retries use the same classification.

## Recording results

Use the build name printed by `submit` in place of `{build_name}`:

```bash
./run_metadata.sh record --build results/{build_name}
```

`record` applies `metadata_selection.busco_threshold` from `build.yaml`
(`(single + duplicated) / total`, default `0.5`):

| Sample result | Action |
| --- | --- |
| Assembly, BUSCO and quantification complete; BUSCO passes | Register in `accepted_samples.tsv` |
| BUSCO complete; new run below the threshold | Add to `excluded_accessions.tsv`, even if quantification failed |
| Incomplete processing | Keep eligible for retry |

Add `--dry-run` to preview actions, stages, errors, BUSCO scores and cleanup sizes
without changing records or deleting files. Use `--runs SRR123 SRR456` to limit
the operation to specific accessions. The selection and build configurations must
use the same exclusion table.

## Abandoning incomplete samples

After deciding to stop retrying, exclude the remaining inactive incomplete runs
and clean up in one command:

```bash
./run_metadata.sh record --build results/{build_name} --exclude-failed
```

This includes unstarted runs and runs with unfinished quantification. Queued,
running, unresolved or conflicting samples stay on hold. Add `--dry-run` for an
optional preview; execution uses the current sample state.

Excluded runs lose downloaded FASTQ/SRA copies and temporary work, including in
builds retaining intermediates. Original inputs, saved products, restart
checkpoints and logs are kept. Cleanup failures leave exclusions in place;
repeat the command to retry cleanup.

Next, run `update` to select replacement runs, then prepare a new build to reuse
compatible completed products. Taxids without a replacement are dropped.
Existing builds keep their original sample selection.

## Optional review plan

For an editable list of decisions, save an optional plan:

```bash
./run_metadata.sh record --build results/{build_name} \
  --exclude-failed --plan work/reviews/{build_name}.tsv
./run_metadata.sh record --apply work/reviews/{build_name}.tsv
```

Before applying, edit only `action` (`accept`, `exclude`, `hold`) and `reason`, or
remove rows to leave them alone. Keep the generated JSON files with the TSV.
Plan creation changes no dataset records and deletes no files. Applying it checks
the saved decisions; changed metadata or selected samples require a new plan.
Reapply the unchanged plan to retry cleanup.

Completed accepted runs are retained if the BUSCO threshold changes. To withdraw
one, explicitly choose `exclude` in a plan. Use a new plan filename for each review.

## Sample taxonomy screening

Enable the optional check in the dataset's `build.yaml`:

```yaml
taxonomy_screen:
  enabled: true
  method: angiosperms353
  max_registered_support: 0.10
  batch_size: 24
```

The default is `enabled: false`. Normal metadata updates and builds run
independently of the screen. The `angiosperms353` method uses registered CDS
assemblies and the build's frozen NCBI taxonomy if available. Install its
[environment](angiosperms353.environment.yaml). Preparation downloads Kew
v4.0 and a reusable taxonomy snapshot. Screening can also fetch missing markers.

```bash
conda env create -n phenoradar-taxonomy -f datasets/angiosperms353.environment.yaml
conda activate phenoradar-taxonomy
./run_metadata.sh prepare-taxonomy
./run_metadata.sh screen-taxonomy --build results/<build> \
  --output work/<build>-taxonomy --threads 8 --ram 16
./run_metadata.sh record --build results/<build> \
  --taxonomy-screen work/<build>-taxonomy/screen.json \
  --exclude-taxonomy --plan work/<build>-taxonomy-review.tsv
./run_metadata.sh record --apply work/<build>-taxonomy-review.tsv
```

Use `--cluster` to run one Slurm or AGE array task per batch. The task recovers
markers, searches references, and counts votes. Concurrency follows
`execution.total_limits` and each task requests
`execution.per_job_resources.taxonomy_screen` (4 CPUs, 16 GB, 8 hours by default).
Without `--cluster`, batches run sequentially. Set `taxonomy_screen.batch_size`
to change the number of samples per batch. Rerun with the same output and
scientific settings to reuse completed batches and recovered markers after
a failure. Batch receipts record recovery and search/voting times. The final job
writes the combined report once all batches succeed.

Use `--taxonomy` to supply a frozen snapshot and `--runs` to screen and review
the same subset. On Apple Silicon, specify a working native `--blat` if the
packaged executable cannot find its system libraries. A completed
screen has `summary.tsv`, per-marker `loci.tsv`, and `screen.json` with inputs,
versions and settings. `registered_support.svg` pools all assessed samples.
`registered_support_by_rank.svg` separates them by the rank used for each
decision. Both show the cutoff and counts of assessed, flagged and unassessed
samples. `exclusion_candidates.tsv` lists the flagged accessions for review.
`record` includes the evidence when passed the screen.
`--exclude-taxonomy` proposes accession exclusions at the evaluated rank in
the existing editable plan. Applying the plan persists those decisions for the next
build. Check the plan before applying, including the existing read/scratch cleanup.

After changing the threshold in `build.yaml`, regenerate the calls and both
plots from saved votes without searching the sequences again:

```bash
./run_metadata.sh report-taxonomy --screen work/<build>-taxonomy/screen.json \
  --output work/<build>-taxonomy-new-threshold
```

For a copied portable database, pass `--database <database-dir>` to `screen-taxonomy`
with the same output and resource options. Use `--taxonomy` to supply a frozen
snapshot. Review `exclusion_candidates.tsv` and publish a separate database
using all the discordant accessions in the completed report:

```bash
./run_metadata.sh filter-taxonomy --database <database-dir> \
  --screen work/<build>-taxonomy/screen.json --output work/<build>-filtered/database
```

The filtered database reuses the original CDS, BUSCO, quantification, proteins,
expression and species-local OrthoDB mappings for retained accessions. The source
database remains available.
To look for
replacement runs, apply the exclusions during metadata selection before the new
build. Turning the screen off leaves previously approved exclusions in effect.

Kew reference species names are matched to scientific names and synonyms in the
frozen NCBI Taxonomy. A unique species TaxID supplies the reference lineage.
Unmatched names use the genus and family, and ambiguous species names can
contribute only at a shared, compatible family or higher rank. The screen
records the reference assignments once and shares them across all batches.
To refresh NCBI Taxonomy, prepare a snapshot at a new path and pass it with
`--taxonomy`. The original sample TaxID and accession remain in the source
database.

Captus runs one Scipio pass to recover Angiosperms353 CDS from single contigs. Loci with distinct
retained copies are skipped. Same-run reference sequences are removed, including
merged Kew run IDs. MMseqs2 aligns the same locus only. For each rank, the
highest-scoring eligible alignment for the registered taxon and all co-leading
alternative taxa are compared at the same A/C/G/T positions, requiring at least
150 positions and 70% of the query CDS. The highest identity wins. If several
taxa share the highest identity, they split one vote equally, including the
registered taxon when it ties for first.

The screen follows the sample’s frozen NCBI lineage from genus toward the
root and evaluates the first named rank with an unmasked reference in any
recovered marker. It emits one result per run. If the selected rank has
fewer than 50 comparable loci, the result is inconclusive. When no rank has
a reference, the result is unassessed. At every evaluated rank, a registered
vote share strictly below `taxonomy_screen.max_registered_support` (default 10%) flags the run for review.
`--max-registered-support` overrides the config for a single run. `registered_support`
and `discordance_support` are the registered and combined alternative vote
shares and add up to 100% of the votes cast. Loci without a valid comparison
cast no vote. `support` reports the leading alternative’s share.

For validation, `--samples` accepts a TSV with `run`, `taxid` and `cds` columns.
`source_runs` and `exclude_accessions` can mask additional known sources.
`--omit-same-species` uses the TSV's `scientific_name` column to mask exact
species-name matches. `--loci` limits the validation to a marker list.

## Choosing where to restart

| Purpose | Command |
| --- | --- |
| Refresh NCBI metadata and curate | `sbatch run_metadata.sh update` |
| Reselect from saved candidates | `sbatch run_metadata.sh update --metadata <curated.tsv>` |
| Retry failed or unfinished processing | `./run_build.sh submit --build results/<existing-build>` |

Submit from the repository root. Metadata jobs request 4 CPUs, 128 GB, and three
days; override these with the `sbatch` options `--mem` and `--time`.
Logs go to `slurm-metadata-<jobid>.out`. Use `./run_metadata.sh` for direct execution.

For AGE, the same script accepts `qsub -cwd -V -pe def_slot 4 -l
s_vmem=32000M,h_rt=259200 run_metadata.sh update` (adjust the PE and memory
resource names for the site). AGE ignores the `#SBATCH` defaults, so provide its
resources explicitly. Database and analysis submissions detect the scheduler
automatically through their normal launchers.

`update` creates `work/datasets/<name>/YYYYMMDDTHHMMSSZ/` at execution start
(UTC) and prints the path. Use `--work <new-path>` for a manual path under `work/`.
Existing run names, including same-second collisions, cause an error.

By default, `update` publishes representative metadata. Add `--dry-run` to
generate a candidate for review; acquisition and selection still run. After
review, use the printed `accept` command to adopt the candidate.

Run `update` to regenerate a missing dataset `metadata.tsv`. If it reports
`review_required`, inspect the candidate's `selection.tsv`. Restore missing runs
from a saved metadata snapshot, or edit `excluded_accessions.tsv` or
`overrides.tsv`, then rerun `update`.

Reselection needs the AMALGKIT-curated table **before representative selection**.
If curation rules changed, first run
`sbatch run_metadata.sh curate --work <new-attempt> --metadata <saved-raw.tsv>`;
after the job finishes, pass its output to `update --metadata`.

Dataset builds use `reuse_from: auto` to reuse verified sample stages. Resume
with `--build results/<printed-build-name>` to keep frozen inputs and completed
work; add `--resources <retry.yaml>` to adjust resources. See
[builds and recovery](../docs/datasets.md).

## Dataset files

| File | Purpose |
| --- | --- |
| `build.yaml` | Metadata selection, QC thresholds, build settings and pinned GeneGalleon software |
| `select_rules.tsv` | AMALGKIT curation rules |
| `excluded_accessions.tsv` | Excluded run accessions and reasons; `accession` required |
| `overrides.tsv` | Optional choices with `taxid`, `run`, and `reason` |
| `accepted_samples.tsv` | QC-approved runs registered by `record`; accession in `run` |
| `selection.tsv` | Selection decisions and changes to review |
| `provenance.json`, `accepted_samples.provenance.json` | Metadata and QC history |

To replace curation rules with the pinned image's rule set, run
`./run_metadata.sh rules --replace`.

Keep dataset settings and reviewed records in Git. Full `metadata.tsv` is ignored.
For an exact historical snapshot, use the completed build's
`database/source_metadata.tsv`.

## Adding another dataset

Copy a dataset's configuration and rules; adjust its query, paths, tissue, BUSCO
lineage, and OrthoDB node. Set `metadata_selection.sample_group` to a group allowed
by `select_rules.tsv`. Start with an exclusion table containing an `accession` header.
Pass `--config datasets/<name>/build.yaml` to metadata and build commands.
