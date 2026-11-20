# Curated datasets

Retrieve and curate RNA-seq metadata with AMALGKIT in a pinned GeneGalleon
container, then select representative runs. Metadata updates and database
builds use separate commands.

| Dataset | Scope |
| --- | --- |
| [angiosperm_leaf](angiosperm_leaf/README.md) | One leaf RNA-seq run per original NCBI taxid |

## Selection policy

Select one run per submitted TaxID by default. Apply exact accession
exclusions and reviewed overrides, then keep the representative in
`accepted_samples.tsv` if its scientific name and TaxID still match. For each
remaining TaxID, choose the eligible run with the largest positive
`total_bases`, breaking ties by accession. New candidates must pass AMALGKIT
curation for the configured `sample_group`.

When an accepted run has a new name or TaxID, it competes again under the
updated identity. The old TaxID can select a replacement. Without a selection
taxonomy snapshot, SRA and BioSample must both report the new TaxID. With a
snapshot, a move to a different species is checked at species rank.
`selection.tsv` shows the old and updated identities. Publishing a candidate
removes the old acceptance. If the run is selected again, `record` registers it
after processing and QC under the updated identity.

Set `metadata_selection.taxonomy_database` in `build.yaml` to select one run per
species in a fixed NCBI taxonomy SQLite snapshot. Set `taxonomy.database` to the
same path for the build. Overrides use resolved species TaxIDs. The selection
records submitted names and TaxIDs as provenance.

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

Taxonomy screening is enabled by default in `build.yaml`:

```yaml
taxonomy_screen:
  enabled: true
  method: angiosperms353
  max_registered_support: 0.10
  batch_size: 24
```

`enabled: true` makes the screening commands available; database builds and
analysis submissions do not schedule the screen automatically. Run it explicitly
after sample processing. The default method,
`angiosperms353`, uses Kew v4.0 and the build's frozen NCBI taxonomy when
available. Install its [environment](angiosperms353.environment.yaml).

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

For other eukaryotes, set `taxonomy_screen.method: orthodb` and activate the
[OrthoDB screen environment](orthodb.environment.yaml). The same commands apply.
The screen uses the build's `busco.lineage` by default. Set
`taxonomy_screen.busco_lineage` to use a different dataset, such as
`vertebrata_odb12.2` for vertebrates. Use `--busco-dataset` to supply a local
BUSCO lineage directory. Otherwise the screen uses a cached copy or downloads it.

The OrthoDB method uses the marker proteins supplied with that BUSCO dataset.
For each species, it keeps the assembly with the most single-copy markers and
uses only those markers. It reuses compatible BUSCO tables or runs
BUSCO on proteins translated from the saved CDS. Each marker compares
amino-acid identity over at least 50 shared residues and 70% of the query.
Ties cast no vote. The screen selects the closest of genus, family, and order
with at least 50 votes and flags registered support strictly below 10%.

Use `--cluster` to run one Slurm or AGE array task per batch. The task recovers
markers, searches references, and counts votes. Concurrency follows
`execution.total_limits` and each task requests
`execution.per_job_resources.taxonomy_screen` (4 CPUs, 16 GB, 8 hours by default).
Without `--cluster`, batches run sequentially. Set `taxonomy_screen.batch_size`
to change the number of samples per batch. Rerun with the same output and
scientific settings to reuse completed batches and recovered markers after
a failure. Batch receipts record recovery and search/voting times. The final job
writes the combined report once all batches succeed.

Taxonomy commands, batch workers, and report jobs limit implicit NumPy/BLAS
library pools to one thread, including inside Captus. This prevents library
initialization from exceeding scheduler memory limits on large nodes. Captus,
BUSCO, and MMseqs still receive the configured CPU count for their own workers.

Use `--taxonomy` to supply a frozen snapshot and `--runs` to screen and review
the same subset. On Apple Silicon, specify a working native `--blat` if the
packaged executable cannot find its system libraries. A completed
screen has `summary.tsv`, per-marker `loci.tsv`, and `screen.json` with inputs,
versions and settings. All screened runs remain in `summary.tsv`, including runs
with fewer than 50 recovered loci: their sequences are still searched and their
available vote shares are retained, but they cannot receive a discordance call.
`taxonomy_summary.tsv` is a readable table with one row per run at its assessment
rank, including recovered/comparable locus counts, reported/leading-alternative
vote percentages, and the assessment. `taxonomy_table_notes.txt` defines the
columns and thresholds. Missing values are `NA`; percentages have two decimal
places. `summary.tsv` retains machine-readable values, and `loci.tsv` contains
marker evidence at the assessment rank. New screens do not generate all-rank
tables. `registered_support.svg` pools all assessed samples;
`registered_support_by_rank.svg` separates them by assessment rank, from genus
at the top to order at the bottom. Both show the cutoff and counts of assessed,
flagged and unassessed samples. Inconclusive runs are counted as unassessed,
not as zero-support samples. `exclusion_candidates.tsv` lists flagged accessions.
`record` includes the evidence when passed the screen.
`--exclude-taxonomy` proposes accession exclusions at the evaluated rank in
the existing editable plan. Applying the plan persists those decisions for the next
build. Check the plan before applying, including the existing read/scratch cleanup.

After changing the threshold in `build.yaml`, regenerate calls, the one-row
summary table and both plots from saved votes without searching sequences again.
Reports of older screens retain their original assessment ranks; replotting does
not apply the new order limit or replace an older inconclusive rank. To apply
the current rank-selection procedure, run a new screen in a new output directory.

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

Kew reference names are matched case-insensitively to scientific names and
synonyms in the same NCBI Taxonomy used for samples. Merged TaxIDs are resolved,
and candidates must have the queried rank and belong to Embryophyta (TaxID 3193).
Species names are tried first, followed by genus names and then Kew family
annotations only when the preceding rank has no eligible matches. A single
candidate supplies its NCBI lineage. Multiple candidates are assigned to their
lowest common ancestor, without further name matching or Kew family filtering.
References are used only at the assigned node and above; unmatched references
are excluded. Assignments above family remain available through order; references
assigned only above order are omitted. The screen records reference assignments
once and shares them across all batches.
References from the same sequencing run or reported species remain eligible by
default. Use explicit masks for validation comparisons when needed.
To refresh NCBI Taxonomy, prepare a snapshot at a new path and pass it with
`--taxonomy`. The original sample TaxID and accession remain in the source
database.

Captus runs one Scipio pass to recover Angiosperms353 CDS from single contigs. Loci with distinct
retained copies are skipped. MMseqs2 aligns the same locus only. For each locus
and rank, the alignment with the highest raw score is retained within each taxon.
Alternative taxa with the highest raw score are compared with the registered
taxon at the same A/C/G/T positions, requiring at least
150 positions and 70% of the query CDS. The highest identity wins. If several
taxa share the highest identity, they split one vote equally, including the
registered taxon when it ties for first.

Both methods assess ranks along the sample's NCBI lineage from genus through
order, stopping at the first rank with at least 50 voting loci. All available
loci are counted at that rank; assessment does not stop after the first 50 loci.
The Kew method includes named intermediate ranks and split votes for ties;
OrthoDB counts only untied loci and considers genus, family, and order. Sequence
searches run once per locus. Only unresolved samples advance to a higher rank,
reusing compressed search hits in task scratch, which is removed after the batch.

Each method emits one result per run. If no rank has enough voting loci, the
last available rank through order is reported, usually order. A lack of any
eligible registered-taxon reference among recovered loci gives `no_reference`;
otherwise the result is `insufficient_comparison`. A run with no recovered loci
is also inconclusive. Zero comparable loci gives an unavailable vote share
(`NA`), not zero support. Inconclusive and unassessed runs remain eligible for
the filtered database. At an adequately sampled assessment rank, a registered
vote share strictly below `taxonomy_screen.max_registered_support` (default 10%)
flags the run for review. Assessment then stops regardless of the outcome.
`--max-registered-support` overrides the config for a single run. `registered_support`
and `discordance_support` are the registered and combined alternative vote
shares and add up to 100% of the votes cast. Loci without a valid comparison
cast no vote. `support` reports the leading alternative’s share.

For validation, `--samples` accepts a TSV with `run`, `taxid` and `cds` columns.
`source_runs` records the query's source accessions. `exclude_accessions`
explicitly masks reference accessions.
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
