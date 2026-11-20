# Sample taxonomy screening

[Documentation](index.md) · [Taxonomic monophyly](taxonomy_check.md) · [Builds and analyses](datasets.md)

Compare assembled sequences with Kew or OrthoDB references to identify samples
that disagree with their reported taxonomy. Screening works with ordinary
builds and portable databases as well as curated datasets. Run it explicitly
after sample processing; builds and analysis submissions do not schedule it.

## Setup

The default Kew method uses Angiosperms353 release 4.0. Create its environment:

```bash
conda env create -n phenoradar-taxonomy -f datasets/angiosperms353.environment.yaml
conda activate phenoradar-taxonomy
```

Set these options in the build config:

```yaml
taxonomy_screen:
  enabled: true
  method: angiosperms353
  max_registered_support: 0.10
  batch_size: 24
```

For other eukaryotes, use `method: orthodb` and
[`datasets/orthodb.environment.yaml`](../datasets/orthodb.environment.yaml).
This method uses marker proteins from the build's `busco.lineage`. Set
`taxonomy_screen.busco_lineage` to choose another lineage, or `--busco-dataset`
to supply its local directory.

## Run the screen

Run from the repository root and select the config explicitly:

```bash
./run_metadata.sh screen-taxonomy --config config/build.yaml \
  --build results/BUILD --output work/BUILD-taxonomy
```

Replace `--build results/BUILD` with `--database results/BUILD/database` to use a
portable database. References are cached automatically. To download them in
advance:

```bash
./run_metadata.sh prepare-taxonomy --config config/build.yaml
```

| Option | Purpose |
| --- | --- |
| `--cluster` | Submit one Slurm or AGE array task per batch |
| `--threads N`, `--ram GB` | Override each batch's CPU count and search memory budget |
| `--taxonomy PATH` | Use an existing NCBI taxonomy SQLite snapshot |
| `--runs SRR123 SRR456` | Screen selected accessions |
| `--min-loci N` | Minimum comparable loci for assessment (default: 50) |
| `--blat PATH` | Select a BLAT executable for Captus |

With `--build`, the screen uses the build's taxonomy snapshot when available;
otherwise it prepares one in the reference cache. Use `--taxonomy` to select the
snapshot for a portable database. Settings under
`execution.per_job_resources.taxonomy_screen` default to 4 CPUs, 16 GB, and
8 hours; cluster concurrency follows [execution limits](running.md#resource-budgets).
Without `--cluster`, batches run sequentially.

Retry with the same output path and scientific settings to reuse completed
batches and recovered markers. The combined report is written after all batches
succeed.

## Read the results

| Output | Contents |
| --- | --- |
| `taxonomy_summary.tsv` | Readable table: one row per run, assessment rank, locus counts, vote percentages, and decision |
| `taxonomy_table_notes.txt` | Column definitions and thresholds |
| `summary.tsv`, `screen.json` | Machine-readable results, inputs, versions, and settings |
| `loci.tsv` | Per-marker evidence at the assessment rank |
| `exclusion_candidates.tsv` | Flagged run accessions |
| `registered_support.svg` | Reported-taxon support across assessed samples |
| `registered_support_by_rank.svg` | Support grouped by assessment rank |

Ranks are assessed from genus through order. By default, the first rank with at
least 50 comparable loci determines the decision, using all comparable loci at that
rank. A reported-taxon vote share strictly below 10% flags a sample; change
`taxonomy_screen.max_registered_support` or use `--max-registered-support` to
set another cutoff.

All screened runs are reported, including samples with fewer than 50 recovered
loci. Their sequences are searched and available vote shares are retained, but
they receive no discordance call. If no rank has enough evidence, the result is
`no_reference` when no eligible reported-taxon reference exists, or
`insufficient_comparison` otherwise. Zero comparable loci gives `NA` support.
Inconclusive samples remain eligible when filtering a database. The readable
table reports percentages to two decimal places; plots distinguish assessed,
flagged, and unassessed samples.

In `summary.tsv`, `registered_support` is the reported-taxon vote share and
`discordance_support` combines alternative taxa. `support` is the leading
alternative's share.

To change the cutoff without repeating sequence searches, edit the config and
regenerate the report from saved votes:

```bash
./run_metadata.sh report-taxonomy --config config/build.yaml \
  --screen work/BUILD-taxonomy/screen.json --output work/BUILD-taxonomy-report
```

This uses the assessment ranks recorded by the screen.

## Exclude flagged samples

After reviewing the report, export a separate database excluding all flagged
accessions:

```bash
./run_metadata.sh filter-taxonomy --config config/build.yaml \
  --database results/BUILD/database --screen work/BUILD-taxonomy/screen.json \
  --output results/BUILD-filtered/database
```

The export reuses retained samples' CDS, BUSCO, quantification, proteins,
expression, and species-local OrthoDB mappings. The source database is preserved.
To exclude accessions from later builds, update the
[accession exclusion list](datasets.md#manually-excluding-unusable-accessions).
Curated datasets can apply the screen through an
[editable review plan](../datasets/README.md#optional-review-plan), which also
records exclusions for selecting replacement runs.

## How comparisons are made

### Kew Angiosperms353

Kew names are matched case-insensitively to NCBI scientific names and synonyms,
with merged TaxIDs resolved and candidates restricted to the queried rank within
Embryophyta (TaxID 3193). Matching tries species, then genus, then Kew family
annotations when the preceding rank has no eligible match. Multiple candidates
are assigned to their lowest common ancestor. References contribute only at
that node or above; unmatched references and assignments above order are omitted.

Captus uses one Scipio pass to recover loci from single contigs. Loci with
multiple distinct retained copies are omitted. MMseqs2 searches corresponding
loci. For each locus and rank, the highest-scoring alignment is retained within
each taxon. The reported taxon is compared with the alternative taxa having the
highest raw score, using the same unambiguous A/C/G/T positions. Comparisons
require at least 150 positions covering 70% of the query CDS. The highest
sequence identity wins one vote; tied taxa share the vote equally. Named
intermediate ranks between genus and order are included.

References from the query's sequencing run or species are eligible. Reference
assignments are shared across batches, and higher-rank comparisons reuse the
same sequence searches.

### OrthoDB/BUSCO

For each reference species, the screen retains the assembly with the most
single-copy BUSCO markers. Query markers come from compatible BUSCO tables, or
from BUSCO runs on proteins translated from the saved CDS. Each marker compares
amino-acid identity over at least 50 shared residues and 70% of the query.
Ties cast no vote. Assessment considers genus, family, and order.

### Validation inputs

Use `--samples PATH` instead of `--build` or `--database` to supply a TSV with
`run`, `taxid`, and `cds` columns. Optional `source_runs` records query accessions
and `exclude_accessions` masks reference accessions. `--omit-same-species` uses
the TSV's `scientific_name` column to mask exact species-name matches; `--loci`
limits comparisons to a marker list.
