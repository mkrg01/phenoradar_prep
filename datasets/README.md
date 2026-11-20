# Curated datasets

Retrieve and curate public RNA-seq metadata with AMALGKIT, then select and
maintain representative runs. Metadata updates and database builds use separate
commands.

| Dataset | Scope |
| --- | --- |
| [angiosperm_leaf](angiosperm_leaf/README.md) | One leaf RNA-seq run per submitted NCBI TaxID |

Complete [installation](../docs/running.md#installation-and-normal-execution)
and run commands from the repository root. `run_metadata.sh` defaults to
`datasets/angiosperm_leaf/build.yaml`; use `--config datasets/<name>/build.yaml`
for another dataset. See [builds and analyses](../docs/datasets.md) for database
submission, reuse, resources, and retries.

## Selection policy

Select one run per submitted TaxID by default. Apply accession exclusions and
reviewed overrides, then retain QC-approved representatives whose scientific
name and TaxID match `accepted_samples.tsv`. For each remaining TaxID, choose the
eligible run with the largest positive `total_bases`, breaking ties by accession.
New candidates must pass AMALGKIT curation for the configured `sample_group`.

To select by species in a fixed NCBI taxonomy snapshot, set
`metadata_selection.taxonomy_database` and the build's `taxonomy.database` to
the same SQLite file. See [taxonomy references](../docs/references.md#taxonomy-reference).

An accepted run whose name or TaxID changes competes again under the updated
identity. SRA/BioSample confirmation uses submitted TaxIDs, or species-level
identities when a taxonomy snapshot is supplied. Review the changes in
`selection.tsv`; `record` registers the run again after processing and QC.

## Choosing where to restart

| Purpose | Command |
| --- | --- |
| Retrieve, curate, and select metadata | `./run_metadata.sh update` |
| Reselect from a saved curated table | `./run_metadata.sh update --metadata CURATED.tsv` |
| Apply changed curation rules to saved raw metadata | `./run_metadata.sh curate --work work/ATTEMPT --metadata RAW.tsv`, then `update --metadata` with its output |

For Slurm, submit metadata jobs with `sbatch run_metadata.sh update` (4 CPUs,
128 GB, three days by default). For AGE, use
`qsub -cwd -V -pe def_slot 4 -l s_vmem=32000M,h_rt=259200 run_metadata.sh update`,
adjusting resource names for the site.

`update` publishes representative metadata and prints a timestamped work path
under `work/datasets/<name>/`. Add `--dry-run` to prepare a candidate for review,
then use the printed `accept` command to publish it. `--work NEW_PATH` selects a
work path; it must not already exist.

If the update reports `review_required`, inspect the candidate's `selection.tsv`
and resolve missing or changed representatives before publishing. Reselection
requires the curated table from before representative selection. Run `update`
to regenerate a missing dataset `metadata.tsv`.

## Recording results

After sample processing, use the build name printed by `submit`:

```bash
./run_metadata.sh record --build results/BUILD
```

`record` applies `metadata_selection.busco_threshold` (default `0.5`, counting
single-copy and duplicated BUSCOs):

| Sample result | Action |
| --- | --- |
| Assembly, BUSCO, and quantification complete; BUSCO passes | Register in `accepted_samples.tsv` |
| BUSCO complete; new run below the threshold | Add to `excluded_accessions.tsv`, even if quantification failed |
| Incomplete processing | Keep eligible for retry |

Use `--dry-run` to preview decisions and cleanup, or `--runs SRR123 SRR456` to
review selected accessions. Metadata selection and the build must use the same
exclusion table.

### Abandoning incomplete samples

After deciding to stop retrying:

```bash
./run_metadata.sh record --build results/BUILD --exclude-failed
```

Inactive incomplete runs are excluded and their downloaded reads and temporary
work are removed. Original inputs, saved products, checkpoints, and logs are
kept. Active or unresolved jobs stay on hold. Rerun to retry failed cleanup.
Then run `update` and prepare a new build to select replacements and reuse
completed products.

## Optional review plan

To review or edit decisions before applying them:

```bash
./run_metadata.sh record --build results/BUILD \
  --exclude-failed --plan work/reviews/BUILD.tsv
./run_metadata.sh record --apply work/reviews/BUILD.tsv
```

Edit only `action` (`accept`, `exclude`, `hold`) and `reason`, or remove rows to
leave them unchanged. Keep the generated JSON files with the TSV. Creating a
plan changes no records and deletes no files; applying it verifies the saved
inputs. Reapply the unchanged plan to retry cleanup.

For [sequence-based taxonomy screening](../docs/taxonomy_screening.md), add
`--taxonomy-screen work/BUILD-taxonomy/screen.json --exclude-taxonomy` when
creating a plan. Applying the plan saves flagged accessions as exclusions for
subsequent metadata selection. Completed accepted runs are otherwise retained
when the BUSCO threshold changes; choose `exclude` explicitly to withdraw one.

## Dataset files

| File | Purpose |
| --- | --- |
| `build.yaml` | Metadata selection, QC, build settings, and pinned software |
| `select_rules.tsv` | AMALGKIT curation rules |
| `excluded_accessions.tsv` | Excluded accessions and reasons; `accession` required |
| `overrides.tsv` | Optional choices with `taxid`, `run`, and `reason` |
| `accepted_samples.tsv` | QC-approved runs registered by `record` |
| `metadata.tsv` | Generated representative metadata; ignored by Git |
| `selection.tsv` | Selection decisions and changes to review |
| `provenance.json`, `accepted_samples.provenance.json` | Metadata and QC history |

Keep settings and reviewed records in Git. A build's `database/source_metadata.tsv`
preserves its input metadata. To use the pinned image's curation rules, run
`./run_metadata.sh rules --replace`.

## Adding another dataset

Copy a dataset's configuration and rules, then adjust the query, paths, tissue,
BUSCO lineage, and OrthoDB node. Choose a `metadata_selection.sample_group`
allowed by `select_rules.tsv` and create an exclusion table with an `accession`
header. Pass `--config datasets/<name>/build.yaml` to metadata and build commands.
