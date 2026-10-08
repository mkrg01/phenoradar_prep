# Curated datasets

Retrieve and curate RNA-seq metadata with AMALGKIT in a pinned GeneGalleon
container, then select representative runs. Metadata updates and database
builds use separate commands.

| Dataset | Scope |
| --- | --- |
| [angiosperm_leaf](angiosperm_leaf/README.md) | One leaf RNA-seq run per original NCBI taxid |

## Selection policy

New candidates must pass AMALGKIT curation for the configured `sample_group`.
For each original `taxid`:

1. Apply exact accession exclusions. Other runs from the same taxid and
   BioProject remain eligible.
2. Keep the run in `accepted_samples.tsv`, even when a larger run appears.
   Exclusions take precedence.
3. Otherwise, choose the greatest positive `total_bases`, breaking ties by
   ascending accession.
4. Apply any optional `overrides.tsv` choice from the eligible candidates.

Review changes to accepted samples in `selection.tsv`. Their adopted names
and sample IDs are retained across metadata updates.

## Recording results

Use the build name printed by `submit` in place of `{build_name}`:

```bash
./run_metadata.sh record --build results/{build_name}
```

`record` applies the BUSCO completeness threshold in `selection.yaml`
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

## Choosing where to restart

| Purpose | Command |
| --- | --- |
| Refresh NCBI metadata and curate | `sbatch run_metadata.sh update` |
| Reselect from saved candidates | `sbatch run_metadata.sh update --metadata <curated.tsv>` |
| Retry failed or unfinished processing | `./run_build.sh submit --build results/<existing-build>` |

Submit from the repository root. Metadata jobs request 4 CPUs, 128 GB, and three
days; override these with the `sbatch` options `--mem` and `--time`.
Logs go to `slurm-metadata-<jobid>.out`. Use `./run_metadata.sh` for direct execution.

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
| `selection.yaml` | Search, tissue, paths, and adoption BUSCO threshold |
| `build.yaml` | Build settings and pinned GeneGalleon software |
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
lineage, and OrthoDB node. Set `selection.yaml`'s `sample_group` to a group allowed
by `select_rules.tsv`. Start with an exclusion table containing an `accession` header.
Pass `--config datasets/<name>/selection.yaml` to metadata commands and the
corresponding `build.yaml` to build commands.
