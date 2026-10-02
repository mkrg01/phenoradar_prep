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

After QC review, `record` registers runs with completed assembly, BUSCO, and
quantification and completeness `(single + duplicated) / total >= busco_threshold`
(default `0.5`). New runs with completed BUSCO below the threshold are added to
`excluded_accessions.tsv`, even before quantification finishes. Execution
failures and unfinished stages remain eligible for retry.

`record --dry-run` previews these decisions; `--runs` limits them to reviewed
accessions. Add other unusable runs to the exclusion table manually. Run
`update` afterward to replace excluded representatives, or drop their taxids
if no candidate remains. Source changes for retained runs are flagged in
`selection.tsv`.

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
| `rules/select_rules.tsv` | Effective AMALGKIT curation rules |
| `excluded_accessions.tsv` | Manual exclusions and recorded BUSCO failures; `accession` required |
| `overrides.tsv` | Optional choices with `taxid`, `run`, and `reason` |
| `accepted_samples.tsv` | Reviewed successes registered by `record` |

Adopted `metadata.tsv`, `selection.tsv`, and `provenance.json` live beside the
dataset settings. Work directories hold stage outputs, candidates, logs, and
`run.json` with inputs, status, and dates. Builds in `results/` preserve the
adopted provenance. Keep dataset definitions and reviewed tables in Git.

## Adding another dataset

Copy a dataset's configuration and rules; adjust its query, paths, tissue, BUSCO
lineage, and OrthoDB node. Match `sample_group` in `selection.yaml` to the
parameter in `rules/select_rules.tsv`. Start with an exclusion table containing
an `accession` header; metadata and acceptance records are created as needed.
Pass `--config datasets/<name>/selection.yaml` to metadata commands and the
corresponding `build.yaml` to build commands.
