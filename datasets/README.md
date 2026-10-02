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
if no candidate remains. Retained runs are matched by accession in the newly
curated metadata. Their reviewed `taxid`, `bioproject`, and `scientific_name`
remain fixed, keeping `species_id` and `analysis_sample_id` stable. Other run
attributes come from the newly retrieved row. `selection.tsv` records the
source taxid, BioProject, and scientific name and flags identity or curation
changes for review. Explicit accession exclusions take precedence over retention.
To adopt an identity correction, review the source values, update the corresponding
identity fields in `accepted_samples.tsv`, and repeat `update` before preparing a
new build. Existing builds retain their frozen identities.

`accepted_samples.tsv` is authoritative even when local `metadata.tsv` is absent
or outdated. If an accepted run is absent from the new source, a saved local row
for the same run can be retained and flagged for review. Without either row,
`update` saves a candidate and reports `review_required` without publishing it.
The accession's taxid remains reserved, so another run cannot silently replace it.
Restore its metadata from the corresponding build's `source_metadata.tsv` or a
saved candidate, or explicitly exclude/override the run, then repeat `update`.
`accept` cannot publish a candidate with unresolved accepted runs.

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
| `select_rules.tsv` | AMALGKIT curation rules |
| `excluded_accessions.tsv` | Manual exclusions and recorded BUSCO failures; `accession` required |
| `overrides.tsv` | Optional choices with `taxid`, `run`, and `reason` |
| `accepted_samples.tsv` | Reviewed successes registered by `record`; accession is the `run` column |
| `selection.tsv` | Representative choices, source identity changes, and review reasons |
| `provenance.json` | Adopted metadata inputs, checksums, and dates |
| `accepted_samples.provenance.json` | QC registration inputs and decisions |

To replace curation rules with the pinned image's rule set, run
`./run_metadata.sh rules --replace`.

Adopted `metadata.tsv`, `selection.tsv`, and `provenance.json` live beside the
dataset settings. Git tracks the documented definitions, decision tables, and
provenance, together with the READMEs and shared environment recipe. Full
`metadata.tsv` and other dataset files are ignored, including acquisition data,
locks, and temporary files. Metadata can be rebuilt from new curated data and
the reviewed tables; a new query does not reproduce every historical attribute.

Work directories hold stage outputs, candidates, logs, and `run.json` with
inputs, status, and dates. Builds in `results/` preserve the adopted provenance,
`database/source_metadata.tsv` before build exclusions, and
`database/metadata.tsv` after exclusions. Preserve these bundles to recover an
exact historical metadata snapshot.

## Adding another dataset

Copy a dataset's configuration and rules; adjust its query, paths, tissue, BUSCO
lineage, and OrthoDB node. Set `selection.yaml`'s `sample_group` to a group allowed
by `select_rules.tsv`. Start with an exclusion table containing an `accession` header.
Pass `--config datasets/<name>/selection.yaml` to metadata commands and the
corresponding `build.yaml` to build commands.
