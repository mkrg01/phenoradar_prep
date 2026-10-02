# Curated datasets

This directory holds dataset definitions and reviewed metadata. Acquisition and
curation use AMALGKIT inside the GeneGalleon image pinned by each dataset's
`build.yaml`. Representative selection uses the shared
[metadata preparation script](../workflow/scripts/metadata_catalog.py).
Metadata preparation and database submission are separate commands.

| Dataset | Selection |
| --- | --- |
| [angiosperm_leaf](angiosperm_leaf/README.md) | One leaf RNA-seq run per original NCBI taxid |

## Choosing where to restart

Choose the entry point according to what needs to change:

| Purpose | Command | Work performed |
| --- | --- | --- |
| Refresh NCBI metadata | `./run_metadata.sh update --work <new-attempt>` | Run AMALGKIT metadata, curate, select representatives, and update dataset metadata |
| Reselect from saved candidates | `./run_metadata.sh update --work <new-attempt> --metadata <curated.tsv>` | Reuse the AMALGKIT-curated source and select representatives using current exclusions, overrides, and adoption records |
| Retry processing with the same samples | `./run_build.sh submit --build results/<existing-build>` | Keep the build's frozen metadata and settings; schedule failed or otherwise unfinished sample stages and missing downstream work |

Both metadata update paths follow the same retention policy: adopted runs stay
selected unless excluded or overridden. The second path needs the saved table
**before representative selection**, not the final dataset `metadata.tsv`, so
that alternative runs remain available. If tissue rules changed, first rerun
`curate --work <new-attempt> --metadata <saved-raw.tsv>`, then pass its curated
table to `update --metadata`; this also avoids another NCBI query.

After either metadata update, prepare a new build name and use `reuse_from: auto`
to retain verified stages for unchanged samples. A processing retry uses the
existing build: completed stages are skipped and missing prerequisites run as
needed. Resource changes can be supplied with `--resources <retry.yaml>`; sample
selection and scientific settings stay frozen.

For new BUSCO failures, run `record --build results/<build>` to add the
below-threshold accessions, then use the reselection path and prepare a new build.
Changes to dataset metadata or exclusions do not alter an already prepared build.

## Layout

```text
datasets/<name>/
  README.md                  # Scope, policy, and update instructions
  selection.yaml             # Search, tissue, selection paths, adoption BUSCO threshold
  build.yaml                 # Database settings and pinned GeneGalleon image
  excluded_accessions.tsv    # Manual exclusions and recorded BUSCO failures
  overrides.tsv              # Optional representative choices with reasons
  rules/
    select_rules.tsv          # Effective rules; edit these before curation
    upstream/select_rules.tsv # Unmodified snapshot from the pinned image
```

The initial directory contains only these inputs and any notes documenting their
origins. `update` creates `metadata.tsv`, `selection.tsv`, and `provenance.json`.
It also writes `accepted_samples.tsv`, initially without success evidence;
`record` adopts reviewed samples after their assembly, BUSCO, and quantification
finish and their BUSCO completeness meets the adoption threshold, without
waiting for the whole database. Use `record --dry-run` to inspect
completion, BUSCO counts, and proposed exclusions before recording QC.
Absent metadata and success tables are treated as an empty initial history.
Exporting or refreshing upstream rules creates `rules/upstream/select_rules.json`.

Keep these definitions, effective rules, exclusions, and reviewed metadata in Git.
Large source metadata, candidates, logs, and provenance for individual attempts
go under `work/datasets/<name>/<attempt>/`, which is ignored by Git. The SIF lives
in the existing `resources/software/` cache. Completed databases live under
`results/<build_name>/database/`, also ignored. No publishing step is configured.

## Selection policy

Each positive original `taxid` gets at most one run. AMALGKIT curation must mark
the row with the requested `sample_group` and `exclusion=no`.

1. Exclude exactly the listed accessions. Other runs from the same taxid and
   BioProject remain eligible; list each run when excluding a whole group.
2. Retain a previously adopted run, even if another run has more `total_bases`.
   Explicit exclusions take precedence. BUSCO counts are recorded QC evidence,
   not an additional retention filter.
3. For taxids without retained success, choose the eligible run with the largest
   positive `total_bases`; equal values resolve by ascending run accession.
4. Apply eligible choices in `overrides.tsv`, keeping their recorded reasons.
   The setting may be omitted or null; an absent or empty file means no overrides.

`accepted_samples.tsv` records adopted representatives; updating metadata does
**not** mark new runs as adopted. Run `record` after reviewing their QC. Completeness
is `(single + duplicated) / total`, and the threshold defaults to `0.5` (inclusive).
New samples with verified assembly and BUSCO products below this threshold are
automatically added to `excluded_accessions.tsv`, even if quantification is still
incomplete. Their reasons include the threshold, and their taxid, BioProject, and
source build are recorded. Existing manual decisions and custom columns are kept;
repeating `record` does not duplicate exclusions. `--runs` limits both adoption
and automatic exclusions; `--dry-run` changes neither table.

Missing or failed assembly/BUSCO stages do not trigger automatic exclusions.
Samples passing BUSCO still need completed quantification before adoption.
Add other unusable runs to `excluded_accessions.tsv` manually. Explicit exclusions
always prevail. `record` leaves metadata unchanged; run `update` afterward to
select the next eligible representative for each excluded taxid, or remove that
taxid if no candidate remains. A saved curated source can be reused without
another NCBI query. Prepare a new build to assemble replacement runs and reuse
other samples' completed stages.

A retained run missing from the latest source, or newly excluded by AMALGKIT,
remains selected with an informational review flag in `selection.tsv`. Use an
explicit exclusion or override to change an adopted choice. `busco_threshold`
applies only to new adoption through `record` or a `baseline` import; it is not
reapplied to already adopted representatives. Even raising the threshold does
not revoke their adoption. Analysis-stage BUSCO thresholds remain separate.

For Methods: metadata are retrieved and curated using AMALGKIT in a pinned
container. One representative run is selected per original NCBI taxid. Explicit
run exclusions take precedence, adopted representatives are retained, and
otherwise the eligible run with the greatest reported total bases is selected
(accession order breaks ties). Adoption requires successful assembly, BUSCO, and
quantification, QC review, and complete BUSCO fraction at or above the configured
threshold (default 50%). New runs with completed BUSCO results below the threshold
are recorded as accession exclusions before the next metadata update; execution
failures remain eligible for retry. Document any manual overrides and QC decisions.

## Adding another dataset

Create a sibling directory, copy the configuration structure, and change its
name, search query, `sample_group`, rules, and all dataset paths. Choose the
appropriate BUSCO lineage and OrthoDB node in `build.yaml`. Create
`excluded_accessions.tsv` with an `accession` header and optional `reason`,
`taxid`, `bioproject`, and `source_build` columns. If overrides are useful, create
`overrides.tsv` with `taxid`, `run`, and `reason` headers. Empty tables are valid.
Export and edit the rules; representative metadata
and success evidence can remain absent until the first update.

Run the same commands with `--config datasets/<name>/selection.yaml`; there is
no dataset-specific Python script to duplicate. `rule_set` chooses the packaged
AMALGKIT rule set when exporting upstream rules. Set its `sample_group` parameter
in the effective rules to match your intended tissue(s), and select a single
target group in `selection.yaml`.

`update --work work/datasets/<name>/<attempt>` fetches, curates, selects, and
publishes metadata in one command. `--dry-run` writes only a candidate;
`--metadata <curated.tsv>` reuses an already curated download after exclusions
or overrides change. The individual `fetch`, `curate`, `select`, and `accept`
commands remain available for inspection and rule development. See the
[angiosperm leaf instructions](angiosperm_leaf/README.md#updating-metadata)
for environment setup and the individual update commands.
The lightweight host environment is defined in `datasets/environment.yaml`;
it is independent of the PhenoRadar analysis container's environment recipes.

## Paths and builds

Run commands from the repository root. Paths in both YAML files resolve against
that root, and metadata can stay in `datasets/`; copying it into `input/` is not
necessary. `run_build.sh prepare` freezes a particular build under `results/`.
Subsequent metadata edits affect future builds, while the prepared build keeps
its saved inputs. Use a new build name for changed inputs or settings.
