# PhenoRadar inputs

[Documentation](index.md)

`phenoradar_inputs` validates and links completed results into
`results/<run_name>/phenoradar_inputs/` for [PhenoRadar](https://github.com/mkrg01/phenoradar).
It does not start analyses. Keep linked source files available; failed validation
preserves the previous collection.

## Collecting results

Complete at least one expression branch, then collect:

```bash
./run_pipeline.sh --configfile config/mydata.yaml \
  --cores 1 --resources mem_gb=4 -- phenoradar_inputs
```

Collection needs `metadata/species_metadata.tsv`, normally created during
preparation; see [backfilling older results](migration.md#backfilling-phenoradar-metadata).
The collection log lists ready, absent, incomplete, and disabled branches.

## Selecting optional inputs

```yaml
phenoradar:
  orthogroups: auto
  kegg: auto
  alignments: auto
  kegg_groups: [module]
  orthogroup_annotations: null
  contrast: null
  tree: null
```

For the first three settings, `auto` includes completed branches, `true` requires
them, and `false` omits them. Corrupt completed results fail validation.
For KO-only data, set `orthogroups: false` and `alignments: false`.

- `kegg_groups`: `[module]`, `[pathway]`, both, or `[]` for expression alone.
- `orthogroup_annotations`: existing headerless OG/taxid/description TSV[.gz] from
  the same OrthoDB release. The collector checks format and expressed-OG overlap.
- `contrast`: one completed branch, `phylogeny/all/contrast`,
  `phylogeny/phenotyped/contrast`, or `phylogeny/representatives/contrast`.
  Pair IDs are joined to base metadata; unpaired species remain. Null leaves them blank.
- `tree`: an existing Newick path whose tips exactly match metadata species.
  Choose the molecular or dated tree explicitly. The collector does not prune trees.

## Published files

```text
results/<run_name>/phenoradar_inputs/
  species_metadata.tsv
  tpm.tsv                         # OG expression
  kegg/ko_tpm_sum.tsv              # KO expression
  kegg/ko_modules.tsv              # Requested KO grouping maps
  kegg/ko_pathways.tsv
  alignments/{og}.faa              # Selected completed alignments
  orthogroup_annotations.tsv[.gz]  # Optional descriptions
  species_tree.nwk                 # Explicitly selected tree
```

Metadata columns are `species`, the `contrast.trait` column (normally `C4`),
`contrast_pair_id`, and `family`. Traits come from `inputs.species_trait` and must
be `0`, `1`, or blank. Producer QC/provenance stay with source results.

## Expression and species checks

Collection requires **one run per species**; select one before collection if the
dataset contains multiple runs. Species/run identities must match the manifest,
with no duplicate species/feature coordinates and finite nonnegative values.
Every selected species needs expression rows. Missing KO expression is reported
as missing, not replaced with zeros; review [KO support](kegg.md#outputs).

OG `tpm.tsv` is normalized to one million per run. KO values sum original TPM and
can overlap across KOs. For KO expression, set these options in PhenoRadar:

```yaml
data:
  tpm_path: results/run001/phenoradar_inputs/kegg/ko_tpm_sum.tsv
  feature_col: ko
  value_col: tpm_sum
  orthogroup_annotation_path: null
```

Replace `run001` with your run name. Alignments retain all copies and columns;
see [gene IDs and sequence format](alignments.md#outputs-and-phenoradar).

## Species exclusions

With nonempty `exclude_species`, run `filter_species` first. Collection uses the
matching completed `filtered/` snapshot and validates its exclusions and sources;
it does not refresh that snapshot itself. Recomputed contrast branches can be
selected above. Set `tree` to the matching pruned tree when needed. Clearing the
exclusion list selects the original data on the next collection.

For downstream setup, see PhenoRadar's
[quick start](https://github.com/mkrg01/phenoradar/blob/main/docs/quickstart.md) and
[data formats](https://github.com/mkrg01/phenoradar/blob/main/docs/data-format.md).
