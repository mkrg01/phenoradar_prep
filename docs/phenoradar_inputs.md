# PhenoRadar inputs

[Documentation](index.md)

`phenoradar_inputs` collects completed results in
`results/<analysis>/phenoradar_inputs/`. It writes a three-column OG TPM table
and links other results. Keep linked source files available.

## Collecting results

`run_analysis.sh` collects automatically after successful `all` execution.
After additional branches finish, collect again:

```bash
./run_analysis.sh submit --analysis analyses/analysis001 --target phenoradar_inputs
```

The collector starts no analyses or downloads. It requires `metadata/samples.tsv`
and `metadata/species_metadata.tsv`; expression is optional. Logs report ready,
absent, and incomplete sections. Validation failures preserve the previous collection.

## Published files

| Path under `phenoradar_inputs/` | Contents |
| --- | --- |
| `species_metadata.tsv` | Base traits/taxonomy metadata |
| `tpm.tsv` | OG expression: `species`, `orthogroup`, `tpm` |
| `alignments/{og}.faa` | Completed all-copy OG alignments |
| `metadata/`, `proteins/`, `orthogroups/` | Original QC, proteins, mappings, expression, and alignments |
| `kegg/` | KO expression, annotation, support/QC, and group maps |
| `phylogeny/{all,phenotyped,representatives}/` | Completed trees, marker alignments, calibrations, pairs, and taxonomy reports |

Unfinished outputs and native solver work are omitted. Choose the downstream tree
explicitly, e.g. `phylogeny/all/species_tree.nwk` or its dated counterpart. Each
contrast branch keeps its own `contrast/species_metadata.tsv`; pair IDs are not
combined with base metadata or other branches.

Existing `orthogroup_annotations.tsv[.gz]` (headerless OG/taxid/description) files
are collected from the result root, `orthogroups/`, or `orthogroups/mapping/`;
descriptions are not downloaded.

## Data checks and downstream use

Collection validates identities, values, checksums, alignments, and tree tips.
The OG export requires one run per species and never averages runs. It removes
the `run` column while preserving source TPM values, already rescaled to one
million per run. See [TPM interpretation](outputs.md#tpm-interpretation).

For OG expression in PhenoRadar:

```yaml
data:
  tpm_path: results/analysis001/phenoradar_inputs/tpm.tsv
  species_col: species
  feature_col: orthogroup
  value_col: tpm
```

For KO expression, use `kegg/ko_tpm_sum.tsv`, `feature_col: ko`,
`value_col: tpm_sum`, and `orthogroup_annotation_path: null`. KO sums retain
original TPM and different [missing-value rules](kegg.md#outputs).
See PhenoRadar's [quick start](https://github.com/mkrg01/phenoradar/blob/main/docs/quickstart.md)
and [data formats](https://github.com/mkrg01/phenoradar/blob/main/docs/data-format.md)
for supported downstream inputs.

## Species exclusions

Normal `analysis.yaml` exclusions apply before analysis; collection needs no extra
filtering. For a [post hoc export](species_filter.md), run `filter_species` first
and collect with the same low-level override. Collection then requires the matching
`filtered/` snapshot and does not mix it with unfiltered results.
