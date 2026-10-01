# PhenoRadar inputs

[Documentation](index.md) · [Output formats](outputs.md)

The collector writes results to
`results/<build>/downstream/<analysis>/phenoradar_inputs/`. It creates a
three-column OG TPM table and links other completed results, so keep the source
files available.

## Collecting results

Collection runs automatically after successful `all` execution. After running
additional branches, collect again:

```bash
./run_analysis.sh submit --analysis results/leaf/downstream/carnivory --target phenoradar_inputs
```

The collector starts no analyses or downloads. It requires `metadata/samples.tsv`
and `metadata/species_metadata.tsv`; expression is optional. Logs identify absent
or incomplete sections. Validation failures preserve the previous collection.

## Published files

| Path under `phenoradar_inputs/` | Contents |
| --- | --- |
| `species_metadata.tsv` | Base traits and taxonomy |
| `tpm.tsv` | OG expression: `species`, `orthogroup`, `tpm` |
| `alignments/{og}.faa` | All-copy OG alignments |
| `metadata/`, `proteins/`, `orthogroups/` | Original data and QC |
| `kegg/` | KO expression, annotation, support, and group maps |
| `phylogeny/{all,phenotyped,representatives}/` | Completed trees and related results |

Choose the tree explicitly, for example `phylogeny/all/species_tree.nwk` or its
dated counterpart. Pair IDs remain in each branch's `contrast/species_metadata.tsv`;
they are not merged into base metadata.

Existing headerless `orthogroup_annotations.tsv[.gz]` files (OG/taxid/description)
are collected from the result root, `orthogroups/`, or `orthogroups/mapping/`.
Descriptions are not downloaded.

## Data checks and downstream use

For OG expression in PhenoRadar:

```yaml
data:
  tpm_path: results/leaf/downstream/carnivory/phenoradar_inputs/tpm.tsv
  species_col: species
  feature_col: orthogroup
  value_col: tpm
```

The export preserves source [TPM values](outputs.md#tpm-interpretation).
`species` identifies a sample with exactly one run; samples are never averaged.
Use `metadata/samples.tsv` to group samples by biological species.

For KO expression, use `kegg/ko_tpm_sum.tsv`, `feature_col: ko`,
`value_col: tpm_sum`, and `orthogroup_annotation_path: null`. Review
[KO normalization and missing values](kegg.md#outputs) before use.

Normal analysis exclusions need no extra collection settings. For a
[post hoc subset](species_filter.md), filter first and collect with the same
override; collection uses the matching `filtered/` snapshot.

See PhenoRadar's [quick start](https://github.com/mkrg01/phenoradar/blob/main/docs/quickstart.md)
and [data formats](https://github.com/mkrg01/phenoradar/blob/main/docs/data-format.md).
