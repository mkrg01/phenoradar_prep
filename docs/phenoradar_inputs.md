# Use analysis outputs in PhenoRadar

[Documentation](index.md) · [Output formats](outputs.md)

Point PhenoRadar at the files in `results/BUILD/downstream/ANALYSIS/`.
Use files from the same analysis or from the same [exported subset](species_filter.md).

## Input files

Paths are relative to the analysis directory.

| Input | Path |
| --- | --- |
| Traits, taxonomy, and available contrast groups | `metadata/species_metadata.tsv` |
| Sample IDs, biological species, and runs | `metadata/samples.tsv` |
| Trait states, missing count, and source checksums | `metadata/species_trait.json` |
| OG expression | `orthogroups/expression/tpm.tsv` or `orthogroups/expression/tpm_sum.tsv` |
| KO expression | `kegg/ko_tpm_sum.tsv` |
| Trimmed OG alignments | `orthogroups/alignments/*.faa` |
| Original per-sequence TPM | `expression/sequence_tpm.tsv.gz` |
| Species tree | `phylogeny/<tree>/species_tree.nwk` |
| Dated tree, when enabled | `phylogeny/<tree>/dating/species_tree.dated.nwk` |

The `all` target creates the PhenoRadar metadata after the enabled contrast
branches finish. `species` identifies a sample in PhenoRadar inputs and in tree
tips. Use `metadata/samples.tsv` to recover biological species and run IDs.
Keep repeated samples of a species in the same validation group.

## Expression settings

```yaml
data:
  metadata_path: results/leaf/downstream/c4/metadata/species_metadata.tsv
  tpm_path: results/leaf/downstream/c4/orthogroups/expression/tpm.tsv
  species_col: species
  feature_col: orthogroup
  value_col: tpm
```

OG `tpm` is rescaled over genes assigned to OGs. To use the original gene TPM
sums, choose `tpm_sum.tsv` and set `value_col: tpm_sum`. For KO expression,
use `kegg/ko_tpm_sum.tsv`, `feature_col: ko`, and `value_col: tpm_sum`.
From v0.8.11, these KO values sum renormalized OG `tpm`, using all OG-assigned
sequences as the normalization basis. They are not renormalized again by KO.
Use matching [orthogroup annotations](https://github.com/mkrg01/phenoradar/blob/main/docs/data-format.md#orthogroup-annotations)
with OG features, or leave `data.orthogroup_annotation_path` unset.

Absent sample–feature coordinates are filled with zero by default in PhenoRadar.
With supported models, set `preprocess.expression.absent_feature_fill: nan` to
keep them missing. Use the long tables to preserve the distinction between an
absent coordinate and an observed zero. See [TPM interpretation](outputs.md#tpm-interpretation)
and the [PhenoRadar data format](https://github.com/mkrg01/phenoradar/blob/main/docs/data-format.md#expression-tsv).

## Amino acid inputs

```yaml
data:
  feature_mode: amino_acid
  metadata_path: results/leaf/downstream/c4/metadata/species_metadata.tsv
  alignment_dir: results/leaf/downstream/c4/orthogroups/alignments
  sequence_tpm_path: results/leaf/downstream/c4/expression/sequence_tpm.tsv.gz
preprocess:
  amino_acid:
    aggregation: tpm_weighted_frequency
```

TPM weighting joins on sample ID and exact FASTA sequence ID. Every imported
alignment record needs one matching TPM, including measured zeros. The
[PhenoRadar amino acid guide](https://github.com/mkrg01/phenoradar/blob/main/docs/amino-acid.md)
covers other aggregation modes and combining amino acid features with OG
expression.

## Validation groups

With one contrast branch, use `contrast_pair_id` for original pairs or
`expanded_pair_clade_id` for [expanded clades](contrast_pairs.md#expanded-pair-clades).
With several branches, use a suffixed column such as `contrast_pair_id_all` or
`expanded_pair_clade_id_phenotyped`. Unassigned samples have blank IDs.

For leave-one-group-out validation on expanded clades:

```yaml
data:
  contrast_pair_col: expanded_pair_clade_id
split:
  group_col: expanded_pair_clade_id
  outer_cv_strategy: logo
```

Using the same column in both settings sends labeled samples with blank groups to
`external_test` and makes any pair-aware ranking use the expanded groups. For
taxonomic holdouts, set `split.group_col` to `family` or `order` and review
missing group labels before running. See [PhenoRadar's split rules](https://github.com/mkrg01/phenoradar/blob/main/docs/data-format.md#metadata-tsv).

## Verify completed outputs

```bash
./run_analysis.sh status --analysis results/leaf/downstream/c4 --verify
```

Verification checks the recorded files and consistency of the completed
analysis. It can read large tables. Use plain `status` for job progress.
