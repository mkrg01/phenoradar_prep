# Use analysis outputs in PhenoRadar

[Documentation](index.md) · [Output formats](outputs.md)

Use completed files directly from
`results/<build>/downstream/<analysis>/`. No `phenoradar_inputs/` directory is
created. Choose an expression table, metadata, and any optional tree, alignment,
or annotation inputs explicitly.

## Input files

Paths below are relative to the analysis directory.

| Input | Path |
| --- | --- |
| Input metadata: traits, taxonomy, and available contrast pairs | `metadata/species_metadata.tsv` |
| Intermediate traits and taxonomy before contrast assignment | `metadata/species_metadata_base.tsv` |
| Sample IDs, biological species, and runs | `metadata/samples.tsv` |
| Original trait labels and numeric encoding | `metadata/trait_encoding.json` |
| OG expression | `orthogroups/expression/tpm.tsv` |
| Per-sequence expression | `expression/sequence_tpm.tsv.gz` (`species`, `sequence_id`, `tpm`) |
| KO expression | `kegg/ko_tpm_sum.tsv` |
| OG alignments | `orthogroups/alignments/{og}.faa` |
| Species tree | `phylogeny/{all,phenotyped,representatives}/species_tree.nwk` |
| Dated species tree, when enabled | `phylogeny/{all,phenotyped}/dating/species_tree.dated.nwk` |
| Branch-specific contrast pairs | `phylogeny/{tree}/contrast/species_metadata.tsv` |
| OG-to-KO evidence | `kegg/og_kos.tsv` (`accepted=1` rows) |
| KO group memberships | `kegg/ko_modules.tsv`, `kegg/ko_pathways.tsv` |

`all` creates `metadata/species_metadata.tsv` after the configured contrast
branches finish. It includes `species`, the trait, `contrast_pair_id`, and the
available `kingdom`, `phylum`, `class`, `order`, `family`, and `genus` names.
Unrecorded ranks remain blank. `species_metadata.json` records the source files
and checksums.

With one configured contrast branch, its pair IDs populate `contrast_pair_id`.
Unpaired samples remain in the table with blank IDs. With multiple contrast
branches, IDs are kept in separate columns such as `contrast_pair_id_all` and
`contrast_pair_id_representatives`; `contrast_pair_id` stays blank. Choose the
desired column using PhenoRadar's `data.contrast_pair_col`. With contrasts
disabled, pair IDs remain blank; saved files from earlier runs are ignored.

Optional files exist only after their branch has completed. QC, representative
selection, and provenance stay beside the results; they are not extra expression
inputs. Choose the tree and pair branch deliberately. Pair IDs are local to each
branch. The final input metadata combines their assignments with taxonomy;
manual joins are unnecessary. Base metadata remains available for targets run
before contrast analysis. Use `family` or `order` as a grouping column when not
using contrast pairs.

## Expression settings

For OG expression:

```yaml
data:
  metadata_path: results/leaf/downstream/carnivory/metadata/species_metadata.tsv
  tpm_path: results/leaf/downstream/carnivory/orthogroups/expression/tpm.tsv
  species_col: species
  feature_col: orthogroup
  value_col: tpm
split:
  group_col: family
```

For KO expression, replace `tpm_path` with
`results/leaf/downstream/carnivory/kegg/ko_tpm_sum.tsv`, set `feature_col: ko`,
`value_col: tpm_sum`, and `orthogroup_annotation_path: null`.

Expression tables omit the redundant `run` column; the sample-to-run mapping
remains in `metadata/samples.tsv`. Older tables with `run` are also accepted by
PhenoRadar, which selects the configured columns.
`species` is a computational sample ID, not a biological-species grouping. Give
each independently assembled run a distinct sample ID; samples are never
averaged. Use `metadata/samples.tsv` to identify biological species.

OG `tpm` is normalized over retained OGs. KO `tpm_sum` sums original OG TPM
without renormalization. Review [TPM interpretation](outputs.md#tpm-interpretation)
and [KO ambiguity and missing values](kegg.md#assignment-and-quantification).

For future expression-weighted amino-acid features, the per-sequence table can
be joined to OG alignments by sample and exact FASTA sequence ID. It contains
original gene TPM, including zero values and unmapped genes; it does not contain
precomputed site features. Current PhenoRadar does not yet implement this MSA
feature input. See [per-sequence TPM](outputs.md#tpm-interpretation) for the contract.

Trait columns contain `0`, `1`, or missing values. Text states are
[encoded once](inputs.md#traits); their source labels remain in
`metadata/trait_encoding.json`.

Existing OrthoDB annotation TSVs can be passed separately using PhenoRadar's
`orthogroup_annotation_path`. Use the file matching the OG reference version;
annotation descriptions are not downloaded or combined automatically. KEGG maps
are separate evidence tables; supplying them as generalized annotation inputs
requires corresponding support in PhenoRadar.

## Completion and verification

`all` validates native outputs before recording completion. For an explicit
integrity check, use the existing command:

```bash
./run_analysis.sh status --analysis results/leaf/downstream/carnivory --verify
```

This checks recorded output hashes and consistency of sample IDs, expression,
completed OG alignments, trees, and available annotations. Plain `status` remains
a lightweight observation and does not read large output tables.

For a [post hoc subset](species_filter.md), point all relevant inputs at
`filtered/` under the same analysis. Do not mix filtered expression or metadata
with the original tree; use the saved pruned tree or matching contrast branch.

The former `phenoradar_inputs` target and collection directory generation have
been removed. Existing directories and completion records are retained and are
not automatically migrated or deleted. Existing collections are not refreshed.
When rerunning or extending an analysis, update PhenoRadar configurations to
native paths.

See PhenoRadar's [quick start](https://github.com/mkrg01/phenoradar/blob/main/docs/quickstart.md)
and [data formats](https://github.com/mkrg01/phenoradar/blob/main/docs/data-format.md).
