# KEGG annotation and KO expression

[Documentation](index.md)

The `kegg` target annotates translated proteins with KofamScan and sums original
input TPM by KEGG Orthology (KO), independently of ODB mapping.

## Run the branch

```yaml
kegg:
  enabled: true
  threads: 4
  mem_gb: 8
  ambiguity: duplicate
```

```bash
./run_pipeline.sh --configfile config/mydata.yaml \
  --cores 24 --resources mem_gb=128 -- kegg
```

The explicit target works without `enabled`; enabling it adds KEGG to `all`.
The container supplies KofamScan, and missing [references](references.md#kofam-and-kegg-reference)
are prepared automatically. Annotation runs once per species. Abundance or
ambiguity-policy changes reuse annotations and repeat aggregation.

## Assignment and quantification

One terminal stop is removed before HMMER input; internal stops are rejected.
KofamScan's acceptance marker determines assignments using unrounded scores.
Hits below threshold or without a defined threshold never receive TPM.

| Status | Meaning |
| --- | --- |
| `unique` | One accepted KO |
| `ambiguous` | Multiple accepted KOs |
| `threshold_missing` | No accepted KO; at least one hit lacks a threshold |
| `below_threshold` | Hits exist but none passes |
| `unannotated` | No reported hit |

`kegg.ambiguity` controls multi-KO genes:

- `duplicate` (default): full TPM contributes to every accepted KO. TPM 10 assigned
  to two KOs adds 10 to each; repeated hits to the same KO count only once.
- `drop`: exclude multi-KO genes from KO sums.
- `error`: reject quantified multi-KO genes, including observed zeros.

Values sum **original TPM without renormalization**. Under `duplicate`, features
overlap and the sum across KOs can exceed total input TPM. Multiple accepted KOs
do not by themselves distinguish multifunctionality from annotation ambiguity.

Use `accepted=1` rows in `gene_kos.tsv` for all assignments. `genes.tsv:selected_ko`
is populated only for unique assignments. QC counts each retained gene/its TPM
once; `quantified_assignments` and `ko_tpm_sum` include contributions across KOs.
Duplicate target IDs and negative/nonfinite TPM values are errors.

## Outputs

Under `results/<run_name>/kegg/`:

| Output | Contents |
| --- | --- |
| `species/<odb_species>/detail.tsv` | Raw scores, thresholds, E-values, and acceptance markers |
| `gene_kos.tsv`, `genes.tsv` | Merged hit assignments and all proteins, including unannotated genes |
| `ko_tpm_sum.tsv`, `ko_tpm_sum_wide.tsv` | Long and wide KO expression |
| `ko_support.tsv` | Annotated/quantified gene counts per KO/run |
| `mapping_qc.tsv`, `runs/<run>.qc.json` | Per-run coverage and TPM diagnostics |
| `ko_modules.tsv`, `ko_pathways.tsv` | KO group memberships |
| `reference_qc.json`, `species/<odb_species>/provenance.json` | Reference and annotation provenance |

Long expression columns are `species`, `run`, `ko`, and `tpm_sum`. Support adds
`annotated_genes` and `quantified_genes`. Measured zeros remain zero; a KO with no
quantified genes has blank TPM in support/wide tables and no numeric long-table
row. Partially quantified KOs sum observed genes only. Runs with no retained KOs
still produce QC; all-zero abundance is valid, with undefined fractions left blank.

These outputs do not establish genomic absence, pathway activity, or MODULE
completeness. Review support and QC before interpreting missing values.

## Passing KO features to PhenoRadar

The [input collector](phenoradar_inputs.md) links expression and requested
module/pathway maps. It requires one run per species; the guide gives the
PhenoRadar KO column settings.

## References

- [KofamScan documentation](https://github.com/takaram/kofam_scan)
- [KofamKOALA paper](https://doi.org/10.1093/bioinformatics/btz859)
- [KEGG API manual](https://www.kegg.jp/kegg/rest/keggapi.html)
- [KEGG MODULE definitions](https://www.kegg.jp/kegg/module.html)
