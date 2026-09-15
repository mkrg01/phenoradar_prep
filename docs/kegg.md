# KEGG annotation and KO expression

[Documentation](index.md)

The `kegg` target annotates proteins with KofamScan and sums original TPM by
KEGG Orthology (KO), independently of ODB mapping.

## Run the branch

```bash
./run_pipeline.sh --cores 24 --resources mem_gb=128 -- kegg
```

Set `kegg.enabled: true` to include KEGG in `all`; the explicit target works
either way. Missing [references](references.md#kofam-and-kegg-reference) are
prepared automatically. Annotation runs once per species and is reused after
abundance or ambiguity-policy changes. See [resources](running.md#resource-budgets).

## Assignment and quantification

One terminal stop is removed before HMMER; internal stops are rejected.
KofamScan's acceptance marker determines assignments using unrounded scores.
Hits below threshold or lacking a threshold receive no TPM.

`kegg.ambiguity` controls genes with multiple accepted KOs:

| Policy | Behavior |
| --- | --- |
| `duplicate` (default) | Full TPM to every accepted KO; repeated hits to one KO count once |
| `drop` | Exclude multi-KO genes |
| `error` | Reject quantified multi-KO genes, including observed zeros |

Values sum **original TPM without renormalization**. With `duplicate`, TPM 10
assigned to two KOs adds 10 to each; the sum across KOs can exceed input TPM.
Multiple assignments alone do not distinguish multifunctionality from ambiguity.

## Outputs

Under `results/<run_name>/kegg/`:

| Output | Contents |
| --- | --- |
| `species/<odb_species>/detail.tsv` | Raw hits, scores, thresholds, E-values, and acceptance markers |
| `gene_kos.tsv` | All hit assignments; use `accepted=1` rows |
| `genes.tsv` | All proteins and annotation status; `selected_ko` only for unique assignments |
| `ko_tpm_sum.tsv`, `ko_tpm_sum_wide.tsv` | Long/wide expression |
| `ko_support.tsv` | Annotated/quantified gene counts per KO/run |
| `mapping_qc.tsv`, `runs/<run>.qc.json` | Coverage and TPM diagnostics |
| `ko_modules.tsv`, `ko_pathways.tsv` | KO group memberships |
| `reference_qc.json`, `species/<odb_species>/provenance.json` | Reference/annotation provenance |

Long expression columns are `species`, `run`, `ko`, and `tpm_sum`.
Measured zeros remain zero. KOs without quantified genes have blank TPM in
support/wide tables and no numeric long-table row; partially quantified KOs sum
observed genes only. All-zero runs are valid; undefined QC fractions stay blank.
Duplicate target IDs and negative/nonfinite TPM are errors.

Coverage QC counts each retained gene/its TPM once; `quantified_assignments` and
`ko_tpm_sum` count contributions across KOs. Review support before interpreting
missing values: these results do not establish genomic absence, pathway activity,
or MODULE completeness.

## Passing KO features to PhenoRadar

The [input collector](phenoradar_inputs.md#data-checks-and-downstream-use) links
completed KO results and documents downstream column settings. Expression runs
remain separate.

## References

[KofamScan](https://github.com/takaram/kofam_scan) ·
[KofamKOALA paper](https://doi.org/10.1093/bioinformatics/btz859) ·
[KEGG API](https://www.kegg.jp/kegg/rest/keggapi.html) ·
[KEGG MODULE](https://www.kegg.jp/kegg/module.html)
