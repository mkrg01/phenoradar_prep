# KEGG annotation and KO expression

[Documentation](index.md) · [Reference data](references.md#kofam-and-kegg-reference)

The KEGG branch annotates one median-length protein per observed OG with
KofamScan, then sums original OG TPM by KO. Representatives come from samples
selected for the analysis. Searches use KOfam's `eukaryote.hal`; KOs outside
that selection are not assessed.

## Run the branch

```bash
./run_analysis.sh submit --analysis results/leaf/downstream/carnivory --target kegg
```

To include KEGG in `all`, set `kegg.enabled: true` before preparing the analysis.
The explicit target works either way. Missing references are prepared automatically.
Annotation jobs search up to 5,000 OGs with 4 CPUs/8 GB by default; see
[resource overrides](running.md#resource-budgets).

## Representative selection

Choose the sequence closest to the median OG protein length. Ties favor the
longer sequence, then fewer nonstandard amino acids, then the smallest gene ID.
One terminal stop is removed; empty sequences, internal stops, and non-letter
characters are excluded. An OG without a valid sequence receives no KO.
Proteins assigned to multiple OGs are rejected.

These labels describe the selected representative. Other OG members are not
searched, so within-OG functional differences are not assessed. Donor IDs and
selection details are saved for review.

## Assignment and quantification

Only KofamScan hits marked as accepted contribute expression. `kegg.ambiguity`
controls OGs with several accepted KOs:

| Policy | Behavior |
| --- | --- |
| `duplicate` (default) | Assign full OG TPM to every accepted KO |
| `drop` | Exclude multi-KO OGs |
| `error` | Reject quantified multi-KO OGs, including observed zeros |

Values sum **original OG `tpm_sum`, without renormalization**. With `duplicate`,
an OG with TPM 10 assigned to two KOs contributes 10 to each, so total KO TPM
can exceed input TPM. Repeated hits to the same KO count once.

## Outputs

Under `results/<build>/downstream/<analysis>/kegg/`:

| Output | Contents |
| --- | --- |
| `ko_tpm_sum.tsv`, `ko_tpm_sum_wide.tsv` | Long/wide expression; long columns: `species`, `run`, `ko`, `tpm_sum` |
| `ko_support.tsv` | Observed contributing OG counts and TPM per KO/run |
| `mapping_qc.tsv` | Per-run counts and original, ODB-mapped, and KO-retained TPM |
| `orthogroups.tsv` | OG annotation status and representative metadata |
| `og_kos.tsv` | OG/KO hits and donor IDs; use `accepted=1` rows |
| `representatives/representatives.tsv` | Selection details, lengths, and membership counts |
| `annotation/<batch>/detail.tsv` | Raw hits, scores, thresholds, and acceptance markers |
| `ko_modules.tsv`, `ko_pathways.tsv` | KO group memberships |
| `provenance.json`, `annotation_provenance.json`, `reference_qc.json` | Method, completion, reference, and checksums |

Measured zeros remain zero. KOs with no quantified OGs have blank TPM in wide/
support tables and no numeric long-table row; partially quantified KOs sum only
observed OGs. All-zero runs are valid. Duplicate OG rows and negative/nonfinite
TPM are errors.

Coverage QC counts each retained OG once. `total_tpm` includes all quantified
transcripts, `odb_mapped_tpm` counts mapped OGs, and `retained_tpm` counts OGs
retained by KO assignment and the ambiguity policy. Contribution totals across
KOs can count an OG more than once.

Use the [collector](phenoradar_inputs.md#data-checks-and-downstream-use) to pass KO
features to PhenoRadar. A [post hoc subset](species_filter.md) keeps the original
OG labels and donor metadata. A new analysis selects new representatives from
its chosen samples.

## References

[KofamScan](https://github.com/takaram/kofam_scan) ·
[KofamKOALA paper](https://doi.org/10.1093/bioinformatics/btz859) ·
[KEGG API](https://www.kegg.jp/kegg/rest/keggapi.html)
