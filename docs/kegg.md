# KEGG annotation and KO expression

[Documentation](index.md)

The `kegg` target selects **one median-length protein per observed OrthoDB
orthogroup (OG)**, annotates it with KofamScan, and assigns its accepted KO(s) to
that OG. KO expression sums the original TPM of the corresponding OGs.
Proteins and ODB mappings come from the completed build; only samples selected
for the analysis contribute representatives.

Searches use KOfam's `eukaryote.hal` selection. KOs outside this selection are
not assessed. These are representative-derived OG labels: individual members
are not searched, and within-OG functional differences are not assessed.

## Run the branch

```bash
./run_analysis.sh submit --analysis results/angiosperm_leaf_20260925/downstream/c4_photosynthesis_20260929 --target kegg
```

Set `kegg.enabled: true` before [preparing analysis](datasets.md#run-an-analysis)
to include KEGG in `all`; the explicit target works either way. Missing
[references](references.md#kofam-and-kegg-reference) are prepared automatically.
The frozen reference must contain `eukaryote.hal`; older snapshots without it
must be rebuilt from matching KOfam inputs.

Representatives are searched in batches of at most 5,000 OGs. Each annotation
job defaults to 4 CPUs/8 GB; representative selection and expression aggregation
each use 1 CPU/8 GB. See [resources](running.md#resource-budgets).
Completed annotation batches with matching inputs are reused on retry within
the same analysis. New analysis IDs do not share annotations automatically.

## Representative selection

Within each OG, calculate the median protein length and select the actual
sequence closest to it. Ties favor the longer sequence, then fewer nonstandard
amino acids, then the lexicographically smallest gene ID. For an even number of
members with different middle lengths, this selects the upper middle length.
A protein is counted once even if multiple expression runs share it; distinct
protein IDs are separate members.

One terminal stop is removed before measuring length. Empty sequences, internal
stops, and other non-letter characters are excluded. An OG with no valid sequence
is recorded as `no_valid_representative` and receives no KO. Proteins without an
ODB assignment do not enter this branch. Multiple ODB assignments for one protein
are rejected.

Only the selected representative is annotated. There are no additional searches
for unannotated OGs or OGs involved in major results, and no annotation of other
members for accuracy sampling. Representative IDs, lengths, membership counts,
and input checksums are retained for provenance.

## Assignment and quantification

KofamScan's acceptance marker determines assignments using unrounded scores.
Hits below threshold or lacking a threshold receive no TPM. An OG whose
representative receives no accepted KO remains unannotated.

`kegg.ambiguity` controls OGs with multiple accepted KOs:

| Policy | Behavior |
| --- | --- |
| `duplicate` (default) | Full OG TPM to every accepted KO; repeated hits to one KO count once |
| `drop` | Exclude multi-KO OGs |
| `error` | Reject quantified multi-KO OGs, including observed zeros |

Values sum **original OG `tpm_sum` without renormalization** from
`orthogroups/expression/runs/<run>.tsv`, reusing saved build expression when
available. The normalized OG `tpm` column is not used. No per-member KO expansion
is needed. With `duplicate`, OG TPM 10 assigned to two KOs adds 10 to each; the
sum across KOs can exceed input TPM.

## Outputs

Under `results/<build>/downstream/<analysis>/kegg/`:

| Output | Contents |
| --- | --- |
| `representatives/representatives.tsv` | Selected protein ID/sample, length, median, member counts, validity, sequence checksum, and batch per OG |
| `representatives/batches/*.faa` | One sequence per selected OG; FASTA identifiers are OG IDs |
| `representatives/provenance.json` | Selection method, ODB identity, input checksums, and batch inventory |
| `annotation/<batch>/detail.tsv` | Raw representative hits, scores, thresholds, E-values, and acceptance markers |
| `annotation/<batch>/{genes,gene_kos}.tsv` | Internal KofamScan tables; `gene_id` contains the OG query ID |
| `orthogroups.tsv` | Representative metadata and OG annotation status; `selected_ko` only for unique assignments |
| `og_kos.tsv` | OG/KO hits with the actual representative protein ID; use `accepted=1` rows |
| `annotation_provenance.json` | OG annotation completion and source batch receipts |
| `ko_tpm_sum.tsv`, `ko_tpm_sum_wide.tsv` | Long/wide KO expression |
| `ko_support.tsv` | Observed contributing OG counts (`quantified_orthogroups`) and TPM per KO/run |
| `mapping_qc.tsv` | Per-run OG counts and original/ODB-mapped/KO-retained TPM |
| `provenance.json` | Expression aggregation method, inputs, and output checksums |
| `ko_modules.tsv`, `ko_pathways.tsv` | KO group memberships |
| `reference_qc.json`, `annotation/<batch>/provenance.json` | Reference/search provenance |

Long expression columns remain `species`, `run`, `ko`, and `tpm_sum`.
Measured zeros remain zero. KOs without quantified OGs have blank TPM in
support/wide tables and no numeric long-table row; partially quantified KOs sum
observed OGs only. All-zero runs are valid; undefined QC fractions stay blank.
Duplicate OG rows and negative/nonfinite TPM are errors.

Coverage QC counts each retained OG/its TPM once; `quantified_assignments` and
`ko_tpm_sum` count contributions across KOs. `total_tpm` includes all quantified
input transcripts, `odb_mapped_tpm` counts mapped OG expression, and
`retained_tpm` counts expression retained by KO assignment and the ambiguity policy.

Use the [input collector](phenoradar_inputs.md#data-checks-and-downstream-use) to
pass KO features to PhenoRadar. A [post hoc species export](species_filter.md)
subsets expression while preserving the original shared OG labels and donor
metadata, even if a donor sample is excluded. It does not select new
representatives. Preparing a new analysis with a different sample selection
selects representatives from that new set.

## References

[KofamScan](https://github.com/takaram/kofam_scan) ·
[KofamKOALA paper](https://doi.org/10.1093/bioinformatics/btz859) ·
[KEGG API](https://www.kegg.jp/kegg/rest/keggapi.html) ·
[KEGG MODULE](https://www.kegg.jp/kegg/module.html)
