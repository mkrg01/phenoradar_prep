# Outputs and TPM interpretation

[Documentation](index.md)

## Directory layout

The configured `run_name` selects `results/<run_name>/`, `work/<run_name>/`,
and `logs/<run_name>/`. Reusable references live in `resources/`.
Slurm writes standard output to `pipeline-<job_id>.out` and standard error to
`pipeline-<job_id>.err` in the repository root.

## Result files

```text
results/<run_name>/
  run.json                          # Configuration and provenance
  metadata/                         # Selection, samples, traits, BUSCO QC
  proteins/                         # Translated species FASTAs
  orthogroups/
    mapping/gene_orthogroups.tsv     # Gene-to-OG assignments
    mapping/mappings.sqlite         # Indexed mappings
    expression/
      tpm_sum.tsv                   # Original TPM sums by OG
      tpm.tsv                       # Rescaled OG TPM
      tpm_sum_wide.tsv
      tpm_wide.tsv
      mapping_qc.tsv
    alignments/                     # Optional all-copy OG alignments
  kegg/                             # Optional KO annotations/expression
  phylogeny/{all,phenotyped,representatives}/
  filtered/                         # Curated species subset
  phenoradar_inputs/                 # Automatically collected downstream inputs/candidates
```

Selection results include `metadata_all.tsv`, `metadata_high_busco.tsv`,
`samples.tsv`, `species_high_busco.txt`, `selection.json`, and
`busco_completeness.svg`. `species_metadata.tsv` provides PhenoRadar metadata.
Per-step logs/benchmarks are under `logs/<run_name>/`; temporary files are under
`work/<run_name>/`. Reference locations are listed in [references](references.md).

See branch guides for [KO expression](kegg.md#outputs),
[alignments](alignments.md#outputs-and-phenoradar), [phylogeny](phylogeny.md#outputs),
[dating](dating.md#outputs), [taxonomy review](taxonomy_check.md#outputs-and-figures),
[contrast pairs](contrast_pairs.md#outputs), [filtered data](species_filter.md#exported-dataset),
and [PhenoRadar inputs](phenoradar_inputs.md#published-files).

## TPM interpretation

Long tables contain `species`, `run`, `orthogroup`, and either `tpm_sum` or `tpm`.
Wide tables contain one row per run, identified by `species` and `run`, with one
column per orthogroup. Missing run/orthogroup combinations are filled with zero.
Multiple runs from a species are not pooled or averaged.

- `tpm_sum` sums the original input TPM values assigned to each orthogroup,
  without rescaling. Unmapped genes are excluded.
- `tpm` rescales retained orthogroup values to sum to one million within each run.
  It therefore describes relative expression within the retained OG set.

Duplicate gene/OG pairs are removed before aggregation. If a gene maps to multiple
orthogroups, `tpm.multimap` controls its treatment:

| Policy | Behavior |
| --- | --- |
| `error` (default) | Stop and report the ambiguous genes |
| `drop` | Exclude genes assigned to multiple OGs |
| `split` | Divide each gene's TPM equally among its assigned OGs |

`mapping_qc.tsv` reports the fraction of input TPM mapped to OGs, ambiguous target
counts, the retained TPM fraction, and other mapping statistics. Aggregation
rejects duplicate target IDs, negative or nonfinite TPM values, and runs with no
positive TPM retained after mapping and ambiguity handling.

## Provenance

Keep `run.json`, branch QC/provenance, the Release's `image.json`, and reference
snapshots with the analysis. See [migration](migration.md#relocated-results) for
relocated results.
