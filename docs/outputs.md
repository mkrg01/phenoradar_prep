# Outputs and TPM interpretation

[Documentation](index.md)

## Directory layout

`run_name` selects directories under `results/`, `work/` (retained work), and
`logs/` (per-step logs/benchmarks). Shared [references](references.md) live in
`resources/`.

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
    alignments/                     # Optional OG alignments
  kegg/                             # Optional KO results
  phylogeny/{all,phenotyped,representatives}/
  filtered/                         # Curated species subset
  phenoradar_inputs/                 # Downstream collection
```

Start with `metadata/selection.json`, `samples.tsv`, and
`busco_completeness.svg` for selection QC; `species_metadata.tsv` is base
PhenoRadar metadata. Optional outputs are described in their
[analysis guides](index.md#choose-an-analysis).

## TPM interpretation

Long tables have `species`, `run`, `orthogroup`, and `tpm_sum` or `tpm`.
Wide tables have one row per run and one column per OG; missing combinations
are zero. Multiple runs per species are kept separately.

- `tpm_sum`: original TPM summed by OG; unmapped genes excluded.
- `tpm`: retained OG values rescaled to one million per run, describing relative
  expression within the retained OG set.

Duplicate gene/OG pairs count once. `tpm.multimap` controls genes mapped to
multiple OGs:

| Policy | Behavior |
| --- | --- |
| `error` (default) | Stop and report ambiguous genes |
| `drop` | Exclude multi-OG genes |
| `split` | Divide TPM equally among assigned OGs |

Review `mapping_qc.tsv` for mapped/retained TPM fractions and ambiguous targets.
Duplicate target IDs, negative/nonfinite TPM, or no positive retained TPM stop
aggregation. [KO expression](kegg.md#outputs) has different normalization and
missing-value rules.

## Provenance

Keep `run.json`, branch QC/provenance, the Release's `image.json`, and reference
snapshots with the analysis.
