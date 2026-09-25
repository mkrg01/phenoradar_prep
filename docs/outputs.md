# Outputs and TPM interpretation

[Documentation](index.md)

## Directory layout

Names supplied to `prepare --name` determine these directories:

| Path | Contents |
| --- | --- |
| `input/` | User-maintained metadata, traits, and other auxiliary inputs |
| `builds/<build>/` | Frozen build settings/metadata, GeneGalleon workspace, and job records |
| `builds/<build>/products/` | Completed portable CDS/BUSCO/quant, proteins, mappings, and provenance |
| `results/build_<build>/` | Build translation and mapping outputs |
| `analyses/<analysis>/` | Frozen analysis settings/inputs and job records |
| `results/<analysis>/` | Downstream results |
| `work/<run>/`, `logs/<run>/` | Retained work and rule logs; run is `build_<build>` or `<analysis>` |
| `resources/` | Shared reference/software caches and species-product registry |

New CDS/BUSCO/quant files originate in the build's
`genegalleon/output/transcriptome_assembly/` workspace. Phase-local `input/`
directories contain staged links/copies; they are separate from top-level `input/`.
Copy the [completed products bundle](datasets.md#copying-a-completed-build-to-another-project)
for reuse elsewhere. Treat generated products as immutable.

## Result files

```text
results/<analysis>/
  run.json                          # Configuration and provenance
  metadata/                         # Selection, samples, traits, BUSCO QC
  proteins/                         # Proteins reused from the build
  orthogroups/
    mapping/                        # Selected-species gene-to-OG mappings
    expression/
      tpm_sum.tsv                   # Original TPM sums by OG
      tpm.tsv                       # Rescaled OG TPM
      tpm_sum_wide.tsv
      tpm_wide.tsv
      mapping_qc.tsv
    alignments/                     # Optional all-copy OG alignments
  kegg/                             # Optional KO results
  phylogeny/{all,phenotyped,representatives}/
  filtered/                         # Optional post hoc export
  phenoradar_inputs/                 # Collected downstream inputs
```

Start with `metadata/selection.json`, `samples.tsv`, and
`busco_completeness.svg` for selection QC. `species_metadata.tsv` is base
PhenoRadar metadata. See the [analysis guides](index.md#choose-an-analysis)
for optional outputs. Keep `run.json`, branch provenance, reference snapshots,
and the release's `image.json` with the analysis.

## TPM interpretation

Long tables contain `species`, `run`, `orthogroup`, and `tpm_sum` or `tpm`.
Wide tables have one row per run and one column per OG; missing combinations
are zero. Build currently requires one run per species.

- `tpm_sum`: original TPM summed by OG, excluding unmapped genes.
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
Duplicate targets, negative/nonfinite TPM, or no positive retained TPM stop
aggregation. [KO expression](kegg.md#outputs) uses different normalization and
missing-value rules.

The [collector](phenoradar_inputs.md) creates `phenoradar_inputs/tpm.tsv` with
exactly `species`, `orthogroup`, and `tpm`, preserving source values. It rejects
multiple runs per species rather than averaging them; run-level tables remain
available for QC.
