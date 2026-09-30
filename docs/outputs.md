# Outputs and TPM interpretation

[Documentation](index.md)

## Directory layout

Build `name` in `config/build.yaml` (or `submit --name`) and analysis
`submit --name` determine these directories:

| Path | Contents |
| --- | --- |
| `input/` | User-maintained metadata, traits, and other auxiliary inputs |
| `results/<build>/` | Frozen build settings/metadata, completion and job records |
| `results/<build>/database/` | Portable CDS/BUSCO/quant, proteins, mappings, per-sample OG expression, and provenance |
| `results/<build>/work/` | Staged inputs, GeneGalleon workspace, translation/mapping/expression computation, and downstream work |
| `results/<build>/downstream/<analysis>/` | Frozen conditions/inputs, job records, and downstream results |
| `results/<build>/logs/{database,downstream/<analysis>}/` | Rule logs and benchmarks |
| `resources/` | Shared reference and software caches |
| `results/<build>/work/cache/` | Internal stage receipts, translations, ODB mappings, and OG expression for this build |
| `results/<build>/reuse.json` | Source database checksums and the source chosen for each reused sample |

`reuse_from` in `config/build.yaml` selects one or more completed databases.
New builds stage matching products and compute missing ones inside their own
`results/<build>/work/`; completed products are published in `database/`.
`migrations/` is not a destination for new products. Existing historical files
remain in place. A completed database can be reused without its original working
caches, either by another build or by downstream analysis.

New CDS/BUSCO/quant files originate in the build's
`work/genegalleon/<sample>/output/transcriptome_assembly/` workspace. Staged input
directories contain staged links/copies; they are separate from top-level `input/`.
Copy the [completed products bundle](datasets.md#copying-a-completed-build-to-another-project)
for reuse elsewhere. Treat generated products as immutable. The database stores
`expression/runs/<run>.tsv` (both `tpm_sum` and normalized `tpm`) and
`expression/runs/<run>.qc.json`. Downstream selects and merges these saved values;
it does not recompute them. Proteins and mapping tables are staged only when a
requested downstream branch needs them.

The layout applies to newly prepared builds and downstream runs. Existing
`builds/`, `analyses/`, and old top-level `results/<run>/`, `work/`, and `logs/`
are not moved automatically. Frozen records and caches can reference those paths.
Old schema-3/4 `products/` bundles remain readable; downstream computes their
missing OG expression, rejecting genes mapped to multiple OGs. Resume old frozen
jobs with their original workflow checkout; code verification still applies.

## Result files

```text
results/<build>/downstream/<analysis>/
  run.json                          # Configuration and provenance
  metadata/                         # Selection, samples, traits, BUSCO QC
  proteins/                         # Proteins reused from the build
  orthogroups/
    mapping/
      snapshot.json                 # Selected species, checksums, mapping QC
      species/{species}.tsv.gz       # gene_id, orthogroup; blank OG means unmapped
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

Mapping tables retain all genes and all distinct OG assignments. Their snapshot
records the input protein hash and QC per species, so unchanged tables can be
reused without rebuilding a global database.

## TPM interpretation

Long tables contain `species` (analysis sample ID), `run`, `orthogroup`, and `tpm_sum` or `tpm`.
Wide tables have one row per run and one column per OG; missing combinations
are zero. Each row is a separately processed sample; multiple samples may belong to the same biological species.

- `tpm_sum`: original TPM summed by OG, excluding unmapped genes.
- `tpm`: retained OG values rescaled to one million per run, describing relative
  expression within the retained OG set.

Duplicate gene/OG pairs count once. If a quantified gene maps to multiple OGs,
aggregation always stops and reports the ambiguous genes. This behavior is not configurable.

Review `mapping_qc.tsv` for mapped/retained TPM fractions and ambiguous targets.
Duplicate targets, negative/nonfinite TPM, or no positive retained TPM stop
aggregation. [KO expression](kegg.md#outputs) uses different normalization and
missing-value rules.

The [collector](phenoradar_inputs.md) creates `phenoradar_inputs/tpm.tsv` with
exactly `species`, `orthogroup`, and `tpm`, preserving source values. Each sample ID must identify exactly one run; it never averages biological samples; run-level tables remain
available for QC.
