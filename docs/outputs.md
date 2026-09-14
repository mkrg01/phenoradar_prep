# Outputs and TPM interpretation

[Documentation](index.md)

## Directory layout

The configured `run_name` selects `results/<run_name>/`, `work/<run_name>/`,
and `logs/<run_name>/`. Reusable references live in `resources/`.

| Shared resource | Fixed location |
| --- | --- |
| Taxonomy snapshot | `resources/taxonomy/taxa.sqlite` |
| OrthoDB v12 node | `resources/orthodb/v12_<node>/` |
| KOfam/KEGG snapshot and downloads | `resources/kegg/snapshot_v1/`, `resources/kegg/downloads/` |
| ASTRAL build | `resources/phylogeny_tools/` |
| TimeTree responses | `resources/timetree_cache/` |
| Launcher-managed cache | `.cache/` |

## Result files

```text
results/<run_name>/
  run.json                          # Resolved configuration, code hashes, Python environment
  metadata/
    metadata_all.tsv
    metadata_high_busco.tsv
    samples.tsv
    species_metadata.tsv            # species, trait, contrast_pair_id, family for PhenoRadar
    species_high_busco.txt
    selection.json                  # Selection counts and input/taxonomy hashes
    busco_completeness.svg
  proteins/                         # Species FASTA files and translation provenance
  orthogroups/
    mapping/
      manifests/                    # Chunk plan and FASTA manifests
      chunks/chunk_000/              # Annotations, hits, summary, provenance, native results
      gene_orthogroups.tsv           # Unique #query / ODB_OG pairs
      mappings.sqlite               # Indexed gene-to-OG mappings
      merge_qc.json
    expression/
      runs/                         # Per-run results and QC JSON
      tpm_sum.tsv                   # Sums of input TPM by orthogroup
      tpm.tsv                       # OG TPM rescaled to one million per run
      tpm_sum_wide.tsv
      tpm_wide.tsv
      mapping_qc.tsv
    alignments/                     # Optional all-copy OG protein alignments
  kegg/                             # Independent optional KO annotations/expression
  phylogeny/
    all/                            # All selected species
    phenotyped/                     # Species with an observed phenotype
    representatives/                # Trait-guided representative selection and inference
  phenoradar_inputs/                 # Selected completed inputs for PhenoRadar
  filtered/                         # Corresponding species-filtered result layout
```

Logs and ODB resource benchmarks are saved under
`logs/<run_name>/orthogroups/mapping/`. Temporary
ODB work is stored under `work/<run_name>/orthogroups/mapping/`.

Branch-specific inventories are in the [KEGG](kegg.md#outputs),
[alignment](alignments.md#outputs-and-phenoradar), [phylogeny](phylogeny.md#outputs),
[dating](dating.md#outputs), [taxonomic review](taxonomy_audit.md#outputs-and-figures),
and [contrast-pair](contrast_pairs.md#outputs) guides. The
[filtered export](species_filter.md#exported-dataset) has its own manifest;
[PhenoRadar collection](phenoradar_inputs.md#published-files) links selected inputs.

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

## Provenance and optional results

Retain `run.json`, branch-specific input checksums, QC, and execution records
with the reference snapshots. See [KO quantification](kegg.md#assignment-and-quantification)
for KO expression semantics and [OG alignments](alignments.md) for sequence outputs.

For results moved from an older directory structure, see the
[migration record guidance](migration.md#relocated-results).
