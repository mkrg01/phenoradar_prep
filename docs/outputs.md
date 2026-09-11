# Outputs and TPM interpretation

[Back to README](../README.md)

Output paths use the `analysis` name from your [configuration](configuration.md).

## Output files

```text
results/<analysis>/
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
`logs/<analysis>/orthogroups/mapping/`. Temporary
ODB work is stored under `work/<analysis>/orthogroups/mapping/`.

The explicit [PhenoRadar input collection](phenoradar_inputs.md) publishes
`results/<analysis>/phenoradar_inputs/`. It links selected completed outputs and
uses base metadata prepared by the metadata step. Optional contrast pair IDs
are left-joined without removing unpaired or unannotated species. KEGG inputs
remain independent of OG mapping.

The optional [MonoPhy review](taxonomy_audit.md) writes `taxonomy_audit/` under
each selected species-tree branch. It includes taxonomic group results,
intruder/outlier events, associated run IDs and one review PDF per rank.
It reads the species tree and taxonomy only; gene trees are not audit inputs.

[`phylogeny_contrast_pairs`](contrast_pairs.md#pairs-after-full-or-phenotyped-inference)
writes `contrast/` under each selected molecular-tree branch, including pair
and species tables, source/assignment records, the observed subtree and summary
figures. Manual exclusion recomputes these results under
`filtered/phylogeny/all/contrast/` and `filtered/phylogeny/phenotyped/contrast/` when
inputs are complete. Original trees and the representative analysis at
`results/<analysis>/phylogeny/representatives/contrast/` remain unchanged.

## TPM interpretation

The optional [`filter_species` export](species_filter.md) writes a curated subset
under `results/<analysis>/filtered/`. Its `manifest.json` records the top-level
`exclude_species` list and source/output hashes. Original outputs remain available.

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

## Optional KEGG outputs

The [KEGG branch](kegg.md) writes `results/<analysis>/kegg/`, separately from the
OG tables above. Its `ko_tpm_sum.tsv` contains sums of original input TPM and is
**not renormalized** to the retained KO set. `ko_support.tsv` records annotated
and quantified gene counts, including KOs with no quantified genes. Such KOs have
blank values in support/wide tables and are omitted from the numeric long table;
observed zero expression remains zero. Run identities are preserved.

By default, `kegg.ambiguity: duplicate` adds a gene's full TPM to every distinct
accepted KO. KO features can therefore overlap and their total can exceed the
input TPM. In KEGG `mapping_qc.tsv`, `retained_targets`, `retained_tpm`, and
`retained_tpm_fraction` count each contributing gene once; `quantified_assignments`
and `ko_tpm_sum` count all gene/KO contributions. `ambiguous_tpm` is included in
retained TPM under `duplicate`; it is excluded only when `drop` is selected.

## Optional OG alignments

The [alignment branch](alignments.md) writes `results/<analysis>/orthogroups/alignments/`:

- `{og}.faa`: untrimmed protein MSA, one row per original gene ID, with all copies.
- `provenance.json`: alignment hashes and links to collection/execution records.

OG membership comes from the filename. Gene IDs use `{species}_g{number}`, so
removing the final `_g{number}` recovers the exact metadata species ID. There is
no separate `members.tsv` or added `species=` header attribute.

Every OG observed in the selected species' ODB mappings is included, even with
one sequence or no variation. Genes assigned to several OGs occur in each, regardless
of `tpm.multimap`. Species, gene, OG and site selection for modelling belongs in
PhenoRadar. No trimming or site-coordinate tables are produced.
