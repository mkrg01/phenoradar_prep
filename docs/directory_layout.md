# Directory layout

[Back to README](../README.md)

Dataset inputs live under `input/`. Results are grouped by purpose, with OG
mapping, expression, and alignments together and each species-tree analysis in
its own branch. Workflow rules, species filtering, and the PhenoRadar collector
use the same path definitions.

```text
input/
  metadata.tsv
  species_trait.tsv
  cds/
  quant/
  busco/
    summary.tsv
    full/                         # Optional per-species phylogeny inputs
  calibrations.tsv                # Optional dating input
  pilot_species.txt               # Optional species selection
config/                           # Settings, not dataset tables
resources/                        # Reusable reference snapshots and tools
results/<analysis>/
  run.json
  metadata/
  proteins/
  orthogroups/
    mapping/
      mappings.sqlite
      gene_orthogroups.tsv
      merge_qc.json
      manifests/
      chunks/
    expression/
    alignments/
  kegg/
  phylogeny/
    all/
    phenotyped/
    representatives/
  phenoradar_inputs/               # Only files selected for PhenoRadar
  filtered/                       # Corresponding species-filtered results
work/<analysis>/
logs/<analysis>/
```

Optional directories are created only when used. `input/` is a convention, not
a requirement to copy large datasets: configured external paths and symbolic
links remain appropriate. Preserve the species/run structure inside `quant/`
and the gene identifiers inside sequence files. The active sample metadata is
`input/metadata.tsv`; a retained source version may be stored alongside it as
`input/metadata.original.tsv`.

`config/` holds analysis settings and `resources/` holds reusable reference
snapshots. Completed results, temporary execution files, and diagnostic records
remain in `results/`, `work/`, and `logs/`, respectively. Work and log branches
follow the same `orthogroups/` and `phylogeny/` organization as their producers.
Runtime and memory benchmark TSVs belong under the corresponding log branch's
`benchmarks/` directory.

## Output organization

OG mapping files are directly inside `orthogroups/mapping/`, without a redundant
`merged/` level. Chunk plans and native mapper results remain in its
`manifests/` and `chunks/` subdirectories. OG expression tables live in
`orthogroups/expression/` and all-copy OG FASTA alignments in
`orthogroups/alignments/`.

Each `phylogeny/` branch stores its species tree directly in the branch
directory, beside marker plans, alignments, gene trees, and `rooting/`. The
`phenotyped/` and `representatives/` branches also keep their own `selection/`
manifests. Dating, taxonomy reports, and contrast-pair results live in
`dating/`, `taxonomy_audit/`, and `contrast/` beside the tree that produced them,
when those analyses are run. The dataset-wide NCBI guide in
`phylogeny/all/rooting/` also supplies representative selection; representatives
resolve their own outgroup in `phylogeny/representatives/rooting/`.

KEGG remains independent of OG mapping. Its KO expression tables and
KO-to-module/pathway tables support KO features and functional groups. No
OG-to-KEGG tables are generated.

`phenoradar_inputs/` is a small view of selected completed results, with links
to expression, sequence, annotation, and tree files. Its concise published
filenames are described in the [PhenoRadar input guide](phenoradar_inputs.md).
Filtered inputs come from the matching `filtered/` dataset.

## Existing-data migration

Old directory aliases are not retained. Existing data are relocated with
filesystem renames, avoiding another copy of the large CDS, expression, protein,
and mapping datasets. Migration records and original versions of rewritten
operational files are kept under `logs/layout_migration/<timestamp>/`.

Selected-sample manifests contain absolute CDS and abundance paths, so those
operational paths must follow the input move. Existing PhenoRadar links must
also follow the moved producer outputs. The migration records source and
destination locations and verifies retained files and links; changing a path
does not imply a new sequence analysis or expression calculation.

Historical `run.json` and provenance records remain unchanged as evidence of
the original execution. They may name former locations. Archived originals and
the migration record explain changed operational files without rewriting old
commands, timestamps, or execution claims. New jobs write records for their
current inputs and locations. Historical filtered snapshots retain their
recorded contents; regenerate filtered results from the relocated sources
before using them in a new analysis. See [outputs](outputs.md) and the
[phylogeny guide](phylogeny.md#outputs-and-validation) for the complete contents
of each result branch.
