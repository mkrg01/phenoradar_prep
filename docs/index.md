# Documentation

## Start from metadata

1. Prepare [metadata and auxiliary inputs](inputs.md).
2. Edit [build/analysis settings](configuration.md) and set up [containers](containers.md).
3. [Build through ODB mapping, then run analysis](datasets.md).
4. Review [outputs and QC](outputs.md) and use the [PhenoRadar inputs](phenoradar_inputs.md).

[Slurm, retries, and targets](running.md) · [Reference data](references.md)

## Choose an analysis

Set options before preparing a named analysis. Each guide shows how to submit
its target separately; missing prerequisites are included automatically.

| Task | Guide |
| --- | --- |
| Annotate proteins and quantify KO features | [KO expression](kegg.md) |
| Align mapped orthogroups, retaining all copies | [OG alignments](alignments.md) |
| Infer species trees from BUSCO markers | [Phylogeny](phylogeny.md) |
| Estimate divergence ages with treePL | [Dating](dating.md) |
| Review placements against NCBI taxonomy | [Taxonomic review](taxonomy_check.md) |
| Select trait contrast pairs | [Contrast pairs](contrast_pairs.md) |

## Prepare downstream inputs

`run_analysis.sh` collects PhenoRadar inputs after successful `all` execution.
After running more branches, [collect again](phenoradar_inputs.md).
For exclusions, normally prepare a new analysis with `exclude_species`;
[post hoc filtering](species_filter.md) can instead export a subset of completed results.

## Maintenance

[Container setup](containers.md) · [Releasing a version](releases.md) ·
[Development and tests](development.md)

[Project README](../README.md)
