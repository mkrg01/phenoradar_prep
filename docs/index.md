# Documentation

Prepare a reusable RNA-seq database, select samples for an analysis, and export
results to PhenoRadar.

## Getting started

1. [Install the workflow](running.md#installation-and-normal-execution).
2. Prepare [input files](inputs.md) and edit the [configuration](configuration.md).
3. [Build a database and run an analysis](datasets.md).
4. Review [outputs and QC](outputs.md), then use the [PhenoRadar inputs](phenoradar_inputs.md).

## Choose an analysis

| Task | Guide |
| --- | --- |
| Align all gene copies in each mapped OG | [OG alignments](alignments.md) |
| Annotate OGs and quantify KO expression | [KEGG](kegg.md) |
| Infer trees from BUSCO markers | [Phylogeny](phylogeny.md) |
| Estimate divergence ages | [Dating](dating.md) |
| Review tree placements against NCBI taxonomy | [Taxonomy check](taxonomy_check.md) |
| Select trait contrast pairs | [Contrast pairs](contrast_pairs.md) |

## Managing runs

- [Slurm resources, targets, and retries](running.md)
- [Reuse or copy a database](datasets.md#reusing-completed-databases)
- [Export a subset of completed results](species_filter.md)
- [Reference data and updates](references.md)
- [Curated metadata datasets and updates](../datasets/README.md)

For contributors: [development and tests](development.md) · [releases](releases.md).

[Project README](../README.md)
