# Documentation

## Prepare a dataset

1. [Install phenoradar_prep](running.md#installation-and-normal-execution).
2. Prepare [input files](inputs.md) and edit the [configuration](configuration.md).
3. [Build a database and run an analysis](datasets.md).
4. Review [outputs and QC](outputs.md) and configure [PhenoRadar](phenoradar_inputs.md).

## Choose an analysis

| Task | Guide |
| --- | --- |
| Align mapped proteins by orthogroup | [OG alignments](alignments.md) |
| Annotate orthogroups and summarize KO expression | [KEGG](kegg.md) |
| Infer a species tree from BUSCO markers | [Phylogeny](phylogeny.md) |
| Estimate divergence ages | [Dating](dating.md) |
| Review tree placements against NCBI taxonomy | [Taxonomy check](taxonomy_check.md) |
| Identify contrasting sister-clade pairs and phylogenetic blocks | [Contrast pairs and blocks](contrast_pairs.md) |

## Manage data and runs

- [Plans, status and retries](datasets.md#plans-and-progress) · [Targets and resources](running.md#targets)
- [Reuse a database](datasets.md#reusing-completed-databases) · [Copy a database](datasets.md#copying-a-completed-build-to-another-project)
- [Export a sample subset](species_filter.md) · [Reference data](references.md)
- [Curated datasets and metadata updates](../datasets/README.md)
- [Development and tests](development.md) · [Releases](releases.md)
