# Documentation

## Start a dataset

1. Check the [input formats and species selection](inputs.md).
2. Create a [configuration](configuration.md) and review the OrthoDB node.
3. Run [preparation and a pilot](running.md#prepare-and-inspect) before the full dataset.
4. Inspect the [expression tables and QC](outputs.md#tpm-interpretation).

[Workflow targets](running.md#targets) · [Reference data](references.md) ·
[Output layout](outputs.md#directory-layout)

## Choose an analysis

| Task | Guide |
| --- | --- |
| Annotate proteins and quantify KEGG Orthology features | [KO expression](kegg.md) |
| Align every mapped orthogroup, keeping all gene copies | [OG protein alignments](alignments.md) |
| Infer species trees from existing BUSCO results | [BUSCO phylogeny](phylogeny.md) |
| Date a species tree with manual or TimeTree calibrations | [Dating](dating.md) |
| Review tree placements against registered taxonomy | [Taxonomic review](taxonomy_audit.md) |
| Select trait contrast pairs from molecular trees or an inferred representative tree | [Contrast pairs](contrast_pairs.md) |

## Prepare downstream inputs

After the required analyses finish, use [manual species exclusion](species_filter.md)
to export a curated subset if needed, then [collect PhenoRadar inputs](phenoradar_inputs.md).

## Maintenance

- [Container releases and Apptainer](containers.md)
- [Development and tests](development.md)
- [Updating older configurations and results](migration.md)

[Project README](../README.md)
