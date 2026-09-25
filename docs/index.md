# Documentation

For manually curated RNA-seq additions/removals and staged Slurm jobs, see
[incremental datasets](datasets.md).

## Start a dataset

1. Check the [input formats and species selection](inputs.md).
2. Edit the [configuration](configuration.md); the matching [container image](containers.md) is selected automatically.
3. [Run the workflow](running.md); optionally start with a [pilot](running.md#pilot-run).
4. Inspect the [expression tables and QC](outputs.md#tpm-interpretation).

[Workflow targets](running.md#targets) · [Reference data](references.md) ·
[Output layout](outputs.md#directory-layout)

## Choose an analysis

| Task | Guide |
| --- | --- |
| Annotate proteins and quantify KEGG Orthology features | [KO expression](kegg.md) |
| Align every mapped orthogroup, keeping all gene copies | [OG protein alignments](alignments.md) |
| Infer species trees from existing BUSCO results | [BUSCO phylogeny](phylogeny.md) |
| Estimate divergence ages with treePL and TimeTree or manual calibrations | [Dating](dating.md) |
| Review tree placements against registered taxonomy | [Taxonomic review](taxonomy_check.md) |
| Select trait contrast pairs from molecular trees or an inferred representative tree | [Contrast pairs](contrast_pairs.md) |

## Prepare downstream inputs

After the required analyses finish, use [manual species exclusion](species_filter.md)
to export a curated subset if needed, then [collect PhenoRadar inputs](phenoradar_inputs.md).

## Maintenance

- [Container setup and releases](containers.md)
- [Releasing a version](releases.md)
- [Development and tests](development.md)

[Project README](../README.md)
