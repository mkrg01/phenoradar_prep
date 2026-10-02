# Inputs and sample selection

[Documentation](index.md) · [Configuration](configuration.md)

A new build starts from AMALGKIT metadata with one run per sample. Multiple runs
may belong to the same biological species; each is assembled and quantified
separately.

## File formats

| File | Required columns or format |
| --- | --- |
| `input/metadata.tsv` | TSV: `scientific_name`, unique `run`, positive NCBI `taxid` |
| `config/excluded_accessions.tsv` | TSV: `accession`, optional `reason` |
| `input/species_trait.tsv` | TSV: `species` and the chosen trait column |
| `input/species_list.txt` | One biological species ID or exact sample ID per line |
| `input/calibrations.tsv` | TSV: `taxa`, `min_age_ma`, `max_age_ma`, `source`; see [dating](dating.md#manual-calibrations) |

Only metadata is always required. Other files depend on your analysis;
configured paths must exist, so set unused optional paths to `null`. The supplied
configs already specify an exclusion list and trait table.

The metadata path is configurable and resolves against the repository root.
For example, `metadata: datasets/angiosperm_leaf/metadata.tsv` reads a curated
dataset directly; copying it to `input/` is unnecessary. The
[dataset catalog](../datasets/README.md) describes acquisition, curation, and
representative updates independently of database execution.

Additional AMALGKIT metadata columns pass through to GeneGalleon. For local
reads, also provide these TSV columns:

| Column | Value |
| --- | --- |
| `private_file` | `yes` |
| `lib_layout` | `single` or `paired` |
| `read1_path` | FASTQ path |
| `read2_path` | Second FASTQ path for paired reads |

Relative FASTQ paths resolve against the metadata directory. Use a new run ID
when read content changes. ODB-mapper paths must avoid spaces and shell metacharacters.

## Identifiers

For `Abelia chinensis` and run `SRR14320411`, the workflow derives:

| Identifier | Value | Meaning |
| --- | --- | --- |
| `species_id` | `Abelia_chinensis` | Biological species |
| `analysis_sample_id` | `Abelia_chinensis_SRR14320411` | Separately processed sample |
| Gene ID | `Abelia_chinensis_SRR14320411_g0` | Assembled gene; matches abundance `target_id` |

Output tables use the column **`species` for the analysis sample ID**. Tree tips
use the same ID. Join through `metadata/samples.tsv` to recover `species_id`,
`scientific_name`, `taxid`, and `run`; do not split sample IDs, since run IDs may
contain underscores. The original AMALGKIT `sample_id` is preserved separately.

In file paths, `{species}` means the analysis sample ID; `{odb_species}` replaces
its hyphens with underscores. Colliding IDs are rejected.
Separate runs remain separate observations and are not automatically independent
biological replicates.

## Species selection

Build must complete every included sample. Remove unwanted runs through metadata
or the [accession exclusion list](datasets.md#manually-excluding-unusable-accessions).

Analysis applies `inputs.species_list`, `exclude_species`, and the per-sample
BUSCO criterion `(single + duplicated) / total >= selection.busco_threshold`
(default `0.5`). Biological species IDs select/exclude all their samples; exact
sample IDs affect only that sample. Metadata columns such as `exclusion` do not
filter samples.

Unknown IDs, invalid counts, duplicate identities, unresolved taxids, or an empty
selection are errors; missing individual taxonomy ranks are allowed. Review
`metadata/selection.json`, `samples.tsv`, and `busco_completeness.svg` in the analysis.

## Traits

Traits come from `inputs.species_trait`, not sample metadata:

```tsv
species	carnivory
Plant alpha	0
Plant beta	1
Plant gamma	NA
```

Set `trait: carnivory` to use this column. Species names normalize to biological
species IDs, and every sample of that species inherits its trait. Duplicate
normalized names are errors. Blanks, `NA`, `NaN`, `nan`, and absent rows are
unknown; extra rows do not add species.

Use `0`, `1`, or blank for PhenoRadar traits. Phenotyped trees accept a single
observed state; representative selection requires two. On all/phenotyped trees,
fewer than two states produces zero contrast pairs. If traits are unused, set
`inputs.species_trait: null` and disable trait-dependent branches.
