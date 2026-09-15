# Input data and species selection

[Documentation](index.md) · [Configuration](configuration.md)

The workflow takes existing CDS assemblies, expression estimates, and BUSCO
results. Dataset files conventionally live under `input/`; configured external
paths and symbolic links can be used to avoid copying large files.

## Upstream data preparation

[AMALGKIT](https://github.com/kfuku52/amalgkit)'s
[`metadata` command](https://github.com/kfuku52/amalgkit/wiki/amalgkit-metadata)
retrieves RNA-seq sample metadata. [GeneGalleon](https://github.com/kfuku52/genegalleon)'s
[transcriptome generation stage](https://github.com/kfuku52/genegalleon/blob/main/docs/main-stages-and-what-they-do.md#gg_transcriptome_generation_entrypointsh)
produces assemblies, longest CDS, expression estimates, and BUSCO results.

| Upstream output | Input in phenoradar_prep | Use here |
| --- | --- | --- |
| AMALGKIT sample metadata | `input/metadata.tsv` | Species/run identities and taxids |
| GeneGalleon longest-CDS FASTA files | `input/cds/{species}_longestCDS.fa.gz` | Translation and downstream sequence analyses |
| GeneGalleon expression estimates | `input/quant/{species}/{run}/{run}_abundance.tsv` | Input TPM for OG and KO aggregation |
| GeneGalleon BUSCO completeness summary | `input/busco/summary.tsv` | Species selection by complete BUSCO fraction |
| GeneGalleon per-species BUSCO full tables | `input/busco/full/{species}.busco.full.tsv` | Existing marker assignments for optional species-tree inference |

## File formats

```text
input/
  metadata.tsv
  species_trait.tsv
  busco/
    summary.tsv
    full/{species}.busco.full.tsv  # For phylogeny
  cds/{species}_longestCDS.fa.gz
  quant/{species}/{run}/{run}_abundance.tsv
  calibrations.tsv                # For manual dating
  pilot_species.txt               # Optional subset
```

| Input | Required columns or format |
| --- | --- |
| Sample metadata | TSV: `scientific_name`, `run`, `taxid` |
| BUSCO summary | TSV: `Species`, `busco_cds_single`, `busco_cds_duplicated`, `busco_cds_fragmented`, `busco_cds_missing`, `busco_cds_total` |
| CDS | Gzip-compressed FASTA, one file per species |
| Abundance | TSV: `target_id`, `tpm`, one file per run |
| Species traits | TSV: `species` and the chosen trait column; required for trait-based analyses |

Per-species BUSCO full tables are needed only for [phylogeny](phylogeny.md#inputs).
Standard filenames and their gzip versions are detected automatically; the
phylogeny guide lists supported layouts. Sequence input uses the same CDS files.
Manual calibration tables are described in [dating](dating.md#manual-calibrations).

## Identifiers

- Run IDs must be unique. Multiple runs per species are supported and must agree
  on taxid; expression is kept per run. BUSCO summary species must be unique.
- Species IDs replace spaces with underscores. Hyphens remain in species IDs but
  become underscores in ODB filenames: `Beta sp-X` becomes `Beta_sp-X` / `Beta_sp_X`.
- FASTA IDs must be unique across species and match abundance `target_id` values.
  [OG alignments](alignments.md) require `{species}_g{number}`, e.g. `Beta_sp-X_g12`.
- ODB-mapper paths must not contain spaces or shell metacharacters.

## Species selection

`selection.species_list` limits the candidate species: `null` considers all species;
a text file path (e.g. `input/pilot_species.txt`) lists one ID per line
(e.g. `Arabidopsis_thaliana`). Unknown IDs are errors.
Candidates then need a complete BUSCO fraction `(single + duplicated) / total`
of at least `selection.busco_threshold` (default `0.5`). Listed species below
the threshold are excluded. An empty selection is an error.
Species without BUSCO summary data are excluded and recorded in `selection.json`.

Selected samples need valid CDS and abundance files, including for phylogeny-only
runs. Taxonomy is looked up for selected species. Invalid counts and duplicate
identities stop preparation. Unresolved taxids fail by default;
`selection.missing_taxonomy: allow` permits them. Missing
individual ranks are allowed. Metadata columns such as `exclusion` do not filter
species; use [manual exclusion](species_filter.md) after analysis.

Preparation writes full/selected metadata and exact input paths in `samples.tsv`.
Selection is by species; the BUSCO histogram counts runs.

## Traits

`inputs.species_trait` is the phenotype source. Sample metadata trait columns
are ignored. Its default path is `input/species_trait.tsv`:

```tsv
species	C4
Plant alpha	0
Plant beta	1
Plant gamma	NA
```

Spaces in species names are normalized to underscores; duplicate normalized
names are rejected. Blanks, `NA`, `NaN`, `nan`, and absent species rows are
unknown. `C4` accepts `0`, `1`, or missing values.

`phylogeny.trait` selects the nonmissing species for phenotyped inference.
`contrast.trait` selects the column for pair assignment and PhenoRadar metadata;
both default to `C4`. Phenotyped inference can use other traits and does not
require two states. Pair assignment supports two observed states, while the
PhenoRadar metadata contract requires `0`, `1`, or blank.

A trait file is required for trait-based analyses. Without it, base PhenoRadar
metadata has blank traits. Extra trait-table species do not enter the dataset.
