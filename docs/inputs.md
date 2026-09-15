# Input data and species selection

[Documentation](index.md) · [Configuration](configuration.md)

The workflow takes existing CDS assemblies, expression estimates, and BUSCO
results. Use external paths or symlinks to avoid copying large files; external
locations need [container bind mounts](containers.md#external-files).

## Upstream data preparation

[AMALGKIT metadata](https://github.com/kfuku52/amalgkit/wiki/amalgkit-metadata)
provides sample identities and taxids.
[GeneGalleon transcriptome generation](https://github.com/kfuku52/genegalleon/blob/main/docs/main-stages-and-what-they-do.md#gg_transcriptome_generation_entrypointsh)
provides longest CDS, expression estimates, and BUSCO results.

## File formats

Default paths are below; change `inputs.*` in the configuration as needed.

| Path under `input/` | Format |
| --- | --- |
| `metadata.tsv` | TSV: `scientific_name`, `run`, `taxid` |
| `busco/summary.tsv` | TSV: `Species`, `busco_cds_single`, `busco_cds_duplicated`, `busco_cds_fragmented`, `busco_cds_missing`, `busco_cds_total` |
| `cds/{species}_longestCDS.fa.gz` | Gzip FASTA, one per species |
| `quant/{species}/{run}/{run}_abundance.tsv` | TSV: `target_id`, `tpm`, one per run |
| `species_trait.tsv` | TSV: `species` and chosen trait column; required for trait-based analyses |

[Phylogeny](phylogeny.md#inputs) additionally needs per-species BUSCO full tables,
by default `busco/full/{species}.busco.full.tsv`.
[Manual dating](dating.md#manual-calibrations) needs a calibration TSV.

## Identifiers

- Run IDs must be unique. Multiple runs per species must agree on taxid;
  expression stays per run. BUSCO summary species must be unique.
- Species IDs replace spaces with underscores. Hyphens remain in IDs but become
  underscores in ODB filenames: `Beta sp-X` becomes `Beta_sp-X` / `Beta_sp_X`.
- FASTA IDs must be unique across species and match abundance `target_id`.
  [OG alignments](alignments.md#input-requirements) require `{species}_g{number}`,
  e.g. `Beta_sp-X_g12`.
- ODB-mapper paths must not contain spaces or shell metacharacters.

## Species selection

`selection.species_list: null` considers all species. To limit candidates, supply
a text file with one exact ID per line, e.g. `Arabidopsis_thaliana`; unknown IDs
are errors. Candidates then need a complete BUSCO fraction
`(single + duplicated) / total >= selection.busco_threshold` (default `0.5`).
Missing BUSCO data or scores below the threshold exclude species; an empty
selection is an error.

Selected samples need valid CDS and abundance files even for phylogeny-only
runs. Invalid counts or duplicate identities stop preparation. Unresolved taxids
fail unless `selection.missing_taxonomy: allow`; missing individual ranks are
allowed. Metadata columns such as `exclusion` do not filter species:
use [manual exclusion](species_filter.md) after analysis.

Results in `metadata/` include `selection.json` (exclusion reasons),
`samples.tsv` (selected samples and exact input paths), and
`busco_completeness.svg` (run counts; selection itself is by species).

## Traits

`inputs.species_trait` supplies traits; sample metadata trait columns are ignored.

```tsv
species	C4
Plant alpha	0
Plant beta	1
Plant gamma	NA
```

Spaces normalize to underscores; duplicate normalized species names are errors.
Blanks, `NA`, `NaN`, `nan`, and absent rows are unknown. `C4` accepts only
`0`, `1`, or missing values. Extra species do not enter the dataset.

`phylogeny.trait` defines the nonmissing species for phenotyped inference;
`contrast.trait` selects the column for pairs and PhenoRadar metadata.
Both default to `C4`. Phenotyped inference accepts other traits, including
single-state traits; pairs support two observed states, and PhenoRadar metadata
requires `0`, `1`, or blank. Without a trait file, base metadata has blank traits.
