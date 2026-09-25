# Input data and species selection

[Documentation](index.md) · [Configuration](configuration.md)

Start with manually curated AMALGKIT `input/metadata.tsv`, one run per species.
The [build phase](datasets.md) generates CDS, BUSCO full tables, and quantification
from NCBI or local reads. The layout below also supports importing existing work.
Keep imported artifacts inside the repository for container access.

## Upstream data preparation

[AMALGKIT metadata](https://github.com/kfuku52/amalgkit/wiki/amalgkit-metadata)
provides sample identities and taxids.
[GeneGalleon transcriptome generation](https://github.com/kfuku52/genegalleon/blob/main/docs/main-stages-and-what-they-do.md#gg_transcriptome_generation_entrypointsh)
provides longest CDS, expression estimates, and BUSCO results.

## File formats

Build creates this layout automatically. Use it when importing existing products;
analysis auxiliary paths are configured under `inputs` in `analysis.yaml`.

| Path under `input/` | Format |
| --- | --- |
| `metadata.tsv` | TSV: `scientific_name`, `run`, `taxid` |
| `busco/summary.tsv` | TSV: `Species`, `busco_cds_single`, `busco_cds_duplicated`, `busco_cds_fragmented`, `busco_cds_missing`, `busco_cds_total` |
| `cds/{species}_longestCDS.fa.gz` | Gzip FASTA, one per species |
| `quant/{species}/{run}/{run}_abundance.tsv` | TSV: `target_id`, `tpm`, one per run |
| `species_trait.tsv` | TSV: `species` and chosen trait column; required for trait-based analyses |
| `species_list.txt` | Optional candidate species IDs, one per line; enabled by `selection.species_list: true` |

[Phylogeny](phylogeny.md#inputs) additionally needs full tables under `busco/full/`.
[Manual dating](dating.md#manual-calibrations) reads `calibrations.tsv`.

## Identifiers

- Build metadata requires one unique run per species and consistent taxids.
  The low-level table utilities also support multiple runs; the build interface does not.
- Species IDs replace spaces with underscores. Hyphens remain in IDs but become
  underscores in ODB filenames: `Beta sp-X` becomes `Beta_sp-X` / `Beta_sp_X`.
- FASTA IDs must be unique across species and match abundance `target_id`.
  [OG alignments](alignments.md#input-requirements) require `{species}_g{number}`,
  e.g. `Beta_sp-X_g12`.
- ODB-mapper paths must not contain spaces or shell metacharacters.

## Species selection

`selection.species_list: false` considers all species. To limit candidates, set
it to `true` in analysis settings and put exact IDs in `inputs.species_list`, one per line,
e.g. `Arabidopsis_thaliana`. Unknown IDs are errors. Candidates need a complete BUSCO fraction
`(single + duplicated) / total >= selection.busco_threshold` (default `0.5`).
Missing BUSCO data or scores below the threshold exclude species; an empty
selection is an error.

Selected samples need valid CDS and abundance files even for phylogeny-only
runs. Invalid counts or duplicate identities stop preparation. Unresolved taxids
fail unless `selection.missing_taxonomy: allow`; missing individual ranks are
allowed. Metadata columns such as `exclusion` do not filter species:
use `exclude_species` in analysis settings, edit metadata, or apply the
[manual run exclusion list](datasets.md#manually-excluding-unusable-accessions)
when preparing a new build.

Results in `metadata/` include `selection.json` (exclusion reasons),
`samples.tsv` (selected samples and exact input paths), and
`busco_completeness.svg` (run counts; selection itself is by species).

## Traits

`input/species_trait.tsv` supplies traits; sample metadata trait columns are ignored.

```tsv
species	C4
Plant alpha	0
Plant beta	1
Plant gamma	NA
```

Spaces normalize to underscores; duplicate normalized species names are errors.
Blanks, `NA`, `NaN`, `nan`, and absent rows are unknown. `C4` accepts only
`0`, `1`, or missing values. Extra species do not enter the dataset.

Top-level `trait` selects one shared column for phenotyped/representative
selection, contrast pairs, and PhenoRadar metadata.
The current configuration uses `carnivory`. Phenotyped inference accepts
single-state traits; pairs support two observed states, and PhenoRadar metadata
requires `0`, `1`, or blank. Without a trait file, base metadata has blank traits.
