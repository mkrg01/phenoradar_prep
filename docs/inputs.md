# Input data and species selection

[Documentation](index.md) · [Configuration](configuration.md)

Start with manually curated AMALGKIT metadata, one run per species. Build obtains
NCBI or local reads and uses GeneGalleon to generate CDS, BUSCO, and quantification.
Completed products go in **`builds/<build>/products/`**, not top-level `input/`.

## File formats

### Files you prepare

Metadata is required for a new build. The other files are optional, depending on
which species you want to include and which analyses you run. Build selects metadata
and exclusions; analysis auxiliary paths are under `inputs` in `analysis.yaml`.

| File | When to use | Format and settings |
| --- | --- | --- |
| `input/metadata.tsv` | **Required:** define the species and RNA-seq runs to build. | TSV: `scientific_name`, `run`, positive NCBI `taxid`; one run per species |
| `config/excluded_accessions.tsv` | **Optional:** exclude unusable or misidentified runs from builds while retaining their metadata. | TSV: `accession`, optional `reason`; `excluded_accessions` in build settings |
| `input/species_trait.tsv` | **Optional:** provide traits for PhenoRadar; required for phenotyped/representative trees and contrast pairs. | TSV: `species` and chosen trait column; `inputs.species_trait` and `trait` |
| `input/species_list.txt` | **Optional:** restrict an analysis to a chosen set of candidate species; BUSCO filtering still applies. | Species IDs, one per line; `inputs.species_list` and `selection.species_list: true` |
| `input/calibrations.tsv` | **Optional:** supply your own age bounds for [dating](dating.md#manual-calibrations) instead of TimeTree calibrations. | TSV: `taxa`, `min_age_ma`, `max_age_ma`, `source`; `inputs.calibrations` and `phylogeny.dating.calibration_source: file` |

Optional files with a configured path must exist; set unused paths to `null`.
Without a species list, use `selection.species_list: false`; without traits,
disable trait-dependent analyses. The supplied configs already specify paths
for the exclusion list and trait table.

Use [AMALGKIT-compatible metadata](https://github.com/kfuku52/amalgkit/wiki/amalgkit-metadata);
additional columns pass through to GeneGalleon. For local reads, also provide
`private_file` = `yes`, `lib_layout` = `single` or `paired`, `read1_path`, and
`read2_path` for paired reads. These are TSV columns. Relative FASTQ paths resolve
against the metadata directory; use a new run ID when read content changes.

### Generated products

| Path under `builds/<build>/products/` | Format |
| --- | --- |
| `metadata.tsv` | Frozen metadata after run exclusions |
| `busco/summary.tsv` | TSV: `Species`, `busco_cds_single`, `busco_cds_duplicated`, `busco_cds_fragmented`, `busco_cds_missing`, `busco_cds_total` |
| `busco/full/{species}.busco.full.tsv` | Full BUSCO table with lineage header, for every build species |
| `cds/{species}_longestCDS.fa.gz` | Gzip FASTA, one per species |
| `quant/{species}/{run}/{run}_abundance.tsv` | TSV: `target_id`, `tpm` |
| `proteins/{odb_species}_protein.fa` | Translated protein FASTA |
| `odb/` | Mapping database, annotations, and reusable snapshot |

Build creates or reuses these automatically. See [output layout](outputs.md#directory-layout)
for intermediate workspaces and [portable builds](datasets.md#copying-a-completed-build-to-another-project)
for copying the complete bundle.

### Importing existing products

Fresh RNA-seq builds need no registration. For legacy products, use the CDS,
BUSCO, and quant paths above relative to an import directory:

```bash
./run_build.sh register --input-dir imports/legacy --metadata input/metadata.tsv
```

BUSCO full tables also accept `{species}_busco.full.tsv`; either form may have
`.gz`. A summary alone cannot complete BUSCO. Keep registered source files inside
the project. Old `input/cds/`, `input/busco/`, and `input/quant/` layouts remain
supported with `--input-dir input`, but new builds do not write there.

Copied bundles use `register --products <directory>` to seed a new build;
analysis consumes them directly. See [importing work](datasets.md#updating-species-and-importing-existing-work)
for ODB snapshot registration.

## Identifiers

- Metadata needs unique species/run IDs and consistent taxids.
- Species IDs replace spaces with underscores. Hyphens remain in IDs but become
  underscores in ODB filenames: `Beta sp-X` becomes `Beta_sp-X` / `Beta_sp_X`.
  Avoid collisions between these forms.
- FASTA IDs must be unique across species and match abundance `target_id`.
  [OG alignments](alignments.md#input-requirements) require `{species}_g{number}`.
- ODB-mapper paths must not contain spaces or shell metacharacters.

## Species selection

Build must complete all included species; missing products do not automatically
exclude them. Remove unwanted runs with the
[manual accession list](datasets.md#manually-excluding-unusable-accessions) or metadata edits.

Analysis applies the optional species list, `exclude_species`, and BUSCO threshold
`(single + duplicated) / total >= selection.busco_threshold` (default `0.5`).
Unknown species IDs, invalid counts, duplicate identities, or an empty selection
are errors. Metadata columns such as `exclusion` do not filter species.
Unresolved taxids fail unless `selection.missing_taxonomy: allow`; missing
individual ranks are allowed.

Review `results/<analysis>/metadata/selection.json`, `samples.tsv`, and
`busco_completeness.svg` for selection reasons, sample paths, and BUSCO QC.

## Traits

`inputs.species_trait` supplies traits; trait columns in sample metadata are ignored.

```tsv
species	carnivory
Plant alpha	0
Plant beta	1
Plant gamma	NA
```

Species names normalize as above; duplicate normalized names are errors. Blanks,
`NA`, `NaN`, `nan`, and absent rows are unknown. Extra rows do not add species.
Set `inputs.species_trait: null` when unused.

Top-level `trait` selects the shared column for tree selection, contrast pairs,
and PhenoRadar metadata. PhenoRadar traits must be `0`, `1`, or blank; the `C4`
column also enforces these values. Phenotyped inference accepts a single observed
state; representative selection requires two. On all/phenotyped trees, fewer than
two states produces zero contrast pairs.
