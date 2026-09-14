# Input data and species selection

[Documentation](index.md) · [Configuration](configuration.md)

The workflow takes existing CDS assemblies, expression estimates, and BUSCO
results. Dataset files conventionally live under `input/`; configured external
paths and symbolic links can be used to avoid copying large files.

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
Manual calibration tables are described in [dating](dating.md#manual-calibrations).

## Identifiers

Run IDs must be unique in metadata; species names must be unique in the BUSCO
summary. Multiple runs from one species are supported, and must agree on its
taxid. Runs remain separate throughout expression aggregation.

`{species}` is the scientific name with spaces replaced by underscores. Hyphens
remain in species IDs, but become underscores in the `odb_species` filenames
used by ODB-mapper. Name collisions are rejected. For example,
`Beta sp-X` becomes species ID `Beta_sp-X` and ODB filename prefix `Beta_sp_X`.

FASTA record IDs are preserved. They must be unique across species and match
the abundance table's `target_id` values. The [OG alignment branch](alignments.md)
additionally requires IDs of the form `{species}_g{number}`, so the species can
be recovered from each sequence. For example, `Beta_sp-X_g12` belongs to
`Beta_sp-X`.

ODB-mapper paths must not contain spaces or shell metacharacters because its
upstream shell scripts cannot handle them.

## Species selection

Species are retained when their complete BUSCO fraction meets
`selection.busco_threshold` (default `0.5`):

```text
complete fraction = (single + duplicated) / total
```

`selection.species_list` can further restrict the dataset. It is a text file
with one exact species ID per line. Every requested species must occur in the
metadata, have BUSCO data, and pass the threshold.

Species absent from the BUSCO summary are excluded automatically. Their rows
remain in `metadata_all.tsv` with `selected=False` and blank BUSCO values.
`selection.json` records `missing_busco_species` and `missing_busco_runs`.
These species require no CDS, abundance file, or taxonomy lookup, and are omitted
from the BUSCO histogram.

Invalid BUSCO counts, duplicate runs, or missing CDS/abundance files for selected
samples stop preparation. This input validation also applies when preparing a
phylogeny-only analysis. Unresolved taxids for species with BUSCO data fail by
default; `selection.missing_taxonomy: allow` permits them. Missing individual
taxonomic ranks are allowed. Metadata columns such as `exclusion` do not filter
the dataset; post-analysis removal uses [manual species exclusion](species_filter.md).

Preparation writes the full metadata with a selection flag, the selected rows,
and `metadata/samples.tsv` with exact CDS and abundance paths. Selection is by
species; the BUSCO histogram counts runs.

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

Base PhenoRadar metadata can be prepared without a trait file; traits are then
blank and the metadata log records this. Trait-based analyses require the file.
Extra species in the trait table do not enter the selected dataset.
