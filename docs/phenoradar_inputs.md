# PhenoRadar inputs

[Documentation](index.md)

`phenoradar_inputs` automatically links available completed results into
`results/<run_name>/phenoradar_inputs/`. No collection settings are needed.
The collection includes candidates for future PhenoRadar inputs as well as
currently usable tables and alignments. It does not start analyses or download
references. Keep linked source files available; failed validation preserves the
previous collection.

## Collecting results

After preparing metadata and running any desired analyses:

```bash
./run_pipeline.sh --cores 1 --resources mem_gb=4 -- phenoradar_inputs
```

Collection needs `metadata/samples.tsv` and `metadata/species_metadata.tsv`,
normally created during preparation; see
[backfilling older results](migration.md#backfilling-phenoradar-metadata).
Expression results are optional. Rerun the target after more analyses finish;
it discovers new outputs and removes links to results that are no longer present.
The collection log lists ready, absent, and incomplete sections.

## Published files

The original result paths are preserved, with short aliases for common inputs:

```text
results/<run_name>/phenoradar_inputs/
  species_metadata.tsv             # Alias of metadata/species_metadata.tsv
  tpm.tsv                          # Alias of orthogroups/expression/tpm.tsv
  alignments/{og}.faa               # Aliases of completed OG alignments
  metadata/                        # Selected samples, traits, taxonomy, and QC
  proteins/                        # Completed per-species translations
  orthogroups/
    expression/                    # Normalized/raw TPM, long/wide tables, and QC
    mapping/                       # Merged gene-to-OG tables and database
    alignments/                    # Alignments, membership, and completion record
  kegg/                            # KO expression, annotations, support, and QC
    ko_modules.tsv                 # Both maps collected when available
    ko_pathways.tsv
    species/                       # Completed per-species KO annotations
  phylogeny/
    all/                           # Full species set
    phenotyped/                    # Species with observed traits
    representatives/               # Representative species
```

Each available phylogeny branch includes completed molecular and dated trees,
gene trees and their marker alignments, selections, rooting information,
TimeTree calibrations, contrast pairs, and taxonomy reports. Completed BUSCO
sequence extraction and alignment steps are collected even before tree inference.
Native solver work
folders and unfinished results are excluded. Completion records and QC accompany
the data. An absent grouping map or optional analysis does not hide other
completed results.

Trees retain their names, such as `phylogeny/all/species_tree.nwk` and
`phylogeny/all/dating/species_tree.dated.nwk`; collection does not select one tree
for downstream analysis. All completed contrast branches retain their own
`contrast/species_metadata.tsv` and pair IDs. Base `species_metadata.tsv` stays
unchanged, so pair IDs from different analyses are never combined.

Existing `orthogroup_annotations.tsv[.gz]` files at the run root, in
`orthogroups/`, or in `orthogroups/mapping/` are also collected. These are
headerless OG/taxid/description tables; arbitrary external annotation paths are
no longer configured here. Descriptions are not downloaded by collection.

## Data checks and downstream use

Collection checks sample identities, expression coordinates and values, OG
alignment inventories, recorded checksums, and species-tree tips. Full trees
must cover the selected species; phenotyped and representative trees may cover
subsets. Corrupt completed inputs stop publication.

Multiple runs per species are preserved as separate expression rows. Collection
does not select or average replicates. Choose the appropriate inputs and handle
replicates according to the downstream analysis; collection alone does not mean
that every file can already be read by PhenoRadar.

OG `tpm.tsv` is normalized to one million per run. KO values sum original TPM and
can overlap across KOs; review [KO support](kegg.md#outputs). For KO expression,
set these options in PhenoRadar:

```yaml
data:
  tpm_path: results/run001/phenoradar_inputs/kegg/ko_tpm_sum.tsv
  feature_col: ko
  value_col: tpm_sum
  orthogroup_annotation_path: null
```

Replace `run001` with your run name. Alignments retain all copies and columns;
see [gene IDs and sequence format](alignments.md#outputs-and-phenoradar).

## Species exclusions

With nonempty `exclude_species`, run `filter_species` first. Collection uses only
the matching completed `filtered/` snapshot and verifies its recorded outputs;
it does not refresh that snapshot or mix in unfiltered analyses. Available
pruned trees and recomputed contrast results keep their paths. Clearing the
exclusion list selects the original data on the next collection.

For downstream setup, see PhenoRadar's
[quick start](https://github.com/mkrg01/phenoradar/blob/main/docs/quickstart.md) and
[data formats](https://github.com/mkrg01/phenoradar/blob/main/docs/data-format.md).
