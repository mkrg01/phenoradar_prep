# PhenoRadar inputs

[Documentation](index.md)

`phenoradar_inputs` links available completed results into
`results/<run_name>/phenoradar_inputs/`. It needs no collection settings and
starts no analyses or downloads. Keep linked source files available.

## Collecting results

After running the desired analyses:

```bash
./run_pipeline.sh --cores 1 --resources mem_gb=4 -- phenoradar_inputs
```

Collection requires `metadata/samples.tsv` and `metadata/species_metadata.tsv`;
expression is optional. Rerun after more analyses finish. The log reports ready,
absent, and incomplete sections; obsolete links are removed. Failed validation
preserves the previous collection.

## Published files

Original result paths are preserved, with aliases for common inputs:

| Path under `phenoradar_inputs/` | Contents |
| --- | --- |
| `species_metadata.tsv` | Alias of `metadata/species_metadata.tsv` |
| `tpm.tsv` | Alias of `orthogroups/expression/tpm.tsv` |
| `alignments/{og}.faa` | Aliases of completed OG alignments |
| `metadata/`, `proteins/` | Selection, traits, taxonomy/QC, and translations |
| `orthogroups/` | Expression, mappings, and alignments |
| `kegg/` | KO expression, annotations, support/QC, and available module/pathway maps |
| `phylogeny/{all,phenotyped,representatives}/` | Completed trees, marker alignments, calibrations, pairs, and taxonomy reports |

QC and completion records accompany results; unfinished outputs and native
solver work folders are omitted. Completed BUSCO extraction/alignment steps
can be collected before inference finishes.

Choose the downstream tree explicitly, e.g. `phylogeny/all/species_tree.nwk`
or `phylogeny/all/dating/species_tree.dated.nwk`. Each contrast branch retains
its own `contrast/species_metadata.tsv`; base metadata is unchanged and pair
IDs from separate analyses are never combined.

Existing headerless OG/taxid/description `orthogroup_annotations.tsv[.gz]` files
are collected from the run root, `orthogroups/`, or `orthogroups/mapping/`.
Collection does not download descriptions.

## Data checks and downstream use

Collection validates identities, expression values, alignment inventories,
checksums, and tree tips. Full trees must cover selected species; phenotyped and
representative trees may cover subsets. Corrupt completed inputs stop collection.

Multiple expression runs remain separate; select or aggregate replicates for
your downstream analysis. Some collected files are future input candidates;
collection does not guarantee PhenoRadar supports every file.

OG `tpm.tsv` is rescaled to one million per run. KO values sum original TPM and
can overlap across KOs; review [support and missing values](kegg.md#outputs).
For KO expression, configure PhenoRadar as follows (replace `run001`):

```yaml
data:
  tpm_path: results/run001/phenoradar_inputs/kegg/ko_tpm_sum.tsv
  feature_col: ko
  value_col: tpm_sum
  orthogroup_annotation_path: null
```

Alignments retain all copies/columns; see [gene IDs](alignments.md#outputs-and-phenoradar).
For downstream setup, see PhenoRadar's
[quick start](https://github.com/mkrg01/phenoradar/blob/main/docs/quickstart.md) and
[data formats](https://github.com/mkrg01/phenoradar/blob/main/docs/data-format.md).

## Species exclusions

With nonempty `exclude_species`, run `filter_species` first. Collection requires
the matching completed `filtered/` snapshot; it neither refreshes it nor mixes
in unfiltered analyses. Clearing exclusions selects original data on recollection.
