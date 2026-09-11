# PhenoRadar inputs

[Back to README](../README.md)

The explicit `phenoradar_inputs` target collects completed results into
`results/<analysis>/phenoradar_inputs/`. It validates the selected dataset and
links existing files instead of copying expression tables or protein alignments.
It does not start missing ODB, KEGG, alignment, filtering, or phylogeny jobs.
Each invocation validates and refreshes the directory, removing stale optional
links. A failed validation preserves the previous publication. Keep the source
results and external references available: this directory contains links, not
an independent data archive.

## Preparing and collecting

`prepare`, the default workflow, `kegg`, and `contrast_pairs` also prepare
`metadata/species_metadata.tsv`. This small
metadata step reads the selected sample/taxonomy tables and the supplied species
trait table. It runs independently of annotation or expression calculations.
For results created before that output existed, prepare it first:

```bash
./run_pipeline.sh --software-deployment-method conda \
  --configfile config/mydata.yaml --cores 1 -- phenoradar_metadata
```

Then collect completed outputs:

```bash
./run_pipeline.sh --software-deployment-method conda \
  --configfile config/mydata.yaml --cores 1 -- phenoradar_inputs
```

The collector is a manual target, separate from the default workflow. Run it
after the desired producer steps finish. Its log lists ready, absent,
incomplete, or disabled branches and the selected source paths at
`logs/<analysis>/phenoradar_inputs.log`.

## Published files

```text
results/<analysis>/phenoradar_inputs/
  species_metadata.tsv
  tpm.tsv                         # When OG expression is selected
  kegg/                           # When KO expression is selected
    ko_tpm_sum.tsv
    ko_modules.tsv                # Default requested grouping
    ko_pathways.tsv               # Only when requested
  alignments/{og}.faa              # When completed alignments are selected
  orthogroup_annotations.tsv[.gz]  # Optional existing OrthoDB descriptions
  species_tree.nwk                 # Optional explicitly selected Newick tree
```

No membership TSV, extra manifest, generated PhenoRadar configuration, gene–KO
copy, module expression score, or OG–KEGG join is added. Existing producer QC and
provenance remain in their original directories.

`species_metadata.tsv` contains `species`, the trait selected by `contrast.trait`
(normally `C4`), `contrast_pair_id`, and `family`. Taxonomic family comes from
the selected taxonomy metadata. Only `inputs.species_trait` supplies phenotypes;
traits from sample metadata are not used. An unavailable trait file leaves blank
traits, with a message in the metadata-step log. Missing annotations and unpaired
species stay in the dataset. PhenoRadar traits must be `0`, `1`, or blank.

## Selecting optional inputs

```yaml
phenoradar:
  orthogroups: auto
  kegg: auto
  alignments: auto
  kegg_groups: [module]
  orthogroup_annotations: null
  contrast: null
  tree: null
```

For `orthogroups`, `kegg`, and `alignments`, `auto` includes a completed branch
and reports/skips an absent or incomplete one; `true` requires completed inputs;
`false` omits it. At least one expression branch must be available. Files claimed
by a completion record must exist and pass validation; corruption is an error.

Set `kegg_groups: [module, pathway]` to include both existing KO membership maps,
or `[]` to publish KO expression alone. These fixed maps require no OG results.
For a KO-only analysis, `orthogroups: false` and `alignments: false` are sufficient
alongside a completed KEGG branch. The collector does not infer module activity
or completeness and does not fit a group model. Module aggregation, treatment
of overlapping/unassigned KOs, and sparse group lasso belong in PhenoRadar.
The current PhenoRadar expression loader accepts KO feature IDs; additional
group-model and sequence-input support require consumer-side implementation.

`orthogroup_annotations` can point to an existing headerless three-column
OG/taxid/description TSV or gzip TSV. The collector checks its format and overlap
with expressed OG IDs, reports annotation coverage, and preserves the source
file. Use the same OrthoDB release as the mappings; this target does not download
an annotation reference or infer its release from a filename.

`contrast` selects exactly one completed result branch:
`phylogeny/representatives/contrast`, `phylogeny/all/contrast`, or
`phylogeny/phenotyped/contrast`.
The collector left-joins its pair IDs onto the base metadata and validates trait
agreement, retaining species absent from the pair assignment. It does not run
pair inference. With `null`, it uses the base metadata's empty pair column.

`tree` is an explicit existing Newick path, relative to the project working
directory or absolute. Tips must match the selected metadata species exactly;
duplicate tips are rejected. Branch lengths are optional, including for an
external taxonomy tree. Choose the full molecular tree or its dated version
deliberately; a phenotyped-only or representative tree cannot cover additional
species in the metadata. The collector does not prune or infer trees.

## Expression and species checks

The collector preserves the existing four-column long expression tables,
including `run`. PhenoRadar ignores that extra column, but sums duplicate
species/feature coordinates. Consequently this export requires **one run per
species**, rejects duplicate coordinates, and never pools or averages runs.
Select compatible runs explicitly before using this target for a dataset with
multiple runs per species.

For each selected expression table, species and run identities must match the
sample manifest; values must be numeric, finite, and nonnegative. Every selected
species needs expression rows. A species with no quantified KO rows is reported
as missing instead of being silently dropped or assigned invented zeros. Missing
individual features retain the source table's semantics; consult KEGG support
and QC before selecting a PhenoRadar missing-value policy.

OG `tpm.tsv` is normalized to one million per run. KO `ko_tpm_sum.tsv` sums original
input TPM without renormalization; the default multi-KO policy can count a gene's
TPM in several KO features. To use the KO table in PhenoRadar, set its existing
data options to `feature_col: ko`, `value_col: tpm_sum`, and
`orthogroup_annotation_path: null`. Keep the dataset's other modelling settings
explicit; this export does not change them.

Protein alignments retain all gene copies and columns. OG IDs come from FASTA
filenames, and species come from the `{species}_g{number}` gene IDs. Original
headers and sequences are linked unchanged; no `species=` attribute is added.

## Species exclusions

With a nonempty top-level `exclude_species`, first run `filter_species`.
The collector uses the matching completed `results/<analysis>/filtered/`
snapshot, checks the recorded source, retained species/runs, exclusions and
selected output hashes, and publishes it at the same `phenoradar_inputs/` path.
It never creates or refreshes the filtered snapshot itself.

Filtering now prepares minimal species metadata from retained sample/taxonomy
rows and the supplied trait source. Recomputed molecular contrast branches can
be selected as above. Set `tree` to the matching pruned tree, such as
`results/<analysis>/filtered/phylogeny/all/species_tree.pruned.nwk`, when using it.
Clearing `exclude_species` selects the original dataset on the next collection.

Producer outputs use the shared [directory layout](directory_layout.md).
The published filenames above stay compact and are independent of the
producer's internal directory structure.
