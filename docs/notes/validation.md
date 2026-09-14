# Recorded dataset checks

[Documentation](../index.md) · [Automated tests](../development.md)

These notes preserve checks reported during development in September 2026.
They are historical observations, not a report of the current dataset status.
The referenced result directories are local artifacts and are not distributed
with the repository.

## BUSCO sequence preparation and trimming

A separate real-input smoke check used Amborella trichopoda, Oryza sativa,
Abelia chinensis and Abeliophyllum distichum from tlight, with 20 shared Complete
BUSCO markers (sorted by ID; each match at least 100 aa). All 80 CDS passed
preparation, with no reading-frame changes; terminal stops were masked. FAMSA
and the previous trimAl/QC retained 18 markers and excluded two with eight
parsimony-informative sites each. With the simplified QC, all 20 are retained:
the two recovered alignments have 189/150 sites, including 134/113 variable
sites, respectively. The original 18 retained alignments and column maps are
unchanged. Original results are in `results/phylogeny_cdskit_trimal_smoke/summary.json`;
the repeated trimAl/QC and counts are in
`results/phylogeny_qc_simplification_smoke/summary.json`. The smoke check's
preselection by match length is not a rule in the production marker selection.
This checked input compatibility through trimming; it did not infer a species tree or assess
full-scale accuracy/performance.

## Contrast representatives

In the recorded dataset, the first skim selects 199 representatives from 2,047
known-trait species. Root selection chooses `Nymphaea_colorata` within those
199 representatives, while the full 5,586-species manifest selects
`Amborella_trichopoda`. The full 199-species contrast inference had not been run
at the time of the check.
A small real-data check selects eight known-trait representatives from nine
input species, roots on `Nymphaea_colorata`, and completes inference with 20
markers, three contrast pairs, and both figure formats. This checks execution
and membership; it does not establish the accuracy of the inferred phylogeny.

## TimeTree and scale

A live TimeTree check made three queries for a four-species test tree and
retrieved calibration candidates and study information. The topology and branch
lengths were supplied for testing, rather than inferred from BUSCO.

The method review reported that inference had not yet been run for all 5,586
species. These small checks did not assess full-scale inference resources,
biological accuracy, calibration sensitivity, or downstream model performance.
The separate [dating-tool evaluation](dating_evaluation.md) includes synthetic
trees at that tip count; those timings cover dating, not BUSCO tree inference.
