# Dating angiosperm species trees

[Documentation](index.md) · [BUSCO phylogeny](phylogeny.md)

Dating currently supports angiosperms only. It applies [AngioCal v1.0](https://doi.org/10.5281/zenodo.3828071)
fossil minima to the existing rooted tree with [treePL](https://doi.org/10.1093/bioinformatics/bts492).
Each selected tree produces one dated tree in millions of years (Ma). treePL uses
smoothing 10 and a root maximum of 154 Ma, following the young-tree settings in
[Zuntini et al. (2024)](https://doi.org/10.1038/s41586-024-07324-0).

Enable dating with `phylogeny.dating.enabled: true` in `analysis.yaml` and select
`phylogeny.trees: [phenotyped]` or `[all]`. Both are supported if requested together.
The input requires the selected tips' NCBI TaxIDs and the frozen taxonomy database.

```bash
./run_analysis.sh submit --analysis results/leaf/downstream/c4 --target phylogeny_calibrations
./run_analysis.sh submit --analysis results/leaf/downstream/c4 --target timetree
```

The first target prepares the fossil report for review. The second also prepares
missing placements, fits the dated tree, and plots it. Enabled dating runs with
`all` too.

NWKIT places stem fossils where the sampled clade and adjacent lineages agree
with the frozen NCBI taxonomy. When the closest sister is absent, the minimum
can apply to an older sampled ancestor. Crown fossils require reviewed anchors.
To include them or override a stem placement, set `inputs.fossil_calibration_map`
to a TSV with columns `fossil_id`, `left_species`, and `right_species`. The
anchors identify the calibrated MRCA in each requested tree. Use distinct species
representing the two sides of the intended divergence. The workflow freezes this
file when preparing the analysis.

In `phylogeny/<set>/angiocal/`, inspect `fossils.tsv` for placements, exclusions
and references, and `calibrations.tsv` for accepted minima. When multiple fossils
map to a node, its strongest minimum is used and all records remain in the report.
Dating requires at least one applicable fossil minimum.

In `phylogeny/<set>/dating/`, `species_tree.dated.nwk` and `node_ages.tsv`
contain the result. `treepl.conf`, `treepl.log` and `provenance.json` record the
fit, checks, root age and whether it reached the upper bound. An estimate at 154
Ma is constrained by that upper bound. The dates are conditional on the input
topology, fossil placements, branch lengths, and smoothing.

The workflow then writes `species_tree.dated.pdf` and `species_tree.dated.svg`
as vector figures and `species_tree.dated.png` at 600 dpi in the same directory.
The plot is 178 mm wide and shows the time axis, the configured trait, and
NCBI families. Each family bar spans consecutive tips assigned to that family.
The plot labels family bars containing at least one tip with trait state `1`.
Leader lines connect the labels to their bars. If space is limited, the plot
prioritizes longer bars.
`species_tree.dated.tips.tsv` lists tips from top to bottom with species names,
families, and original trait labels. The plot uses the estimated node ages and
requests one CPU and 1 GB of memory.

Placement and dating each use one CPU and request 4 GB by default. The dating
environment installs checksum-verified copies of treePL, NWKIT, and the AngioCal
workbook once. Jobs reuse the frozen taxonomy snapshot and the environment's copy
of the workbook, so neither step requires external requests.

Use the original release for an existing frozen analysis. Prepare a new
analysis with this configuration to use AngioCal and treePL.
