# Dating angiosperm species trees

[Documentation](index.md) · [BUSCO phylogeny](phylogeny.md)

The dating branch fits divergence ages to a rooted, bifurcating angiosperm
tree using [AngioCal v1.0](https://zenodo.org/records/3828071) fossil minima
and [treePL](https://academic.oup.com/bioinformatics/article/28/20/2689/203074).
The input tree has branch lengths in substitutions per site. The dated tree
reports branch lengths in millions of years (Ma).

## Configuration and execution

The supplied analysis config enables dating for the `all` tree. Dating
supports `all`, `phenotyped`, or both under `phylogeny.trees` and uses the
frozen NCBI TaxIDs and taxonomy. Set `phylogeny.dating.enabled: false` before
preparing an analysis with non-angiosperms or `representatives`.

```bash
./run_analysis.sh submit --analysis results/leaf/downstream/c4 --target phylogeny_calibrations
./run_analysis.sh submit --analysis results/leaf/downstream/c4 --target timetree
```

Run the first target to review fossil placements before fitting ages. The
`timetree` target includes missing prerequisites. Dating also runs with
`all` when enabled.

## Review fossil placements

NWKIT maps stem fossils where the sampled clade and adjacent lineages agree
with the frozen NCBI taxonomy. An unsampled sister lineage can move the
minimum to an older sampled ancestor. Crown fossils require reviewed anchors.
To add a crown placement or override a stem placement, set
`inputs.fossil_calibration_map` to a TSV with `fossil_id`, `left_species`,
and `right_species`. Use two species spanning the intended divergence on each
selected tree. The file is frozen when the analysis is prepared.

Check `phylogeny/<set>/angiocal/fossils.tsv` for placements, exclusions, and
sources. The dating fit requires at least one applicable fossil minimum.
When several fossils calibrate the same node, the oldest minimum is used and
the full evidence is kept in the report.

## Interpret the dated tree

treePL fits one tree with a fixed smoothing value of 10 and a maximum age of
154 Ma at the root of the **sampled** tree. Zuntini et al. tested multiple
smoothing values and used 154 Ma as an upper bound on crown angiosperms in
their young-tree scenario. Check which clade your sampled root represents
before interpreting its age. [Zuntini et al. (2024)](https://doi.org/10.1038/s41586-024-07324-0)

In `phylogeny/<set>/dating/`, use `species_tree.dated.nwk` and
`node_ages.tsv` for ages. Check `provenance.json` and `treepl.log` for the fit
and whether the root reached its bound. Ages are conditional on the input
tree, calibrations, and smoothing. Figures (`PDF`, `SVG`, and `PNG`) show the
time axis, configured trait, and families. The adjacent `*.tips.tsv` lists tips
in figure order.
