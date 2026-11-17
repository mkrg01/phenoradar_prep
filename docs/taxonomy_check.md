# Taxonomy check

[Documentation](index.md) · [Species exclusion](species_filter.md)

MonoPhy compares rooted BUSCO trees with NCBI taxonomy and links conflicting
placements to RNA-seq runs. Reports support review; they never exclude samples
automatically.

## Configuration and execution

Taxonomy checks are enabled by default in `config/analysis.yaml` and run with
`--target all`. The default ranks are family, subfamily, tribe, subtribe, and
genus. To run the check on its own:

```bash
./run_analysis.sh submit --analysis results/leaf/downstream/carnivory --target taxonomy_check
```

The target includes missing tree inference. Supported tree sets are `all` and
`phenotyped`. Review jobs default to 1 CPU and 8 GB.

`outlierlevel` is the required focal-group fraction inside a candidate core
clade, in `(0, 1]`; it is not confidence. `collapse_monophyletic` affects figures only.

## Read the reports

Under `results/<build>/downstream/<analysis>/phylogeny/<set>/taxonomy_check/`:

| File | Use |
| --- | --- |
| `taxon_results.tsv` | Start here: monophyly and intruder/outlier counts by rank/group |
| `candidates.tsv` | Flagged tips, roles, and associated runs |
| `samples.tsv`, `rank_status.tsv` | Sample/rank review |
| `taxonomy.tsv`, `group_members.tsv` | Registered taxonomy and assessed membership |
| `ranks/<rank>/` | Assessment tree, native results, and PDF/SVG figures |
| `summary.json`, `engine.json`, `R_session.txt` | Settings and provenance |

An **outlier** belongs to the focal group but falls outside its selected core
clade; an **intruder** belongs to another group but falls inside. A tip may have
both roles for different groups. `focal_taxon` is the assessed group, not a
corrected identity.

Figures use triangles for intruders, squares for outliers, diamonds for both,
and gray for missing ranks. Unflagged monophyletic groups may be collapsed.

Missing ranks and single-tip groups cannot be assessed for their own monophyly.
There is no support filter or probability of mislabeling. Sampling, taxonomy,
paralogy, or tree error can cause conflicts; flags alone cannot establish
misidentification or contamination.

## References

[MonoPhy paper](https://doi.org/10.7717/peerj-cs.56) ·
[CRAN package](https://CRAN.R-project.org/package=MonoPhy) ·
[Source](https://github.com/oschwery/MonoPhy)
