# Evaluation of unmodified dating tools

Evaluation date: 2026-09-10. The preferred installation policy is to use
unmodified upstream software, with standard Conda packages where possible.
Following this evaluation, the workflow was migrated from treePL to LSD2.
The treePL observations below describe the previous implementation; current
configuration, installation and output details are in [phylogeny.md](phylogeny.md).

## Unmodified treePL

The cached official source archive for commit
`f41af04ae7cc830deadbe83a1217ed9feca60c86` was extracted afresh and built without
the project patch. Its SHA-256 was checked against the existing pinned value
`9d77514f03fa1583c1b207183e60fe11e2812218cbf83a77e37c643424d423d4`.
The build reused previously compiled bundled NLopt and ADOL-C libraries.

On synthetic balanced trees with 4 and 32 tips, heterogeneous branch lengths,
and a root fixed at 100 Ma:

- Prime, ordinary leave-one-out CV, and fixed-smoothing dating completed.
- Final trees passed topology, ultrametricity and root-calibration checks.
- Random CV returned NaN scores on the four-tip tree despite exit status zero.
- Random CV returned finite scores on the 32-tip tree. This does not eliminate
  the uninitialized-index defect visible in the source.

These observations do not establish that upstream treePL is unusable. They also
do not validate it for the full dataset. Ordinary CV avoids the particular
random-sampling failure, but other source defects, optimization settings and
large-tree cost still require assessment. The former wrapper additionally required a project-specific convergence marker
absent from upstream output. That was a wrapper requirement, rather than evidence
that unmodified upstream treePL could never be used.

## LSD2 versus ape::chronos

LSD2 2.4.4 was built with the upstream `make -j2` command and no source changes,
from cached source at commit `c61110f3a4fa05325b45c97b2134792ff9d55d4c`.
Archive SHA-256:
`9bbeaa0f8f35783c1d8dec74df6c93a804dbca808fa04484f9123de4e7258b53`.
The compiler was GNU C++ 13.3.0. Chronos used the installed conda-forge
`r-ape=5.8_1=r44h3704496_2`, with R 4.4.3.

Inputs were synthetic coalescent trees generated using `ape::rcoal`, seed
12345, rescaled to a root age of 100 Ma. Each branch was multiplied by an
independent lognormal factor with log standard deviation 0.5 and by 0.001 to
produce substitution lengths. Both methods received the same tree and
calibrations: a fixed 100 Ma root and one internal node bounded at 90–110% of
its simulated age. All tips were fixed at the present. The 32-, 128- and
512-tip input Newick files were verified to be byte-identical across methods.

Chronos used `model="correlated"`, `lambda=1` and default optimization controls.
LSD2 used its default single-rate estimation, variance weighting, fixed input
root, no outlier removal, and no confidence-interval calculation:

```bash
lsd2 -i input.nwk -d lsd-dates.txt -o dated.result \
  -s 32000 -l -1 -u 0 -U 0 -v 1
```

The site count here is an illustrative benchmark setting, not an estimate for
real CASTLES output. Internal bounds were translated to negative dates with
their order reversed. Runs used one CPU thread and a 4 GiB address-space cap.
No process hit the memory cap. Timing excludes compilation and input generation
for LSD2; Chronos process timing also includes R startup and synthetic input
generation. The measured Chronos fitting times alone were 4.28 and 6.77 seconds
for 32 and 128 tips, respectively, so startup does not explain the difference.

| Tips | LSD2 approximate wall time | Chronos wall time/status |
| --- | --- | --- |
| 32 | 0.01 s | 4.57 s; convergence flag true, but outer-iteration-limit warning |
| 128 | 0.01 s | 7.03 s; function-evaluation limit, convergence flag false |
| 512 | 0.02 s | Stopped at the 40 s trial limit |
| 5,586 | 0.72 s | Not run |

The very short LSD2 times include process-launch/measurement overhead. Its
5,586-tip run used approximately 9.1 MiB peak RSS. These are preliminary
runtime observations on synthetic inputs, not a comparison at equal accuracy
or a prediction for every empirical dataset. A 40-second timeout does not
establish that Chronos cannot finish a larger analysis.

LSD2 preserved all tips and the rooted topology at every tested size. All tip
date annotations were zero; the root annotation was -100. Calibration bounds
were satisfied within the precision of the native date output. Time-branch
root-to-tip distances differed from 100 Ma by at most about 0.000102 Ma across
these tests, consistent with its six-significant-digit output. Thus strict
downstream floating-point ultrametric checks require attention to export
precision, rather than interpreting rounded branch sums as biological
differences in tip ages.

The **time-unit output is `dated.result.date.nexus`**. The plain `.nwk` output
contains fitted substitution lengths and must not be used as a Ma tree.
The implemented workflow adapter reads native node dates and exports their
differences as time branches without changing LSD2 itself. It validates hard
calibrations within native rounding precision and retains native output and
every adjustment in the provenance.

LSD2 was selected on runtime and upstream-source grounds. Its default single estimated
rate differs from Chronos/treePL branch-specific rate models. LSD2 supports
user-specified rate partitions, but its `-q` relaxed-clock parameter controls
confidence-interval simulations rather than enabling branch-specific rate
estimation. Faster execution alone does not establish comparable dating
accuracy under broad angiosperm rate variation.

## IQ-TREE and installation choice

IQ-TREE contains LSD2 and is available through Bioconda. Its documented dating
interface takes an alignment, even when supplying a previously inferred tree.
With the actual Bioconda IQ-TREE 3.1.3 executable, this tree-only invocation:

```bash
iqtree3 -te input.nwk --date dates.txt --date-tip 0 --prefix dated -T 1
```

used IQ-TREE's date-file format (tip/date lines, comma-separated MRCA taxa,
colon-separated ranges, and no LSD2 count header). It exited successfully but
invoked phylogenetic diversity analysis and produced
no dated tree. Inspection of the official 3.1.3 source confirmed that its
phylogenetic-analysis path, including LSD2 dating, is entered when alignment or
partition input is present. This check does not rule out other future interfaces;
it establishes that this direct replacement does not date the existing tree.
The tested Conda recipe embeds LSD2 from the same `c61110f...` source used here.

No standard standalone `lsd2` package was found in Bioconda or conda-forge at the
evaluation date. The workflow therefore builds standalone LSD2 automatically in
a Snakemake Conda environment using the verified official source archive,
Conda GNU C++ 13, and the upstream makefile. There are no local LSD2 source
patches. trimAl, cdskit and nwkit use standard Conda packages. ASTRAL-IV keeps
its official 128-bit build for more than 5,000 species.

The actual three new stage environments were created successfully, including
the Snakemake post-deploy build. A real-tool workflow test with `--use-conda`
completed inference, manual-calibration dating and cached TimeTree dating;
unchanged reruns and recovery of missing dating output reused earlier inference.

The completed adapter was also run on the same 5,586-tip synthetic input using
the Conda-built LSD2 executable. All tips and the rooted topology were retained;
root-to-tip lengths ranged from 99.99999999999994 to 100.00000000000006 Ma.
The adapter and fitting together took approximately 2.66 seconds in this check.
No calibrated age needed a rounding adjustment in this case. This measures a
synthetic dating stage, not full-scale inference or accuracy on the real data.

## Evidence and references

Detailed build records, inputs, outputs, timing records and evaluation JSON
files from this session are in
`/tmp/phenoradar-dating-eval-_1zg203z/`. This directory is temporary; the findings
and methodological limits above are retained in the repository.
The final 5,586-tip adapter check is in
`/tmp/phenoradar-lsd-5586-validation-jsygphrb/summary.json`.

- [treePL official run instructions](https://github.com/blackrim/treePL/wiki/Quick-run)
- [ape Chronos reference](https://stat.ethz.ch/CRAN/web/packages/ape/refman/ape.html#chronos)
- [conda-forge ape package](https://anaconda.org/conda-forge/r-ape)
- [LSD2 source and documentation](https://github.com/tothuhien/lsd2)
- [LSD method paper](https://doi.org/10.1093/sysbio/syv068)
- [IQ-TREE dating interface](https://iqtree.github.io/doc/Dating)
- [IQ-TREE 3.1.3 entry point](https://github.com/iqtree/iqtree3/blob/v3.1.3/main/main.cpp)
- [Snakemake Conda post-deploy scripts](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#integrated-package-management)
