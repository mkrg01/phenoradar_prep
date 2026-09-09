# Pinned treePL build

`prepare_phylogeny_tools.py` verifies upstream commit
[`f41af04ae7cc830deadbe83a1217ed9feca60c86`](https://github.com/blackrim/treePL/tree/f41af04ae7cc830deadbe83a1217ed9feca60c86)
and applies [treepl.patch](treepl.patch) before compilation. The patch and binary
checksums are part of `resources/phylogeny_tools/treepl.json`. The source archive
includes NLopt 2.4.2 and ADOL-C 2.6.3; these are built with the same recipe and
linked statically. ADOL-C is built with OpenMP, without Boost/ColPack.

The patch fixes the following observed source defects:

- `check_possible_cv_nodes` indexes `parent_nds_ints` using the uninitialized
  variable `cvpar`; use its actual argument `cvnode`. This preserves the intended
  exclusion of sister tips within a random CV sample.
- Random CV rounds a 10% sample down to zero for four-tip trees; use at least
  one tip. Larger inputs retain the original sample size and ten groups.
- Initialize the three `prime` detail flags to false. Otherwise a search with no
  improving candidate can print indeterminate flags.
- Pass `bestcvopt` and `moredetailcvad` to CV optimization, as configured by
  `optcvad`, rather than reusing the non-CV autodiff optimizer settings.
- Return `*this` from three reference-returning autodiff compound operators
  (`-=`, `*=`, `/=`), which previously reached the end without returning.
- Apply output precision to the Newick string streams themselves. Upstream's
  unattached `setprecision(10)` did not change their six-decimal default. Use 17
  digits in native output and objective logs to avoid damaging short branches
  and to compare repeated final objectives.
- Print `treepl_final_converged` from the existing optimizer's return value.
  This adds observability; it does not change the optimization or stopping rule.

No likelihood, rate penalty, calibration constraint or optimization stopping
criterion is changed. The wrapper deliberately uses bounded native iterations
and independent restarts, rather than `thorough`; see [dating documentation](../../docs/phylogeny.md).
The native flag tests progress within a perturbation/optimization cycle and is
recorded separately from agreement between independent restarts. Neither is a
proof of a global optimum. Synthetic tests check time units, constraints,
random CV (including samples containing multiple tips), and workflow integration.
They do not establish empirical dating accuracy.
