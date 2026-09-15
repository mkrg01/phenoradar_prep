#!/usr/bin/env python3
"""Date CASTLES substitution lengths with unmodified upstream treePL and validated native cross-validation."""
import argparse
import json
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile

from common import atomic_writer, file_record, now, read_tsv, write_json, write_tsv
from infer_phylogeny import executable, read_tree

DEFAULT_SETTINGS = {
    # Native option names and defaults from the pinned upstream source.
    # None selects native CV instead of passing a fixed smooth value.
    "smooth": None, "cvstart": 1000.0, "cvstop": 0.1, "cvmultstep": 0.1,
    "lfiter": 3, "pliter": 5, "cviter": 3, "thorough": False,
}
TREEPL_COMMIT = "f41af04ae7cc830deadbe83a1217ed9feca60c86"
TREEPL_ARCHIVE_SHA256 = "9d77514f03fa1583c1b207183e60fe11e2812218cbf83a77e37c643424d423d4"


def validate_settings(settings):
    if not isinstance(settings, dict):
        raise ValueError("treePL settings must be a mapping")
    unknown = set(settings) - set(DEFAULT_SETTINGS)
    if unknown:
        raise ValueError(f"unknown treePL settings: {sorted(unknown)}")
    cfg = {**DEFAULT_SETTINGS, **settings}
    if type(cfg["thorough"]) is not bool:
        raise ValueError("treePL thorough must be a boolean")
    for key in ["lfiter", "pliter", "cviter"]:
        if type(cfg[key]) is not int or not 0 < cfg[key] < 2**31:
            raise ValueError(f"treePL {key} must be a positive signed 32-bit integer")
    for key in ["cvstart", "cvstop", "cvmultstep", "smooth"]:
        value = cfg[key]
        if key == "smooth" and value is None:
            continue
        if type(value) not in {int, float} or not math.isfinite(value) or value <= 0:
            raise ValueError(f"treePL {key} must be finite and positive")
    if not 0 < cfg["cvmultstep"] < 1 or cfg["cvstart"] < cfg["cvstop"]:
        raise ValueError("treePL requires cvstart >= cvstop and 0 < cvmultstep < 1")
    if len(smoothing_grid(cfg)) > 100:
        raise ValueError("treePL smoothing grid exceeds 100 values")
    return cfg


def smoothing_grid(cfg):
    values, value = [], cfg["cvstart"]
    while value >= cfg["cvstop"]:
        values.append(value)
        if len(values) > 100:
            break
        value *= cfg["cvmultstep"]
    return values


def calibration_rows(tree, path):
    names = set(tree.leaf_names())
    constraints, used = [], set()
    for row in read_tsv(path):
        if not {"taxa", "min_age_ma", "max_age_ma", "source"} <= row.keys():
            raise ValueError("calibrations require taxa, min_age_ma, max_age_ma, source columns")
        taxa = row["taxa"].split(",")
        if len(set(taxa)) != len(taxa) or len(taxa) < 2 or set(taxa) - names:
            raise ValueError(f"calibration taxa must name at least two distinct tree tips: {row['taxa']}")
        minimum, maximum = float(row["min_age_ma"]), float(row["max_age_ma"])
        if not (math.isfinite(minimum) and math.isfinite(maximum) and 0 < minimum <= maximum):
            raise ValueError("calibration ages must be finite positive Ma bounds with min <= max")
        if not row["source"].strip():
            raise ValueError("every calibration must record its source")
        node = tree.common_ancestor(taxa)
        if node.name in used:
            raise ValueError(f"multiple calibrations resolve to the same MRCA: {row['taxa']}")
        used.add(node.name)
        constraints.append({**row, "node": node.name, "min_age_ma": minimum, "max_age_ma": maximum})
    if not constraints:
        raise ValueError("absolute dating requires at least one explicit calibration; no arbitrary root age is supplied")
    bounds = {c["node"]: c for c in constraints}
    for node in tree.traverse("postorder"):
        lower = max((child.props["age_lower"] for child in node.children), default=0)
        bound = bounds.get(node.name)
        if bound:
            lower = max(lower, bound["min_age_ma"])
            if lower > bound["max_age_ma"]:
                raise ValueError(f"inconsistent ancestor/descendant calibration bounds at {node.name}")
        node.add_prop("age_lower", lower)
    return constraints



def parse_prime(path):
    options, found = [], set()
    with Path(path).open() as handle:
        for line in handle:
            line = line.strip()
            match = re.fullmatch(r"(opt|optad|optcvad) = ([0-5])", line)
            if match:
                if match[1] in found:
                    raise ValueError("treePL prime produced duplicate optimizer settings")
                found.add(match[1]); options.append(line)
            elif line in {"moredetail", "moredetailad", "moredetailcvad"}:
                options.append(line)
    if found != {"opt", "optad", "optcvad"}:
        raise ValueError(f"treePL prime did not report all optimizers; see {path}")
    return options


def parse_cv(path, grid):
    rows = []
    for line in Path(path).read_text().splitlines():
        match = re.fullmatch(r"chisq: \(([^)]+)\) (\S+)", line.strip())
        if not match:
            raise ValueError(f"malformed treePL CV result: {line}")
        smooth, score = map(float, match.groups())
        if not all(math.isfinite(v) and v >= 0 for v in (smooth, score)):
            raise ValueError("treePL CV contains nonfinite/negative scores")
        rows.append((smooth, score))
    if len(rows) != len(grid) or any(not math.isclose(s, g, rel_tol=1e-6) for (s, _), g in zip(rows, grid)):
        raise ValueError("treePL CV did not complete the requested smoothing grid")
    return rows


def final_objective(path, expected_fits=1):
    values = []
    for line in Path(path).read_text().splitlines():
        if "thorough optimization hit" in line:
            raise ValueError(f"treePL reached its iteration limit; see {path}")
        if line.startswith("after opt calc: "):
            values.append(float(line.split(":", 1)[1]))
    if (len(values) != expected_fits
            or any(not math.isfinite(value) or value >= 1e15 for value in values)):
        raise ValueError(f"treePL final optimization is incomplete; see {path}")
    return values[-1]


def final_smoothing(path, expected_fits=1):
    values = [float(line.split(":", 1)[1]) for line in Path(path).read_text().splitlines()
              if line.startswith("smoothing:")]
    if len(values) != expected_fits or any(not math.isfinite(v) or v <= 0 for v in values):
        raise ValueError(f"treePL final smoothing is missing or invalid; see {path}")
    return values[-1]


def validate_dated(path, original, constraints):
    """Validate native six-decimal lengths without altering any branch length."""
    dated = read_tree(path, original.leaf_names())
    dated.dist = 0.0  # The root has no parent edge.
    intern = {}
    def codes(tree):
        result = {}
        for node in tree.traverse("postorder"):
            key = ("tip", node.name) if node.is_leaf else ("node", tuple(sorted(result[c] for c in node.children)))
            result[node] = intern.setdefault(key, len(intern))
        return result
    a, b = codes(original), codes(dated)
    if a[original] != b[dated] or len(a) != len(b):
        raise ValueError("treePL changed the rooted topology")
    names = {code: node.name for node, code in a.items()}
    for node, code in b.items():
        node.name = names[code]
    bounds = {row["node"]: row for row in constraints}
    lower, upper, counts, ages = {}, {}, {}, {}
    for node in dated.traverse("postorder"):
        if node.is_leaf:
            lo = hi = age = 0.0
            counts[node] = 1
        else:
            # Upstream prints each edge to six decimal places. Check whether
            # those intervals admit an ultrametric tree and the supplied bounds;
            # this is validation only, with no branch-length reconciliation.
            tolerance = {c: 0.500001e-6 + 2 * math.ulp(c.dist) for c in node.children}
            lo = max(lower[c] + max(0.0, c.dist - tolerance[c]) for c in node.children)
            hi = min(upper[c] + c.dist + tolerance[c] for c in node.children)
            counts[node] = sum(counts[c] for c in node.children)
            age = sum((ages[c] + c.dist) * counts[c] for c in node.children) / counts[node]
        if node.name in bounds:
            lo = max(lo, bounds[node.name]["min_age_ma"])
            hi = min(hi, bounds[node.name]["max_age_ma"])
        if lo > hi + 8 * math.ulp(max(1.0, abs(lo), abs(hi))):
            raise ValueError(f"treePL output violates calibration or ultrametricity beyond native rounding at {node.name}")
        lower[node], upper[node], ages[node] = lo, hi, age
    height = ages[dated]
    if not math.isfinite(height) or height <= 0:
        raise ValueError("treePL returned an invalid time-tree height")
    return dated, {n.name: age for n, age in ages.items() if not n.is_leaf}, height


def newick_text(tree):
    from ete4.parser.newick import make_parser
    return tree.write(parser=make_parser(1, dist="%.17g"), format_root_node=True) + "\n"


def treepl_build(command):
    """Accept only the checked, unmodified build installed by dating.post-deploy.sh."""
    binary = file_record(command)
    manifest = Path(command).resolve().parent.parent / "share/treepl/build.json"
    if not manifest.is_file():
        raise ValueError("treePL build provenance missing; run workflow/envs/dating.post-deploy.sh")
    build = json.loads(manifest.read_text())
    if (build.get("source_patch_applied") is not False or build.get("commit") != TREEPL_COMMIT
            or build.get("archive_sha256") != TREEPL_ARCHIVE_SHA256
            or build.get("executable_sha256") != binary["sha256"]):
        raise ValueError("treePL executable is not the verified unmodified build; recreate the dating environment")
    return build, binary


def date(tree, provenance, calibrations, outdir, command, settings=None, threads=1, seed=12345):
    cfg = validate_settings({} if settings is None else settings)
    if type(threads) is not int or threads != 1:
        raise ValueError("the validated treePL workflow requires threads=1")
    if type(seed) is not int or not 0 < seed < 2**31:
        raise ValueError("treePL requires a positive signed 32-bit seed")
    source_tree, tree = tree, read_tree(tree)
    source_qc = json.loads(Path(provenance).read_text())
    if source_qc.get("branch_length_unit") != "substitutions_per_site" or not source_qc.get("outgroup"):
        raise ValueError("dating requires a rooted tree with substitution-per-site branch lengths")
    if source_qc["outgroup"] not in set(tree.leaf_names()):
        raise ValueError("dating outgroup is absent from the species tree")
    numsites = source_qc.get("total_gene_sites")
    if type(numsites) is not int or not 0 < numsites < 2**31:
        raise ValueError("treePL requires total_gene_sites as a positive 32-bit integer in species-tree provenance")
    if any(len(n.children) != 2 for n in tree.traverse() if not n.is_leaf):
        raise ValueError("treePL requires a rooted bifurcating species tree")
    if not any(n.dist > 0 for n in tree.traverse() if not n.is_root):
        raise ValueError("dating requires positive substitution lengths")
    if any(not re.fullmatch(r"[A-Za-z0-9_.-]+", n) for n in tree.leaf_names()):
        raise ValueError("treePL tip labels must contain only letters, digits, underscores, dots or hyphens")
    used = set(tree.leaf_names())
    for number, node in enumerate(tree.traverse(), 1):
        if not node.is_leaf:
            node.name = f"N{number:06d}"
            while node.name in used:
                node.name = "N" + node.name
            used.add(node.name)
    constraints = calibration_rows(tree, calibrations)
    command = executable(command)
    build, binary = treepl_build(command)
    out = Path(outdir).resolve()
    runs = out / "treepl_runs"
    runs.mkdir(parents=True, exist_ok=True)
    run_root = Path(tempfile.mkdtemp(prefix="run-", dir=runs))
    input_text = newick_text(tree)
    # Report upstream's native floor; do not change the input tree ourselves.
    floors = [{"node": n.name, "original_substitutions_per_site": n.dist,
               "treepl_substitutions_per_site": 1 / numsites}
              for n in tree.traverse() if not n.is_root and n.dist < 1 / numsites]
    commands = []
    def run(name, options):
        work = run_root / name
        work.mkdir()
        (work / "input.nwk").write_text(input_text)
        lines = ["treefile = input.nwk", "outfile = dated.nwk", f"numsites = {numsites}",
                 "nthreads = 1", f"seed = {seed}",
                 *(["thorough"] if cfg["thorough"] else []),
                 *[f"{key} = {cfg[key]}" for key in ("lfiter", "pliter", "cviter")]]
        for row in constraints:
            lines.extend([f"mrca = {row['node']} " + row["taxa"].replace(",", " "),
                          f"min = {row['node']} {row['min_age_ma']:.17g}",
                          f"max = {row['node']} {row['max_age_ma']:.17g}"])
        (work / "config.txt").write_text("\n".join(lines + options) + "\n")
        invocation = {"argv": [command, "config.txt"], "cwd": str(work), "seed": seed}
        commands.append(invocation)
        write_json(work / "command.json", invocation)
        print(f"treePL {name}: {work}", flush=True)
        with (work / "run.log").open("w") as log:
            subprocess.run(invocation["argv"], cwd=work, stdout=log, stderr=subprocess.STDOUT,
                           check=True, env={**os.environ, "OMP_NUM_THREADS": "1"})
        if "thorough optimization hit" in (work / "run.log").read_text():
            raise ValueError(f"treePL reached its iteration limit; see {work / 'run.log'}")
        return work
    prime = run("prime", ["prime"])
    optimizers = parse_prime(prime / "run.log")
    grid = smoothing_grid(cfg) if cfg["smooth"] is None else []
    if grid:
        # Native CV chooses smoothing and performs the final fit in this same
        # process. No external repetitions, score averaging, or fit selection.
        # Leave-one-out CV avoids the pinned source's randomcv sampling defects.
        options = ["cv", "cvoutfile = cv.out", f"cvstart = {cfg['cvstart']:.17g}",
                   f"cvstop = {cfg['cvstop']:.17g}", f"cvmultstep = {cfg['cvmultstep']:.17g}"]
    else:
        options = [f"smooth = {cfg['smooth']:.17g}"]
    work = run("fit", optimizers + options)
    values = parse_cv(work / "cv.out", grid) if grid else []
    expected_fits = len(grid) + 1
    objective = final_objective(work / "run.log", expected_fits)
    smoothing = final_smoothing(work / "run.log", expected_fits)
    allowed = grid or [cfg["smooth"]]
    if not any(math.isclose(smoothing, s, rel_tol=1e-7) for s in allowed):
        raise ValueError("treePL final smoothing is outside the requested values")
    dated, ages, height = validate_dated(work / "dated.nwk", tree, constraints)
    diagnostics = []
    if grid and any(math.isclose(smoothing, s, rel_tol=1e-7) for s in (grid[0], grid[-1])):
        diagnostics.append("CV optimum is on the grid boundary; extend the grid before scientific interpretation")
    scores = [{"smoothing": grid[i], "chisq": score} for i, (_, score) in enumerate(values)]
    def exported_config(path):
        return (path.read_text().replace("treefile = input.nwk", "treefile = treepl.input.nwk")
                .replace("outfile = dated.nwk", "outfile = treepl.dated.nwk")
                .replace("cvoutfile = cv.out", "cvoutfile = treepl.cv.out"))
    # Only restore node labels for the pipeline tree; native edge values stay intact.
    # Publish after validation and write provenance last.
    for filename, content in [("species_tree.dated.nwk", newick_text(dated)),
                              ("treepl.input.nwk", input_text),
                              ("treepl.prime.config.txt", exported_config(prime / "config.txt")),
                              ("treepl.prime.log", (prime / "run.log").read_text()),
                              ("treepl.config.txt", exported_config(work / "config.txt")),
                              ("treepl.cv.out", (work / "cv.out").read_text() if grid else ""),
                              ("treepl.dated.nwk", (work / "dated.nwk").read_text()),
                              ("treepl.log", (work / "run.log").read_text())]:
        with atomic_writer(out / filename) as handle:
            handle.write(content)
    write_tsv(out / "calibrations.resolved.tsv", list(constraints[0]), constraints)
    write_tsv(out / "node_ages.tsv", ["node", "age_ma"],
              [{"node": node, "age_ma": age} for node, age in ages.items()])
    write_tsv(out / "cross_validation.tsv", ["smoothing", "chisq"], scores)
    write_tsv(out / "branch_length_adjustments.tsv", ["node", "original_substitutions_per_site", "treepl_substitutions_per_site"], floors)
    # Retired reports must not be exported alongside the current fit.
    for obsolete in ("optimization_replicates.tsv", "rounding_adjustments.tsv"):
        (out / obsolete).unlink(missing_ok=True)
    write_json(out / "provenance.json", {
        "created_at": now(), "commands": commands, "executable": binary, "build": build,
        "tree": file_record(source_tree), "source_provenance": file_record(provenance),
        "calibrations": file_record(calibrations), "branch_length_unit": "million_years",
        "method": "treePL penalized likelihood", "rate_model": "branch-specific rates with additive smoothing penalty",
        "settings": cfg, "optimizer_options": optimizers,
        "numsites": numsites, "numsites_source": "sum of retained trimmed gene sites",
        "numsites_is_effective_sample_size": False, "confidence_intervals": False,
        "root_preserved": True, "topology_preserved": True, "tip_age_ma": 0,
        "concatenation_used": False, "root_age_ma": height, "smoothing": smoothing,
        "cv_method": "native leave-one-out" if grid else "fixed smoothing",
        "smoothing_selection": "treePL native CV" if grid else "user-specified smooth",
        "objective": objective, "workdir": str(work),
        "short_branches_adjusted": len(floors), "native_time_branch_lengths_preserved": True,
        "node_age_calculation": "mean descendant-tip distance from native rounded branch lengths",
        "optimization_validation": "native completion and output checks; global convergence is not established",
        "diagnostics": diagnostics, "review_status": "needs_review" if diagnostics else "checks_passed"})
    for message in diagnostics:
        print("treePL diagnostic: " + message, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("tree", "provenance", "calibrations", "outdir", "command"):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--settings", type=json.loads, default={})
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--seed", type=int, default=12345)
    date(**vars(parser.parse_args()))
