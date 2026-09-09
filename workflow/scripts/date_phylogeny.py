#!/usr/bin/env python3
"""Date CASTLES substitution lengths with calibrated, cross-validated treePL."""
import argparse
import json
import math
import os
from pathlib import Path
import re
import statistics
import subprocess
import tempfile

from common import atomic_writer, file_record, now, read_tsv, write_json, write_tsv
from infer_phylogeny import executable, read_tree

DEFAULT_SETTINGS = {
    "smoothing": None, "cv_start": 1000.0, "cv_stop": 0.001, "cv_multiplier": 0.1,
    "cv_replicates": 3, "replicates": 3, "numsites": None, "optimization_iterations": 2,
}


def validate_settings(settings):
    if not isinstance(settings, dict):
        raise ValueError("treePL settings must be a mapping")
    unknown = set(settings) - set(DEFAULT_SETTINGS)
    if unknown:
        raise ValueError(f"unknown treePL settings: {sorted(unknown)}")
    cfg = {**DEFAULT_SETTINGS, **settings}
    for key in ["cv_replicates", "replicates", "optimization_iterations"]:
        if type(cfg[key]) is not int or cfg[key] < 2:
            raise ValueError(f"treePL {key} must be an integer >= 2")
        if cfg[key] > 1000:
            raise ValueError(f"treePL {key} must not exceed 1000")
    for key in ["cv_start", "cv_stop", "cv_multiplier", "smoothing"]:
        value = cfg[key]
        if key == "smoothing" and value is None:
            continue
        if type(value) not in {int, float} or not math.isfinite(value) or value <= 0:
            raise ValueError(f"treePL {key} must be finite and positive")
    if not 0 < cfg["cv_multiplier"] < 1 or cfg["cv_start"] <= cfg["cv_stop"]:
        raise ValueError("treePL requires cv_start > cv_stop and 0 < cv_multiplier < 1")
    if cfg["numsites"] is not None and (type(cfg["numsites"]) is not int or not 0 < cfg["numsites"] < 2**31):
        raise ValueError("treePL numsites must be null or a positive 32-bit integer")
    if len(smoothing_grid(cfg)) > 100:
        raise ValueError("treePL smoothing grid exceeds 100 values")
    return cfg


def smoothing_grid(cfg):
    values, value = [], cfg["cv_start"]
    while value >= cfg["cv_stop"]:
        values.append(value)
        if len(values) > 100:
            break
        value *= cfg["cv_multiplier"]
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


def final_objective(path):
    objective, converged = None, None
    with Path(path).open() as handle:
        for line in handle:
            if "thorough optimization hit" in line:
                raise ValueError(f"treePL reached its iteration limit; see {path}")
            if line.startswith("treepl_final_converged: "):
                converged = line.strip().endswith(": 1")
            if line.startswith("after opt calc: "):
                objective = float(line.split(":", 1)[1])
    if converged is None or objective is None or not math.isfinite(objective) or objective >= 1e15:
        raise ValueError(f"treePL final optimization is incomplete; use the project build and inspect {path}")
    return objective, converged


def validate_dated(path, original, constraints):
    dated = read_tree(path, original.leaf_names())
    # treePL may omit/replace internal labels. Match children by descendant
    # signatures, in linear memory, and restore our stable node names.
    # Intern signatures jointly to avoid quadratic storage on ladder-like trees.
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
    names = {code: n.name for n, code in a.items()}
    for node, code in b.items():
        node.name = names[code]
    depths = {dated: 0.0}
    for node in dated.traverse():
        if not node.is_root:
            depths[node] = depths[node.up] + node.dist
    tips = [depths[n] for n in dated.leaves()]
    height = statistics.mean(tips)
    tolerance = max(1e-8, height * 1e-7)
    if not math.isfinite(height) or height <= 0 or max(tips) - min(tips) > tolerance:
        raise ValueError("treePL time tree is not ultrametric for extant tips")
    ages = {n.name: height - d for n, d in depths.items() if not n.is_leaf}
    for row in constraints:
        age = ages[row["node"]]
        if not row["min_age_ma"] - tolerance <= age <= row["max_age_ma"] + tolerance:
            raise ValueError(f"treePL output violates calibration at {row['node']}")
    return dated, ages, height


def newick_text(tree):
    from ete4.parser.newick import make_parser
    return tree.write(parser=make_parser(1, dist="%.17g"), format_root_node=True) + "\n"


def date(tree, provenance, calibrations, outdir, command, settings=None, threads=1, seed=12345):
    cfg = validate_settings({} if settings is None else settings)
    if type(threads) is not int or threads < 1 or type(seed) is not int or not 0 < seed < 2**31 - 2000:
        raise ValueError("treePL requires positive threads and a positive 32-bit seed")
    source_tree, tree = tree, read_tree(tree)
    source_qc = json.loads(Path(provenance).read_text())
    if source_qc.get("branch_length_unit") != "substitutions_per_site" or not source_qc.get("outgroup"):
        raise ValueError("dating requires a rooted tree with substitution-per-site branch lengths")
    numsites = cfg["numsites"] if cfg["numsites"] is not None else source_qc.get("total_gene_sites")
    if type(numsites) is not int or not 0 < numsites < 2**31:
        raise ValueError("treePL requires total_gene_sites in species-tree provenance (or explicit treepl.numsites)")
    if len(tree.children) != 2 or any(len(n.children) != 2 for n in tree.traverse() if not n.is_leaf):
        raise ValueError("treePL requires a rooted bifurcating species tree")
    # Restrict config tokens; paths can still contain spaces because all treePL
    # inputs use simple relative filenames in an isolated working directory.
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
    out = Path(outdir).resolve(); out.mkdir(parents=True, exist_ok=True)
    runs = out / "treepl_runs"; runs.mkdir(exist_ok=True)
    run_root = Path(tempfile.mkdtemp(prefix="run-", dir=runs))
    input_text = newick_text(tree)
    floors = [{"node": n.name, "original_substitutions_per_site": n.dist, "treepl_substitutions_per_site": 1 / numsites}
              for n in tree.traverse() if not n.is_root and n.dist < 1 / numsites]
    write_tsv(out / "branch_length_adjustments.tsv", ["node", "original_substitutions_per_site", "treepl_substitutions_per_site"], floors)
    commands = []
    def run(name, options, run_seed):
        work = run_root / name; work.mkdir()
        (work / "input.nwk").write_text(input_text)
        lines = ["treefile = input.nwk", "outfile = dated.nwk", f"numsites = {numsites}",
                 f"nthreads = {threads}", f"seed = {run_seed}",
                 *[f"{stage}iter = {cfg['optimization_iterations']}" for stage in ["lf", "pl", "cv"]]]
        for row in constraints:
            lines.extend([f"mrca = {row['node']} " + row["taxa"].replace(",", " "),
                          f"min = {row['node']} {row['min_age_ma']:.17g}", f"max = {row['node']} {row['max_age_ma']:.17g}"])
        (work / "config.txt").write_text("\n".join(lines + options) + "\n")
        argv = [command, "config.txt"]
        commands.append({"argv": argv, "cwd": str(work), "seed": run_seed})
        print(f"treePL {name}: {work}", flush=True)
        with (work / "run.log").open("w") as handle:
            subprocess.run(argv, cwd=work, stdout=handle, stderr=subprocess.STDOUT, check=True,
                           env={**os.environ, "OMP_NUM_THREADS": str(threads)})
        return work
    prime = run("prime", ["prime"], seed)
    optimizers = parse_prime(prime / "run.log")
    scores, cv_winners, means = [], [], {}
    grid = smoothing_grid(cfg)
    if cfg["smoothing"] is None:
        for replicate in range(cfg["cv_replicates"]):
            work = run(f"cv_{replicate + 1:02d}", optimizers + ["randomcv", "cvoutfile = cv.out",
                       f"cvstart = {cfg['cv_start']:.17g}", f"cvstop = {cfg['cv_stop']:.17g}",
                       f"cvmultstep = {cfg['cv_multiplier']:.17g}"], seed + replicate + 1)
            values = parse_cv(work / "cv.out", grid)
            scores.extend({"replicate": replicate + 1, "smoothing": grid[i], "chisq": score}
                          for i, (_, score) in enumerate(values))
            cv_winners.append(min(values, key=lambda x: x[1])[0])
        means = {s: statistics.mean(r["chisq"] for r in scores if r["smoothing"] == s) for s in grid}
        smoothing = min(means, key=means.get)
    else:
        smoothing = cfg["smoothing"]
    write_tsv(out / "cross_validation.tsv", ["replicate", "smoothing", "chisq"], scores)
    results = []
    for replicate in range(cfg["replicates"]):
        work = run(f"final_{replicate + 1:02d}", optimizers + [f"smooth = {smoothing:.17g}"], seed + 101 + replicate)
        objective, native_converged = final_objective(work / "run.log")
        dated, ages, height = validate_dated(work / "dated.nwk", tree, constraints)
        results.append({"replicate": replicate + 1, "objective": objective, "root_age_ma": height,
                        "native_converged": native_converged,
                        "work": str(work), "ages": ages})
    best = min(results, key=lambda r: r["objective"])
    work = Path(best["work"])
    dated, ages, height = validate_dated(work / "dated.nwk", tree, constraints)
    age_ranges = [{"node": node, "age_ma": age,
                   "replicate_min_age_ma": min(r["ages"][node] for r in results),
                   "replicate_max_age_ma": max(r["ages"][node] for r in results)} for node, age in ages.items()]
    max_spread = max(r["replicate_max_age_ma"] - r["replicate_min_age_ma"] for r in age_ranges) / height
    diagnostics = []
    objective_spread = (max(r["objective"] for r in results) - best["objective"]) / max(1.0, abs(best["objective"]))
    if objective_spread > 1e-4:
        diagnostics.append("optimization objectives differ across replicates; increase optimization_iterations and inspect results")
    if cfg["smoothing"] is None and smoothing in (grid[0], grid[-1]):
        diagnostics.append("CV optimum is on the grid boundary; extend the grid before scientific interpretation")
    if len(set(cv_winners)) > 1:
        diagnostics.append("CV replicates select different smoothing values; inspect cross_validation.tsv")
    if max_spread > 0.05:
        diagnostics.append("node ages vary by more than 5% of root age across optimization replicates")
    for filename, content in [("species_tree.dated.nwk", newick_text(dated)),
                              ("treepl.input.nwk", input_text),
                              ("treepl.config.txt", (work / "config.txt").read_text()
                               .replace("treefile = input.nwk", "treefile = treepl.input.nwk")
                               .replace("outfile = dated.nwk", "outfile = treepl.dated.nwk")),
                              ("treepl.dated.nwk", (work / "dated.nwk").read_text())]:
        with atomic_writer(out / filename) as handle:
            handle.write(content)
    write_tsv(out / "calibrations.resolved.tsv", list(constraints[0]), constraints)
    write_tsv(out / "node_ages.tsv", list(age_ranges[0]), age_ranges)
    fields = ["replicate", "objective", "root_age_ma", "native_converged", "work"]
    write_tsv(out / "optimization_replicates.tsv", fields, [{k: r[k] for k in fields} for r in results])
    write_json(out / "provenance.json", {"created_at": now(), "commands": commands,
               "executable": file_record(command), "tree": file_record(source_tree),
               "source_provenance": file_record(provenance), "calibrations": file_record(calibrations),
               "branch_length_unit": "million_years", "method": "treePL penalized likelihood",
               "rate_model": "branch-specific rates with additive smoothing penalty", "confidence_intervals": False,
               "settings": cfg, "numsites": numsites,
               "numsites_source": "explicit override" if cfg["numsites"] is not None else "sum of retained trimmed gene sites",
               "numsites_is_effective_sample_size": False, "smoothing": smoothing,
               "cv_method": "random subsample and replicate (10 groups per run)" if scores else "fixed smoothing",
               "cv_replicate_winners": cv_winners, "selected_replicate": best["replicate"],
               "cv_mean_scores": [{"smoothing": s, "mean_chisq": score} for s, score in means.items()],
               "selected_workdir": str(work), "max_node_age_spread_fraction_of_root": max_spread,
               "relative_objective_spread": objective_spread,
               "optimization_validation": "independent restart stability; global convergence is not established",
               "short_branches_adjusted": len(floors), "diagnostics": diagnostics,
               "review_status": "needs_review" if diagnostics else "checks_passed",
               "concatenation_used": False, "root_age_ma": height})
    for message in diagnostics:
        print("treePL diagnostic: " + message, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["tree", "provenance", "calibrations", "outdir", "command"]:
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--settings", type=json.loads, default={})
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--seed", type=int, default=12345)
    date(**vars(parser.parse_args()))
