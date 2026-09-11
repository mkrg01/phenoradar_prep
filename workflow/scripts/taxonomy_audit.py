#!/usr/bin/env python3
"""Run MonoPhy on an existing species tree and link its review flags to runs."""
import argparse
from collections import Counter, defaultdict
from functools import lru_cache
import json
import math
import os
from pathlib import Path
import shutil
import sqlite3
import subprocess
import tempfile

from common import file_record, now, read_tsv, write_json, write_tsv

DEFAULTS = {"enabled": False, "ranks": ["family", "subfamily", "tribe", "subtribe", "genus"],
            "outlierlevel": 0.5, "collapse_monophyletic": True, "mem_gb": 8}
RANKS = set("superkingdom kingdom subkingdom superphylum phylum subphylum superclass class subclass "
            "infraclass superorder order suborder infraorder parvorder superfamily family subfamily "
            "tribe subtribe genus subgenus section subsection series subseries species subspecies "
            "varietas forma".split())
MONOPHY_VERSION = "1.3.2"


def validate_settings(settings):
    if not isinstance(settings, dict):
        raise ValueError("taxonomy_audit must be a mapping")
    removed = set(settings) & {"min_reference_species", "max_plot_species"}
    if removed:
        raise ValueError("removed taxonomy_audit settings; delete " + ", ".join(sorted(removed)) +
                         "; detection now uses MonoPhy and plots have no tip limit")
    unknown = set(settings) - set(DEFAULTS)
    if unknown:
        raise ValueError("unknown taxonomy_audit settings: " + ", ".join(sorted(unknown)))
    result = {**DEFAULTS, **settings}
    ranks = result["ranks"]
    if (not isinstance(ranks, list) or not ranks or
            any(not isinstance(r, str) or r not in RANKS for r in ranks) or len(set(ranks)) != len(ranks)):
        raise ValueError("taxonomy_audit.ranks must be a nonempty list of unique named NCBI ranks")
    for key in ["enabled", "collapse_monophyletic"]:
        if type(result[key]) is not bool:
            raise ValueError(f"taxonomy_audit.{key} must be true or false")
    threshold = result["outlierlevel"]
    if type(threshold) not in (int, float) or not math.isfinite(threshold) or not 0 < threshold <= 1:
        raise ValueError("taxonomy_audit.outlierlevel must be a number in (0, 1]")
    if type(result["mem_gb"]) is not int or result["mem_gb"] < 1:
        raise ValueError("taxonomy_audit.mem_gb must be an integer >= 1")
    return result


class Taxonomy:
    """Read only the project's frozen ETE database; no NCBI/network calls."""
    def __init__(self, path):
        self.db = sqlite3.connect(Path(path).resolve().as_uri() + "?mode=ro", uri=True)

    @lru_cache(maxsize=None)
    def canonical(self, taxid):
        seen = set()
        while taxid not in seen:
            seen.add(taxid)
            row = self.db.execute("SELECT taxid_new FROM merged WHERE taxid_old=?", (taxid,)).fetchone()
            if row is None:
                return taxid
            taxid = row[0]
        raise ValueError("cycle in merged taxonomy IDs")

    @lru_cache(maxsize=None)
    def node(self, taxid):
        return self.db.execute("SELECT taxid,spname,rank,track FROM species WHERE taxid=?",
                               (self.canonical(taxid),)).fetchone()

    @lru_cache(maxsize=None)
    def lineage(self, taxid):
        node = self.node(taxid)
        if node is None:
            return {}
        result = {}
        for value in node[3].split(","):
            ancestor = self.node(int(value))
            if ancestor is None:
                raise ValueError(f"incomplete frozen taxonomy lineage for {taxid}")
            if ancestor[2] in RANKS:
                if ancestor[2] in result:
                    raise ValueError(f"repeated rank in frozen taxonomy lineage for {taxid}")
                result[ancestor[2]] = (ancestor[0], ancestor[1])
        return result


def load_samples(path):
    rows = read_tsv(path)
    species, runs = {}, set()
    for row in rows:
        if not all(row.get(k, "").strip() for k in ["species", "scientific_name", "taxid", "run"]):
            raise ValueError("samples require species, scientific_name, taxid and run")
        if row["run"] in runs:
            raise ValueError("duplicated run in samples: " + row["run"])
        runs.add(row["run"])
        int(row["taxid"])
        name = row["species"]
        if name in species and any(row[k] != species[name][k] for k in ["scientific_name", "taxid"]):
            raise ValueError("conflicting species metadata: " + name)
        species[name] = row
    if not species:
        raise ValueError("empty sample manifest")
    return rows, dict(sorted(species.items()))


def parse_tree(text, allowed, exact=False):
    from ete4 import Tree
    if not text.strip() or text.count(";") != 1:
        raise ValueError("expected one nonempty Newick tree")
    tree = Tree(text.strip(), parser=1)
    names = list(tree.leaf_names())
    if len(names) < 4 or len(names) != len(set(names)):
        raise ValueError("tree needs at least four unique species labels")
    if set(names) - set(allowed) or (exact and set(names) != set(allowed)):
        raise ValueError("tree species differ from sample manifest")
    for node in tree.traverse():
        if node.dist is not None and (not math.isfinite(node.dist) or node.dist < 0):
            raise ValueError("invalid branch length")
    return tree



EVENT_FIELDS = ["candidate_id", "species", "scientific_name", "run_ids", "rank", "role",
                "registered_taxid", "registered_taxon", "focal_taxid", "focal_taxon", "plot"]


def audit(tree, tree_qc, samples, taxonomy, outdir, settings=None):
    settings = validate_settings({} if settings is None else settings)
    inputs = {k: Path(v).resolve() for k, v in dict(tree=tree, tree_qc=tree_qc, samples=samples,
               taxonomy=taxonomy).items()}
    if Path(outdir).is_symlink():
        raise ValueError("audit output directory must not be a symlink")
    out = Path(outdir).resolve()
    if any(out == p or out in p.parents for p in inputs.values()):
        raise ValueError("audit output directory overlaps its inputs")
    if out.exists() and any(out.iterdir()):
        marker = out / "summary.json"
        if not marker.is_file() or json.loads(marker.read_text()).get("report_type") != "taxonomy_audit":
            raise ValueError("refusing to replace a directory not owned by taxonomy_audit")
    rows, species = load_samples(samples)
    names = list(species)
    source = parse_tree(inputs["tree"].read_text(), names, exact=True)
    tree_report = json.loads(inputs["tree_qc"].read_text())
    if tree_report.get("species") != len(names):
        raise ValueError("species-tree QC disagrees with its leaf count")
    outgroup = tree_report.get("outgroup")
    if not isinstance(outgroup, str) or outgroup not in species:
        raise ValueError("species-tree QC must record an outgroup present in the sample manifest")
    if len(source.children) != 2 or not any(n.is_leaf and n.name == outgroup for n in source.children):
        raise ValueError("species-tree root does not match the recorded single-species outgroup")
    records = {key: file_record(path) for key, path in inputs.items()}
    tax = Taxonomy(taxonomy)
    try:
        taxa = {name: tax.lineage(int(row["taxid"])) for name, row in species.items()}
    finally:
        tax.db.close()
    ranks = settings["ranks"]
    taxonomy_rows, lookup = [], {}
    for name in names:
        row = {"tip": name}
        for rank in ranks:
            value = taxa[name].get(rank)
            label = f"{value[1]} [{value[0]}]" if value else ""
            row[rank] = label
            if value:
                lookup[(rank, label)] = value
        taxonomy_rows.append(row)
    rscript = shutil.which("Rscript")
    if not rscript:
        raise RuntimeError("Rscript not found; activate workflow/envs/monophy.yaml and run its post-deploy script")
    runs = defaultdict(list)
    for row in rows:
        runs[row["species"]].append(row["run"])
    out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".taxonomy-audit-", dir=out.parent) as tmp:
        stage = Path(tmp) / "report"
        stage.mkdir()
        write_tsv(stage / "taxonomy.tsv", ["tip", *ranks], taxonomy_rows)
        request = dict(tree=str(inputs["tree"]),
                       outgroup=outgroup, outdir=str(stage), settings=settings, version=MONOPHY_VERSION,
                       scientific_names={name: row["scientific_name"] for name, row in species.items()})
        request_path = Path(tmp) / "request.json"
        write_json(request_path, request)
        script = Path(__file__).with_suffix(".R")
        environment = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
                       "MKL_NUM_THREADS": "1", "R_LIBS_USER": "NULL", "R_LIBS": "NULL",
                       "TZ": os.environ.get("TZ") or "UTC"}
        subprocess.run([rscript, "--vanilla", str(script), str(request_path)], check=True, env=environment)
        events = read_tsv(stage / "events.tsv")
        candidates, by_species = [], defaultdict(list)
        for event in events:
            name, rank, cid = event["species"], event["rank"], event["candidate_id"]
            expected = taxa[name][rank]
            focal = lookup[(rank, event["focal_taxon"])]
            candidates.append(dict(candidate_id=cid, species=name, scientific_name=species[name]["scientific_name"],
                run_ids=";".join(sorted(runs[name])), rank=rank, role=event["role"],
                registered_taxid=expected[0], registered_taxon=expected[1], focal_taxid=focal[0], focal_taxon=focal[1],
                plot=f"ranks/{rank}/tree.pdf"))
            by_species[name].append(cid)
        write_tsv(stage / "candidates.tsv", EVENT_FIELDS, candidates)
        rank_rows = read_tsv(stage / "rank_status.tsv")
        by_rank = defaultdict(list)
        for row in rank_rows:
            by_rank[row["species"]].append(row)
        sample_report = []
        for row in rows:
            name = row["species"]
            statuses = {r["status"] for r in by_rank[name]}
            status = ("review_flag" if by_species[name] else "non_monophyletic_group" if "non_monophyletic" in statuses
                      else "no_flag" if "monophyletic" in statuses else "not_assessable")
            sample_report.append({**row, "status": status, "candidate_ids": ";".join(by_species[name])})
        write_tsv(stage / "samples.tsv", list(sample_report[0]), sample_report)
        engine = json.loads((stage / "engine.json").read_text())
        summary = dict(report_type="taxonomy_audit", schema_version=2, method="MonoPhy", created_at=now(),
            report_only=True, settings=settings, species=len(names), runs=len(rows),
            outgroup=outgroup, candidate_species=len({c["species"] for c in candidates}), candidate_events=len(candidates),
            rank_status_counts=dict(Counter(r["status"] for r in rank_rows)), engine=engine, inputs=records,
            code={"python": file_record(__file__), "R": file_record(script),
                  "common": file_record(Path(__file__).with_name("common.py"))},
            limitations=[
                "MonoPhy reports taxonomy/tree discordance, not confirmed misidentification or corrected identities.",
                "An event is relative to its focal group; the same tip can be both an intruder and an outlier.",
                "Other species labels are not independently authenticated; grouped mislabels can implicate correct tips.",
                "Missing-rank tips are omitted separately at each rank; see rank_status and per-rank trees.",
                "One-tip groups cannot be assessed for their own monophyly, but can intrude other groups.",
                "The existing species-tree root is checked against the recorded outgroup.",
                "Monophyly depends on rooting and resolution; no branch-support cutoff or mislabel probability is used.",
                "All runs associated with a per-species CDS set are listed; a flag cannot localize an erroneous run.",
                "Plots are cladograms with registered-taxonomy tip colors; only unflagged monophyletic groups may be collapsed."])
        write_json(stage / "summary.json", summary)
        # Replace only a previously owned report, with rollback on publication failure.
        backup = Path(tmp) / "previous"
        if out.exists():
            out.rename(backup)
        try:
            stage.rename(out)
        except BaseException:
            if backup.exists():
                backup.rename(out)
            raise
    print(json.dumps({"method": "MonoPhy", "candidate_species": summary["candidate_species"], "output": str(out)}), flush=True)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["tree", "tree-qc", "samples", "taxonomy", "outdir"]:
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--settings", default="{}", help="JSON; see taxonomy_audit in config/config.yaml")
    args = vars(parser.parse_args())
    args["settings"] = json.loads(args["settings"])
    audit(**args)


if __name__ == "__main__":
    main()
