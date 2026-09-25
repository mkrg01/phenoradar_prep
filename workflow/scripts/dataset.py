#!/usr/bin/env python3
"""Build reusable species products from manual metadata through ODB mapping."""
import argparse
import copy
import csv
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import yaml

from common import file_record, now, read_tsv, write_json, write_tsv
from configuration import validate_analysis, validate_keys
from accession_exclusions import partition, read_exclusions
from phase_config import read_yaml
from dataset_software import resolve as resolve_software, validate as validate_software
from dataset_assets import (COUNTS, SAFE, digest, identities, import_existing, link_file, locked,
                            normalize_private_paths, record, register_busco, register_quant, register_reference, resolve, verify)

STAGES = ("assembly", "busco", "quant")
UNTIL = (*STAGES, "mapping")
MANAGED = {"mode_transcriptome_assembly", "kallisto_reference", "remove_amalgkit_fastq_after_completion", "delete_tmp_dir"}


def absolute(root, value):
    path = Path(value)
    return path.resolve() if path.is_absolute() else (Path(root) / path).resolve()


def inside(root, path):
    path = Path(path).resolve()
    if not path.is_relative_to(Path(root).resolve()):
        raise ValueError(f"dataset storage must be inside the project for container mounts: {path}")
    return path


def settings(root, config, analysis_config=None):
    from phase_config import validate_slurm
    root = Path(root).resolve()
    cfg = read_yaml(config)
    if analysis_config is not None:
        raise ValueError("build does not accept analysis overrides; use run_analysis.sh")
    unknown = set(cfg) - {"metadata", "store", "translation", "busco", "odb", "genegalleon", "slurm", "excluded_accessions"}
    if unknown: raise ValueError(f"unknown build settings: {sorted(unknown)}; migrate old dataset/config files first")
    excluded = cfg.get("excluded_accessions")
    if excluded is not None and (not isinstance(excluded, str) or not excluded.strip()):
        raise ValueError("excluded_accessions must be null or a TSV path")
    cfg["excluded_accessions"] = str(absolute(root, excluded)) if excluded is not None else None
    analysis = read_yaml(root / "workflow/pipeline_defaults.yaml")
    analysis["translation"] = cfg["translation"]
    if set(cfg["busco"]) != {"lineage"} or not isinstance(cfg["busco"]["lineage"], str) or not cfg["busco"]["lineage"].strip():
        raise ValueError("busco.lineage must be a nonempty string")
    if set(cfg["odb"]) - {"node", "cache_dir", "existing_results", "chunk_size"}:
        raise ValueError("unknown build.odb settings")
    analysis["odb"].update(cfg["odb"], incremental=True)
    analysis["phylogeny"]["lineage"] = cfg["busco"]["lineage"]
    analysis["phylogeny"]["trees"] = []
    for key in ("contrast_pairs", "dating", "taxonomy_check"):
        analysis["phylogeny"][key]["enabled"] = False
    analysis["selection"] = {"species_list": False, "busco_threshold": 0, "missing_taxonomy": "allow"}
    analysis["exclude_species"] = []
    validate_keys(analysis); validate_analysis(analysis)
    if type(analysis["translation"].get("table")) is not int or analysis["translation"]["table"] < 1:
        raise ValueError("translation.table must be a positive integer")
    if type(analysis["odb"]["node"]) is not int or analysis["odb"]["node"] < 1:
        raise ValueError("odb.node must be a positive integer")
    gg = cfg["genegalleon"]
    validate_software(gg)
    for key, value in gg.get("settings", {}).items():
        if not re.fullmatch(r"[a-z][a-z0-9_]*", key) or key.startswith("run_") or key in MANAGED:
            raise ValueError(f"managed/invalid GeneGalleon setting: {key}")
        if not isinstance(value, (str, int, float, bool)):
            raise ValueError(f"GeneGalleon setting must be scalar: {key}")
    validate_slurm(cfg["slurm"])
    cfg["store"] = str(inside(root, absolute(root, cfg["store"])))
    gg["cache_dir"] = str(inside(root, absolute(root, gg.get("cache_dir", "resources/software/genegalleon"))))
    for key in ("repository", "image"):
        if gg.get(key): gg[key] = str(absolute(root, gg[key]))
    cfg["conditions"] = stage_conditions(cfg)
    return cfg, analysis


def stage_conditions(cfg):
    gg = cfg["genegalleon"]
    software = {k: gg.get(k) for k in ("version", "revision", "image_uri", "image_sha256")}
    # Overrides are also bound by content, not merely their path.
    from dataset_software import source_records
    if gg.get("repository"):
        software["repository"] = [{"path": str(Path(p["path"]).relative_to(gg["repository"])), "sha256": p["sha256"]}
                                  for p in source_records(gg["repository"])]
    image = gg.get("image") or (str(Path(gg["repository"]) / "genegalleon.sif") if gg.get("repository") else None)
    if image and Path(image).is_file(): software["image"] = file_record(image)["sha256"]
    base = {"software": software, "settings": gg.get("settings", {}), "translation": cfg["translation"]}
    return {s: digest(dict(base, **({"lineage": cfg["busco"]["lineage"]} if s == "busco" else {}))) for s in STAGES}


def check_conditions(products, conditions):
    if not conditions: return
    for stage, key in zip(STAGES, ("reference", "busco", "quant")):
        product = products[key]
        previous = product.get("provenance", {}).get("condition") if product else None
        if previous and previous != conditions[stage]:
            raise ValueError(f"{stage} settings differ from registered product; choose a separate store for a deliberate rebuild")


def inspect_items(store, items, analysis, requested=None, conditions=None):
    result = []
    for item in items:
        try:
            products = resolve(store, item, analysis["phylogeny"]["lineage"], need_full=True)
            check_conditions(products, conditions)
            status = {stage: "reuse" if products[key] else "pending"
                      for stage, key in zip(STAGES, ("reference", "busco", "quant"))}
            result.append({"species": item["species"], "run": item["row"]["run"], **status,
                           "reason": "", "mapping": "pending_inputs" if not products["reference"] else "check_after_translation", "reference_id": products["reference"]["reference_id"] if products["reference"] else ""})
        except (ValueError, OSError, KeyError) as error:
            result.append({"species": item["species"], "run": item["row"]["run"],
                           **{stage: "conflict" for stage in STAGES}, "mapping": "conflict", "reason": str(error), "reference_id": ""})
    return result


def excluded_report(excluded):
    return [{"species": row["species"], "run": row["run"], **{s: "excluded" for s in STAGES},
             "reason": "excluded_accession: " + row["reason"], "mapping": "excluded", "reference_id": ""}
            for row in excluded]


def plan(root, config, metadata=None, analysis_config=None):
    cfg, analysis = settings(root, config, analysis_config)
    metadata = absolute(root, metadata or cfg["metadata"])
    source = record(metadata)
    policy_record = record(cfg["excluded_accessions"]) if cfg["excluded_accessions"] else None
    policy = read_exclusions(cfg["excluded_accessions"])
    fields, items = identities(metadata)
    items, excluded = partition(items, policy)
    normalize_private_paths(items, metadata)
    verify(source)
    if policy_record: verify(policy_record)
    selection = {"source_metadata": source, "exclusion_policy": policy_record, "excluded": excluded}
    report = inspect_items(cfg["store"], items, analysis, conditions=cfg["conditions"]) + excluded_report(excluded)
    return cfg, analysis, metadata, fields, items, selection, report


def implementation(root):
    # Bind the batch to code, not mutable branch names. Heavy data are recorded separately.
    paths = [*Path(root, "workflow").rglob("*.py"), *Path(root, "workflow").rglob("*.smk"),
             Path(root, "workflow/Snakefile"), Path(root, "run_pipeline.sh"),
             *Path(root, "workflow").rglob("*.yaml"), *Path(root, "workflow").rglob("*.sh")]
    if Path(root, "VERSION").exists(): paths.append(Path(root, "VERSION"))
    return [record(p) for p in sorted(paths)]


def prepare(root, name, config, metadata=None, analysis_config=None):
    root = Path(root).resolve()
    if not SAFE.fullmatch(name): raise ValueError("dataset name must be a simple directory name")
    target = root / "builds" / name
    if target.exists() or (root / "results" / ("build_" + name)).exists():
        raise ValueError("dataset/run name already exists; resume it or choose a new name")
    cfg, analysis, metadata, fields, items, selection, report = plan(root, config, metadata, analysis_config)
    if not items: raise ValueError("all metadata runs are excluded; no build was prepared")
    conflicts = [r for r in report if r["assembly"] == "conflict"]
    if conflicts: raise ValueError("resolve conflicts before preparing: " + json.dumps(conflicts))
    target.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".prepare-", dir=target.parent))
    try:
        frozen = staging / "metadata.tsv"
        # Normalize local FASTQ paths once; do not infer biological metadata.
        for item in items:
            row = item["row"]
            if row.get("private_file", "").lower() == "yes":
                if row.get("lib_layout") not in {"single", "paired"}: raise ValueError("private runs require lib_layout")
                required = ["read1_path"] + (["read2_path"] if row["lib_layout"] == "paired" else [])
                for key in required:
                    if not row.get(key): raise ValueError(f"private run requires {key}")
                    row[key] = str(absolute(metadata.parent, row[key]))
        write_tsv(frozen, fields, [i["row"] for i in items])
        gg_records, software_lock = [], None
        needs_upstream = any(r[s] == "pending" for r in report for s in STAGES)
        # Pin chosen references so a later import cannot change this dataset's reference.
        for item, state in zip(items, report):
            item["reference_id"] = state["reference_id"]
        raw_records = {}
        for item, state in zip(items, report):
            if item["row"].get("private_file", "").lower() == "yes" and any(state[s] == "pending" for s in ("assembly", "quant")):
                raw_records[item["species"]] = {k: record(item["row"][k]) for k in ("read1_path", "read2_path") if item["row"].get(k)}
        auxiliary = {}
        for filename, entry in (("source_metadata.tsv", selection["source_metadata"]),
                                ("excluded_accessions.tsv", selection["exclusion_policy"])):
            if entry:
                shutil.copy2(verify(entry), staging / filename)
                copied = record(staging / filename)
                if copied["sha256"] != entry["sha256"]: raise ValueError(f"input changed while freezing: {filename}")
                auxiliary[filename] = dict(copied, path=str(target / filename))
        write_tsv(staging / "excluded_runs.tsv", ["species", "run", "reason"], selection["excluded"])
        auxiliary["excluded_runs.tsv"] = dict(record(staging / "excluded_runs.tsv"), path=str(target / "excluded_runs.tsv"))
        if cfg["excluded_accessions"]: cfg["excluded_accessions"] = str(target / "excluded_accessions.tsv")
        analysis["run_name"] = "build_" + name
        analysis["input_root"] = str(target / "input")
        analysis["odb"]["cache_dir"] = str(inside(root, absolute(root, analysis["odb"]["cache_dir"])))
        if analysis["odb"]["existing_results"]:
            analysis["odb"]["existing_results"] = str(absolute(root, analysis["odb"]["existing_results"]))
        if needs_upstream:
            cfg["genegalleon"], gg_records, software_lock = resolve_software(cfg["genegalleon"])
        manifest = {"schema_version": 2, "kind": "build", "name": name, "created_at": now(), "root": str(root),
                    "config": cfg, "analysis": analysis, "fields": fields, "items": items,
                    "metadata": dict(record(frozen), path=str(target / "metadata.tsv")),
                    "auxiliary": auxiliary, "excluded": selection["excluded"], "raw_inputs": raw_records, "implementation": implementation(root),
                    "genegalleon": gg_records, "software_lock": software_lock}
        write_json(staging / "build.json", manifest)
        write_tsv(staging / "plan.tsv", list(report[0]), report)
        with (staging / "pipeline.yaml").open("w") as handle:
            yaml.safe_dump(analysis, handle, sort_keys=False)
        write_json(staging / "checksums.json", {p.name: file_record(p)["sha256"] for p in
                                              (staging / "build.json", staging / "pipeline.yaml")})
        os.rename(staging, target)
    finally:
        if staging.exists(): shutil.rmtree(staging)
    return target


def load(path, check_code=False):
    path = Path(path).resolve()
    hashes = json.loads((path / "checksums.json").read_text())
    for name, expected in hashes.items():
        if file_record(path / name)["sha256"] != expected: raise ValueError(f"frozen dataset changed: {name}")
    manifest = json.loads((path / "build.json").read_text())
    verify(manifest["metadata"])
    for entry in manifest["auxiliary"].values(): verify(entry)
    if check_code:
        for entry in manifest["implementation"] + manifest["genegalleon"]: verify(entry)
    return manifest


def item_products(manifest, item):
    bound = copy.deepcopy(item)
    if item.get("reference_id"): bound["row"]["reference_id"] = item["reference_id"]
    products = resolve(manifest["config"]["store"], bound, manifest["analysis"]["phylogeny"]["lineage"], need_full=True)
    check_conditions(products, manifest["config"].get("conditions"))
    return products


def status(path):
    manifest = load(path)
    items = copy.deepcopy(manifest["items"])
    for item in items:
        if item.get("reference_id"): item["row"]["reference_id"] = item["reference_id"]
    report = inspect_items(manifest["config"]["store"], items, manifest["analysis"], conditions=manifest["config"].get("conditions"))
    for item, row in zip(items, report):
        row["jobs"] = {}
        for stage in STAGES:
            receipt = Path(path) / "jobs/status" / f"{item['species']}.{stage}.json"
            if receipt.exists(): row["jobs"][stage] = json.loads(receipt.read_text())
        if row["assembly"] != "conflict":
            product = item_products(manifest, item)
            assessment = product.get("assessment") or product.get("busco")
            if assessment:
                c = assessment["counts"]
                row["busco_complete_fraction"] = (c[COUNTS[0]] + c[COUNTS[1]]) / c[COUNTS[-1]]
    completion = Path(path) / "completed.json"
    if completion.exists():
        from build_products import load_complete
        load_complete(completion)
        for row in report: row["mapping"] = "reuse"
    else:
        mapping_plan = Path(manifest["root"]) / "results" / manifest["analysis"]["run_name"] / "orthogroups/mapping/incremental_plan/plan.json"
        if mapping_plan.exists():
            planned = json.loads(mapping_plan.read_text())
            for row in report:
                if row["assembly"] != "conflict":
                    row["mapping"] = "planned_reuse" if row["species"] in planned["reused_species"] else "planned_mapping"
    return report + excluded_report(manifest.get("excluded", []))


def workspace(path):
    return Path(path) / "genegalleon"


def stage_workspace(path, manifest):
    work = workspace(path)
    for directory in ("input/amalgkit_metadata", "input/species_cds", "output", "downloads"):
        (work / directory).mkdir(parents=True, exist_ok=True)
    # Keep the complete ordered list fixed for every array and retry.
    for item in manifest["items"]:
        row = workspace_metadata(work, item, manifest["raw_inputs"].get(item["species"], {}), create=True)
        destination = work / "input/amalgkit_metadata" / (item["species"] + "_metadata.tsv")
        if destination.exists():
            if read_tsv(destination) != [row]: raise ValueError("staged GeneGalleon metadata changed")
        else:
            write_tsv(destination, manifest["fields"], [row])
    return work


def workspace_metadata(work, item, raw_inputs, create=False):
    row = dict(item["row"])
    for key, entry in raw_inputs.items():
        source = Path(entry["path"])
        relative = Path("input/reads") / item["species"] / (key + "".join(source.suffixes))
        staged = work / relative
        if create and not staged.exists():
            link_file(verify(entry), staged)
        else:
            verify(dict(entry, path=str(staged)))
        row[key] = "/workspace/" + str(relative)
    return row


def gg_environment(manifest, item, products, stage, work, task_id):
    # Do not inherit unrelated GeneGalleon overrides from the submission shell.
    env = {k: v for k, v in os.environ.items() if not k.startswith(("GG_TRANSCRIPTOME_", "GG_COMMON_"))}
    gg = manifest["config"]["genegalleon"]
    overrides = dict(gg.get("settings", {}))
    ref = products["reference"]
    overrides.update({
        "mode_transcriptome_assembly": "metadata", "run_amalgkit_metadata_or_integrate": 0,
        "run_amalgkit_getfastq": int(stage in {"assembly", "quant"}),
        "run_assembly": int(stage == "assembly" and not ref),
        "run_longestcds": int(stage == "assembly" and not ref),
        "run_longestcds_fx2tab": 0, "run_longestcds_mmseqs2taxonomy": 0,
        "run_longestcds_contamination_removal": 0, "run_busco_isoforms": 0,
        "run_busco_longest_cds": int(stage == "busco"),
        "run_busco_contamination_removed_longest_cds": 0,
        "run_assembly_stat": int(stage == "assembly"),
        "run_amalgkit_quant": int(stage == "quant"), "run_amalgkit_merge": int(stage == "quant"),
        "run_multispecies_summary": 0, "remove_amalgkit_fastq_after_completion": 0,
        "kallisto_reference": "species_cds", "delete_tmp_dir": 0,
    })
    for key, value in overrides.items():
        env["GG_TRANSCRIPTOME_" + key.upper()] = str(int(value)) if isinstance(value, bool) else str(value)
    # GeneGalleon uses shell sort for array indexing; match Python's ASCII order
    # both outside and inside the container, independent of the login locale.
    env.update({"LC_ALL": "C", "SINGULARITYENV_LC_ALL": "C", "APPTAINERENV_LC_ALL": "C",
                "gg_workspace_dir": str(work), "gg_container_image_path": gg["image"],
                "GG_COMMON_BUSCO_LINEAGE": manifest["analysis"]["phylogeny"]["lineage"],
                "GG_COMMON_GENETIC_CODE": str(manifest["analysis"]["translation"]["table"]),
                "GG_ARRAY_TASK_ID": str(task_id), "SLURM_ARRAY_TASK_ID": str(task_id)})
    return env


def worker(path, stage, task_id):
    path = Path(path).resolve()
    manifest = load(path, check_code=True)
    if stage not in STAGES or not 1 <= task_id <= len(manifest["items"]): raise ValueError("invalid stage/task index")
    item = manifest["items"][task_id - 1]
    species = item["species"]
    receipt_dir = path / "jobs/status"
    receipt_path = receipt_dir / f"{species}.{stage}.json"
    # Serialize all work on a species, including jobs accidentally submitted twice.
    with locked(Path(manifest["config"]["store"]) / species / ".worker.lock"):
        products = item_products(manifest, item)
        key = dict(zip(STAGES, ("reference", "busco", "quant")))[stage]
        if products[key]:
            write_json(receipt_path, {"state": "reused", "at": now()})
            return
        work = workspace(path)
        if not (work / "input/amalgkit_metadata" / f"{species}_metadata.tsv").is_file():
            raise ValueError("prepare job workspace with submit --dry-run or submit before running workers")
        expected_row = workspace_metadata(work, item, manifest["raw_inputs"].get(species, {}))
        if read_tsv(work / "input/amalgkit_metadata" / f"{species}_metadata.tsv") != [expected_row]:
            raise ValueError("staged GeneGalleon metadata changed")
        ref = products["reference"]
        if stage != "assembly" and not ref: raise ValueError(f"assembly prerequisite incomplete: {species}")
        out = work / "output/transcriptome_assembly"
        if ref:
            link_file(verify(ref["cds"]), out / "longest_cds" / f"{species}_longestCDS.fa.gz")
            link_file(verify(ref["cds"]), work / "input/species_cds" / f"{species}_longestCDS.fa.gz")
        # Incomplete native outputs must not satisfy GeneGalleon's existence checks.
        # Preserve them for diagnosis, and rerun only this unfinished stage.
        patterns = {"assembly": [("assembled_transcripts_with_isoforms", f"{species}_isoform.fa.gz"), ("corset_clusters", f"{species}_corset.clusters.tsv"), ("corset_counts", f"{species}_corset.counts.tsv"), ("assembly_stat", f"{species}_assembly_stat.tsv"), ("longest_cds", f"{species}_longestCDS.fa.gz"), ("longest_cds_transcript", f"{species}_longestCDS.transcript.fa.gz")],
                    "busco": [("busco_full_longest_cds", f"{species}_busco.full.tsv"), ("busco_short_longest_cds", f"{species}_busco.short.txt")],
                    "quant": [("amalgkit_quant", species), ("amalgkit_merge", species)]}
        previous = receipt_path.exists()
        if previous and json.loads(receipt_path.read_text())["state"] in {"running", "failed"}:
            quarantine = path / "jobs/incomplete" / species / stage / str(len(list((path / "jobs/incomplete" / species / stage).glob("*"))))
            for directory, pattern in patterns[stage]:
                for source in (out / directory).glob(pattern):
                    quarantine.mkdir(parents=True, exist_ok=True)
                    source.rename(quarantine / (directory + "-" + source.name))
        write_json(receipt_path, {"state": "running", "started_at": now(), "species": species, "run": item["row"]["run"],
                                      "stage": stage, "job_id": os.environ.get("SLURM_JOB_ID")})
        try:
            ordered = sorted(i["species"] + "_metadata.tsv" for i in manifest["items"])
            actual = sorted(p.name for p in (work / "input/amalgkit_metadata").iterdir()
                            if p.is_file() and not p.name.startswith("."))
            if actual != ordered: raise ValueError("staged GeneGalleon metadata file set changed")
            gg_index = ordered.index(species + "_metadata.tsv") + 1
            env = gg_environment(manifest, item, products, stage, work, gg_index)
            repository = Path(manifest["config"]["genegalleon"]["repository"])
            command = ["bash", str(repository / "workflow/gg_transcriptome_generation_entrypoint.sh")]
            subprocess.run(command, cwd=repository, env=env, check=True)
            provenance = {"source": "genegalleon", "dataset": str(path), "stage": stage,
                          "settings": {k: v for k, v in env.items() if k.startswith(("GG_TRANSCRIPTOME_", "GG_COMMON_"))},
                          "software": manifest["genegalleon"], "run": item["row"]["run"],
                          "raw_inputs": manifest["raw_inputs"].get(species, {}),
                          "condition": manifest["config"]["conditions"][stage]}
            store = manifest["config"]["store"]
            if stage == "assembly":
                ref = register_reference(store, item, out / "longest_cds" / f"{species}_longestCDS.fa.gz", provenance)
            elif stage == "busco":
                register_busco(store, ref, full=out / "busco_full_longest_cds" / f"{species}_busco.full.tsv",
                               short=out / "busco_short_longest_cds" / f"{species}_busco.short.txt",
                               lineage=manifest["analysis"]["phylogeny"]["lineage"], provenance=provenance)
            else:
                run = item["row"]["run"]
                abundance = out / "amalgkit_quant" / species / run / f"{run}_abundance.tsv"
                # Native counts, effective lengths, H5/JSON and merged tables remain in this persistent workspace.
                for suffix in ("eff_length", "est_counts", "tpm", "metadata"):
                    p = out / "amalgkit_merge" / species / f"{species}_{suffix}.tsv"
                    if not p.is_file() or not p.stat().st_size: raise ValueError(f"quant/merge did not finish: {p}")
                register_quant(store, ref, item, abundance, provenance)
            write_json(receipt_path, {"state": "complete", "finished_at": now(), "reference_id": ref["reference_id"],
                                      "job_id": os.environ.get("SLURM_JOB_ID")})
        except BaseException as error:
            write_json(receipt_path, {"state": "failed", "at": now(), "error": str(error), "species": species,
                                      "run": item["row"]["run"], "stage": stage, "job_id": os.environ.get("SLURM_JOB_ID")})
            raise


def materialize(path):
    path = Path(path).resolve()
    manifest = load(path, check_code=True)
    target = path / "input"
    with locked(path / ".materialize.lock"):
        if (path / "input_receipt.json").exists():
            receipt = json.loads((path / "input_receipt.json").read_text())
            for entry in receipt["files"]: verify(entry)
            return target
        rows, summaries, selected = [], [], []
        excluded = manifest.get("excluded", [])
        products = {}
        for item in manifest["items"]:
            p = item_products(manifest, item)
            missing = [k for k in ("reference", "busco", "quant") if not p[k]]
            if missing: raise ValueError(f"dataset incomplete: {item['species']}: {', '.join(missing)}")
            selected.append(item); products[item["species"]] = p
        if not selected: raise ValueError("build metadata has no species")
        staging = Path(tempfile.mkdtemp(prefix=".input-", dir=path))
        try:
            for item in selected:
                species, run = item["species"], item["row"]["run"]
                p = products[species]
                link_file(verify(p["reference"]["cds"]), staging / "cds" / f"{species}_longestCDS.fa.gz")
                link_file(verify(p["quant"]["abundance"]), staging / "quant" / species / run / f"{run}_abundance.tsv")
                if p["busco"].get("full"):
                    source = verify(p["busco"]["full"])
                    dest = staging / "busco/full" / f"{species}.busco.full.tsv"
                    if source.suffix == ".gz":
                        import gzip
                        dest.parent.mkdir(parents=True, exist_ok=True)
                        with gzip.open(source, "rb") as src, dest.open("wb") as dst: shutil.copyfileobj(src, dst)
                    else: link_file(source, dest)
                rows.append(item["row"])
                summaries.append({"Species": item["row"]["scientific_name"], **p["busco"]["counts"]})
            write_tsv(staging / "metadata.tsv", manifest["fields"], rows)
            write_tsv(staging / "busco/summary.tsv", ["Species", *COUNTS], summaries)
            if target.exists(): raise ValueError("unpublished input directory exists; inspect it before retrying")
            os.rename(staging, target)
            receipt = {"created_at": now(), "species": [i["species"] for i in selected], "excluded": excluded,
                       "files": [record(p) for p in sorted(target.rglob("*")) if p.is_file()]}
            write_json(path / "input_receipt.json", receipt)
        finally:
            if staging.exists(): shutil.rmtree(staging)
    return target


def run_mapping(path, execution=None):
    from build_products import complete
    from phase_config import write_profile
    path = Path(path).resolve()
    manifest = load(path, check_code=True)
    materialize(path)
    root = Path(manifest["root"])
    slurm = load_execution(execution) if execution else manifest["config"]["slurm"]
    profile = write_profile(path / "jobs" / (Path(execution).stem if execution else "local") / "profile", slurm)
    command = [str(root / "run_pipeline.sh"), "--slurm", "--profile", str(profile),
               "--jobs", str(slurm["jobs"]), "--configfile", str(path / "pipeline.yaml")]
    subprocess.run([*command, "--", "mapping"], cwd=root, check=True)
    complete(path)


def load_execution(path):
    wrapper = json.loads(Path(path).read_text())
    if digest(wrapper["slurm"]) != wrapper["sha256"]: raise ValueError("submission resources changed")
    return wrapper["slurm"]


def submit(path, until="mapping", species=None, dry_run=False, resources=None):
    path = Path(path).resolve()
    manifest = load(path, check_code=True)
    wanted = set(Path(species).read_text().splitlines()) if species else None
    if wanted is not None and (not wanted or wanted - {i["species"] for i in manifest["items"]}):
        raise ValueError("pilot species must be present in frozen metadata")
    report = status(path)
    if (path / "completed.json").exists():
        print("Build already complete; all recorded products verified.")
        return []
    if any(r["assembly"] == "conflict" for r in report): raise ValueError("resolve reported input conflicts before submitting")
    stop = STAGES.index(until) if until in STAGES else len(STAGES) - 1
    pending = {stage: [index for index, state in enumerate(report, 1)
                       if state[stage] == "pending" and (wanted is None or state["species"] in wanted)]
               for stage in STAGES[:stop + 1]}
    from phase_config import execution_settings
    slurm = execution_settings(manifest["config"]["slurm"], resources)
    jobs = path / "jobs"
    jobs.mkdir(exist_ok=True)
    (jobs / "logs").mkdir(exist_ok=True)
    with locked(jobs / ".submit.lock"):
        records = sorted(p for p in jobs.glob("submission_*.json") if not p.name.endswith(".resources.json"))
        if not dry_run:
            for record_path in records:
                old = json.loads(record_path.read_text())
                if any(r.get("state") in {"submitting", "unknown"} for r in old["jobs"]):
                    raise ValueError(f"unresolved submission in {record_path}; inspect Slurm and record its job ID before retrying")
                ids = [r["job_id"] for r in old["jobs"] if r.get("job_id")]
                if ids:
                    import pwd
                    queued = subprocess.check_output(["squeue", "--noheader", "--user", pwd.getpwuid(os.getuid()).pw_name,
                                                      "--format=%i"], text=True).splitlines()
                    active = [job.strip() for job in queued if job.strip().split("_", 1)[0] in ids]
                    if active: raise ValueError("dataset still has queued/running jobs; inspect or cancel them before resubmitting: " + ", ".join(active))
        if any(pending.values()): stage_workspace(path, manifest)
        batch_number = 1 + max([int(p.name.split("_")[1].split(".")[0]) for p in jobs.glob("submission_*.resources.json")] + [0])
        execution = jobs / f"submission_{batch_number:04d}.resources.json"
        write_json(execution, {"slurm": slurm, "sha256": digest(slurm)})
        batch = {"created_at": now(), "until": until, "pilot_species": sorted(wanted) if wanted else None, "resources": record(execution), "jobs": []}
        receipt = jobs / f"submission_{batch_number:04d}.json"
        previous = None
        commands = []
        scheduled = []
        array_size = slurm["array_size"]
        for stage, indices in pending.items():
            for offset in sorted({((i - 1) // array_size) * array_size for i in indices}):
                scheduled.append((stage, [i for i in indices if offset < i <= offset + array_size], offset))
        if until == "mapping" and wanted is None:
            scheduled.append(("mapping", None, 0))
        for stage, indices, offset in scheduled:
            job_resources = slurm["stages"]["controller" if stage == "mapping" else stage]
            label = stage + (f"_{offset + 1}_{offset + array_size}" if offset else "")
            script = jobs / f"{batch_number:04d}_{label}.sh"
            python = str(Path(sys.executable).resolve())
            cli = Path(manifest["root"]) / "workflow/scripts/dataset.py"
            worker_command = [python, str(cli), "mapping" if stage == "mapping" else "worker", "--dataset", str(path)]
            if stage == "mapping":
                worker_command += ["--execution", str(execution)]
                invocation = shlex.join(worker_command)
            else:
                worker_command += ["--stage", stage]
                invocation = shlex.join(worker_command) + f' --task-id "$((SLURM_ARRAY_TASK_ID + {offset}))"'
            script.write_text("#!/usr/bin/env bash\nset -euo pipefail\n" + "exec " + invocation + "\n")
            script.chmod(0o755)
            cmd = ["sbatch", "--parsable", "--nodes=1", "--ntasks=1", f"--cpus-per-task={job_resources['cpus']}",
                   f"--mem={job_resources['mem_mb']}M", f"--time={job_resources['time']}", "--chdir=" + manifest["root"],
                   f"--job-name={manifest['name']}_{label}", f"--output={jobs}/logs/{batch_number:04d}_{label}_%A_%a.out",
                   f"--error={jobs}/logs/{batch_number:04d}_{label}_%A_%a.err"]
            for key in ("partition", "account"):
                if slurm.get(key): cmd.append(f"--{key}={slurm[key]}")
            if indices: cmd.append("--array=" + ",".join(str(i - offset) for i in indices) + "%" + str(slurm["concurrency"]))
            if previous:
                cmd.extend(["--dependency=afterok:" + previous, "--kill-on-invalid-dep=yes"])
            cmd.append(str(script))
            commands.append(cmd)
            job = {"stage": stage, "indices": indices, "array_offset": offset, "command": cmd, "job_id": None, "state": "submitting"}
            batch["jobs"].append(job)
            if not dry_run:
                # Record submission intent before invoking sbatch. Preserve already submitted IDs on later failure.
                write_json(receipt, batch)
                try:
                    output = subprocess.check_output(cmd, text=True).strip()
                except (OSError, subprocess.CalledProcessError):
                    job["state"] = "rejected"
                    write_json(receipt, batch)
                    raise
                identifier = output.split(";")[0]
                if not identifier.isdigit():
                    job["state"] = "unknown"
                    write_json(receipt, batch)
                    raise ValueError(f"unrecognized sbatch result; inspect queue before retrying: {output}")
                job["job_id"] = identifier
                job["state"] = "submitted"
                write_json(receipt, batch)
                previous = identifier
            else:
                previous = f"JOB_ID_{label}"
        if dry_run:
            for command in commands: print(shlex.join(command))
        elif not commands:
            print("No pending work for this selection and endpoint.")
        return commands


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("plan", "prepare", "register", "fetch-software"):
        command = sub.add_parser(name)
        command.add_argument("--root", default=".")
        command.add_argument("--config", default="config/build.yaml")
        if name != "fetch-software": command.add_argument("--metadata")
        if name == "prepare": command.add_argument("--name", required=True)
        if name == "register":
            command.add_argument("--input-dir")
            command.add_argument("--products", help="Portable completed products to register in this project")
            command.add_argument("--odb-results", action="append", default=[], help="Existing ODB snapshot to add to the automatic cache")
            command.add_argument("--odb-only", action="store_true", help="Register only --odb-results, without scanning species artifacts")
    for name in ("status", "submit", "materialize", "worker", "mapping", "complete"):
        command = sub.add_parser(name)
        command.add_argument("--build", "--dataset", dest="dataset", required=True)
        if name == "submit": command.add_argument("--until", choices=UNTIL, default="mapping")
        if name == "mapping": command.add_argument("--execution")
        if name == "submit":
            command.add_argument("--resources")
            command.add_argument("--species-list")
            command.add_argument("--dry-run", action="store_true")
        if name == "worker":
            command.add_argument("--stage", choices=STAGES, required=True)
            command.add_argument("--task-id", type=int, required=True)
    args = parser.parse_args()
    if args.command in {"plan", "prepare", "register", "fetch-software"}:
        root = Path(args.root).resolve()
        config = absolute(root, args.config)
        if args.command == "plan":
            *_, report = plan(root, config, args.metadata)
            writer = csv.DictWriter(sys.stdout, fieldnames=list(report[0]), delimiter="\t")
            writer.writeheader(); writer.writerows(report)
            return int(any(r["assembly"] == "conflict" for r in report))
        if args.command == "prepare":
            print(prepare(root, args.name, config, args.metadata))
        elif args.command == "fetch-software":
            cfg, _ = settings(root, config)
            resolved, _, software_lock = resolve_software(cfg["genegalleon"])
            print(json.dumps({"repository": resolved["repository"], "image": resolved["image"],
                              "source": software_lock["source"].get("identity", {"kind": "local_source"}),
                              "container": software_lock["container"].get("identity", {"kind": "local_image"})}, indent=2))
        else:
            cfg, analysis = settings(root, config)
            if args.products and (args.input_dir or args.metadata or args.odb_only or args.odb_results):
                parser.error("--products cannot be combined with --input-dir, --metadata, --odb-only or --odb-results")
            if args.odb_only and (not args.odb_results or args.input_dir or args.metadata):
                parser.error("--odb-only requires --odb-results and cannot use --input-dir or --metadata")
            from incremental_odb import import_snapshot
            cache = inside(root, absolute(root, cfg["odb"]["cache_dir"]))
            exclusions = read_exclusions(cfg["excluded_accessions"])
            odb = []
            if args.products:
                from portable_build import completion_path, register_products
                source = inside(root, completion_path(absolute(root, args.products)))
                rows, snapshot = register_products(source, cfg["store"], cache, cfg["busco"]["lineage"],
                                                  cfg["translation"], cfg["odb"]["node"], exclusions)
                odb.append(str(snapshot))
            elif args.odb_only:
                rows = []
            else:
                rows = import_existing(cfg["store"], absolute(root, args.input_dir or "input"),
                                       absolute(root, args.metadata or cfg["metadata"]), cfg["busco"]["lineage"],
                                       excluded_runs=exclusions)
            for source in args.odb_results:
                odb.append(str(import_snapshot(absolute(root, source), cache, node=cfg["odb"]["node"])))
            print(json.dumps({"registered": sum(r["status"] == "registered" for r in rows),
                              "no_cds": [r["species"] for r in rows if r["status"] == "no_cds"],
                              "excluded": [r for r in rows if r["status"] == "excluded"],
                              "odb_snapshots": odb}, indent=2))
    elif args.command == "status":
        report = status(args.dataset)
        print(json.dumps(report, indent=2))
        return int(any(r["assembly"] == "conflict" for r in report))
    elif args.command == "submit": submit(args.dataset, args.until, args.species_list, args.dry_run, args.resources)
    elif args.command == "materialize": print(materialize(args.dataset))
    elif args.command == "worker": worker(args.dataset, args.stage, args.task_id)
    elif args.command == "complete":
        from build_products import complete
        print(complete(args.dataset))
    else: run_mapping(args.dataset, execution=args.execution)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        print(f"build: {error}", file=sys.stderr)
        sys.exit(1)
