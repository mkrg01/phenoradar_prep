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
from layout import run_layout
from sample_identity import select_samples
from dataset_software import resolve as resolve_software, validate as validate_software
from dataset_assets import (COUNTS, SAFE, digest, identities, link_file, locked,
                            normalize_private_paths, record, register_busco, register_quant, register_reference, resolve, verify)

STAGES = ("assembly", "busco", "quant")
UNTIL = (*STAGES, "mapping", "database")
FIXED_GENEGALLEON_SETTINGS = {"amalgkit_rrna_filter": "no", "amalgkit_contam_filter": "no"}
MANAGED = {"mode_transcriptome_assembly", "kallisto_reference", "remove_amalgkit_fastq_after_completion", "delete_tmp_dir"} | FIXED_GENEGALLEON_SETTINGS.keys()


def absolute(root, value):
    path = Path(value)
    return path.resolve() if path.is_absolute() else (Path(root) / path).resolve()


def inside(root, path):
    path = Path(path).resolve()
    if not path.is_relative_to(Path(root).resolve()):
        raise ValueError(f"dataset storage must be inside the project for container mounts: {path}")
    return path


def settings(root, config, analysis_config=None, name=None):
    from phase_config import validate_slurm
    root = Path(root).resolve()
    cfg = read_yaml(config)
    if analysis_config is not None:
        raise ValueError("build does not accept analysis overrides; use run_analysis.sh")
    unknown = set(cfg) - {"name", "metadata", "reuse_from", "translation", "busco", "odb", "genegalleon", "slurm", "excluded_accessions", "storage"}
    if unknown: raise ValueError(f"unknown build settings: {sorted(unknown)}; use config/build.yaml for build settings")
    excluded = cfg.get("excluded_accessions")
    if excluded is not None and (not isinstance(excluded, str) or not excluded.strip()):
        raise ValueError("excluded_accessions must be null or a TSV path")
    cfg["excluded_accessions"] = str(absolute(root, excluded)) if excluded is not None else None
    analysis = read_yaml(root / "workflow/pipeline_defaults.yaml")
    cfg.setdefault("storage", {"keep_intermediates": False})
    analysis["storage"] = cfg["storage"]
    analysis["translation"] = cfg["translation"]
    if set(cfg["busco"]) != {"lineage"} or not isinstance(cfg["busco"]["lineage"], str) or not cfg["busco"]["lineage"].strip():
        raise ValueError("busco.lineage must be a nonempty string")
    if set(cfg["odb"]) - {"ncbi_tax_id", "chunk_size"}:
        raise ValueError("unknown build.odb settings")
    analysis["odb"].update(cfg["odb"], incremental=True)
    analysis["phylogeny"]["lineage"] = cfg["busco"]["lineage"]
    analysis["phylogeny"]["trees"] = []
    for key in ("contrast_pairs", "dating", "taxonomy_check"):
        analysis["phylogeny"][key]["enabled"] = False
    analysis["selection"] = {"species_list": False, "busco_threshold": 0}
    analysis["exclude_species"] = []
    validate_keys(analysis); validate_analysis(analysis)
    cfg["storage"].setdefault("keep_intermediates", False)
    if type(analysis["translation"].get("table")) is not int or analysis["translation"]["table"] < 1:
        raise ValueError("translation.table must be a positive integer")
    gg = cfg["genegalleon"]
    validate_software(gg)
    for key, value in gg.get("settings", {}).items():
        if not re.fullmatch(r"[a-z][a-z0-9_]*", key) or key.startswith("run_") or key in MANAGED:
            raise ValueError(f"managed/invalid GeneGalleon setting: {key}")
        if not isinstance(value, (str, int, float, bool)):
            raise ValueError(f"GeneGalleon setting must be scalar: {key}")
    cfg["slurm"] = validate_slurm(cfg["slurm"])
    from database_reuse import normalize
    cfg["reuse_from"] = normalize(root, cfg.get("reuse_from"))
    build_name = name or cfg.get("name")
    cache = build_directory(root, build_name, config) / "work/cache" if build_name else root / "results/.plan/work/cache"
    cfg["store"] = str(cache / "products")  # Internal, frozen stage receipts; not a user setting.
    analysis["translation_cache"] = str(cache / "proteins")
    analysis["expression_cache"] = str(cache / "expression")
    analysis["odb"]["cache_dir"] = str(cache / "odb")
    gg["cache_dir"] = str(inside(root, absolute(root, gg.get("cache_dir", "resources/software/genegalleon"))))
    for key in ("repository", "image"):
        if gg.get(key): gg[key] = str(absolute(root, gg[key]))
    cfg["conditions"] = stage_conditions(cfg)
    return cfg, analysis


def stage_conditions(cfg):
    gg = cfg["genegalleon"]
    software = {k: gg.get(k) for k in ("version", "revision", "image_uri")}
    # Preserve condition hashes for products built with the former null setting.
    software["image_sha256"] = None
    # Overrides are also bound by content, not merely their path.
    from dataset_software import source_records
    if gg.get("repository"):
        software["repository"] = [{"path": str(Path(p["path"]).relative_to(gg["repository"])), "sha256": p["sha256"]}
                                  for p in source_records(gg["repository"])]
    image = gg.get("image") or (str(Path(gg["repository"]) / "genegalleon.sif") if gg.get("repository") else None)
    if image and Path(image).is_file(): software["image"] = file_record(image)["sha256"]
    effective_settings = {**gg.get("settings", {}), **FIXED_GENEGALLEON_SETTINGS}
    base = {"software": software, "settings": effective_settings, "translation": cfg["translation"]}
    return {s: digest(dict(base, **({"lineage": cfg["busco"]["lineage"]} if s == "busco" else {}))) for s in STAGES}


def check_conditions(products, conditions):
    if not conditions: return
    for stage, key in zip(STAGES, ("reference", "busco", "quant")):
        product = products[key]
        previous = product.get("provenance", {}).get("condition") if product else None
        if previous and previous != conditions[stage]:
            raise ValueError(f"{stage} settings differ from registered product; set reuse_from: null and choose a new build name for a deliberate rebuild")


def inspect_items(store, items, analysis, requested=None, conditions=None, reusable=None, errors=None):
    result = []
    for item in items:
        try:
            if item["species"] in (errors or {}):
                raise ValueError(errors[item["species"]])
            if item["species"] in (reusable or {}):
                product = reusable[item["species"]]['product']
                result.append({"species": item["species"], "run": item["row"]["run"],
                               **{stage: "reuse" for stage in STAGES}, "reason": "", "mapping": "reuse",
                               "reference_id": product['cds']['sha256']})
                continue
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


def plan(root, config, metadata=None, analysis_config=None, name=None):
    cfg, analysis = settings(root, config, analysis_config, name)
    metadata = absolute(root, metadata or cfg["metadata"])
    source = record(metadata)
    policy_record = record(cfg["excluded_accessions"]) if cfg["excluded_accessions"] else None
    policy = read_exclusions(cfg["excluded_accessions"])
    fields, items = identities(metadata)
    items, excluded = partition(items, policy)
    normalize_private_paths(items, metadata)
    verify(source)
    if policy_record: verify(policy_record)
    from database_reuse import select
    reusable, errors, sources = select(cfg['reuse_from'], items, cfg)
    selection = {"source_metadata": source, "exclusion_policy": policy_record, "excluded": excluded,
                 "reusable": reusable, "reuse_sources": sources}
    report = inspect_items(cfg["store"], items, analysis, conditions=cfg["conditions"],
                           reusable=reusable, errors=errors) + excluded_report(excluded)
    return cfg, analysis, metadata, fields, items, selection, report


def implementation(root):
    # Bind the batch to code, not mutable branch names. Heavy data are recorded separately.
    paths = [*Path(root, "workflow").rglob("*.py"), *Path(root, "workflow").rglob("*.smk"),
             Path(root, "workflow/Snakefile"), Path(root, "run_pipeline.sh"),
             *Path(root, "workflow").rglob("*.yaml"), *Path(root, "workflow").rglob("*.sh")]
    if Path(root, "VERSION").exists(): paths.append(Path(root, "VERSION"))
    return [record(p) for p in sorted(paths)]


def build_directory(root, name, config):
    root = Path(root).resolve()
    if name is None:
        name = read_yaml(config).get("name")
    if name is None:
        raise ValueError("build name is required: set name in build config or pass --name")
    if not isinstance(name, str) or not SAFE.fullmatch(name):
        raise ValueError("build name must be a simple directory name")
    return root / "results" / name


def prepare(root, name, config, metadata=None, analysis_config=None):
    root = Path(root).resolve()
    target = build_directory(root, name, config)
    name = target.name
    if target.exists():
        raise ValueError("dataset/run name already exists; resume it or choose a new name")
    cfg, analysis, metadata, fields, items, selection, report = plan(root, config, metadata, analysis_config, name)
    cfg["name"] = name
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
        # Pin chosen references so later builds cannot change this dataset's reference.
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
        analysis["run_name"] = name
        relative = target.relative_to(root)
        analysis["output_root"] = str(relative / "work/database")
        analysis["work_root"] = str(relative / "work/mapping")
        analysis["log_root"] = str(relative / "logs/database")
        analysis["input_root"] = str(target / "work/input")
        from database_reuse import stage
        stage(staging, target, selection['reusable'], cfg)
        reuse_receipt = {"sources": selection['reuse_sources'],
                         "samples": {s: entry['source'] for s, entry in selection['reusable'].items()}}
        write_json(staging / 'reuse.json', reuse_receipt)
        auxiliary['reuse.json'] = dict(record(staging / 'reuse.json'), path=str(target / 'reuse.json'))
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
        mapping_plan = Path(manifest["root"]) / run_layout(manifest["analysis"])[0] / "orthogroups/mapping/incremental_plan/plan.json"
        if mapping_plan.exists():
            planned = json.loads(mapping_plan.read_text())
            for row in report:
                if row["assembly"] != "conflict":
                    row["mapping"] = "planned_reuse" if row["species"] in planned["reused_species"] else "planned_mapping"
    return report + excluded_report(manifest.get("excluded", []))


def workspace(path, manifest, item=None):
    base = Path(path) / ("work/genegalleon" if manifest["analysis"].get("output_root") else "genegalleon")
    return base / item["species"] if item else base


def stage_workspace(path, manifest):
    # GeneGalleon owns species-named files. Isolate every run before invoking it.
    for item in manifest["items"]:
        work = workspace(path, manifest, item)
        for directory in ("input/amalgkit_metadata", "input/species_cds", "output", "downloads"):
            (work / directory).mkdir(parents=True, exist_ok=True)
        row = workspace_metadata(work, item, manifest["raw_inputs"].get(item["species"], {}), create=True)
        destination = work / "input/amalgkit_metadata" / (item["row"]["species_id"] + "_metadata.tsv")
        if destination.exists():
            if read_tsv(destination) != [row]: raise ValueError("staged GeneGalleon metadata changed")
        else:
            write_tsv(destination, manifest["fields"], [row])
    return workspace(path, manifest)


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
    overrides = {**gg.get("settings", {}), **FIXED_GENEGALLEON_SETTINGS}
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


def cleanup_genegalleon(path, manifest, item, stage, products, *, apply=True):
    """Retire only owned scratch; products and native BUSCO receipts stay valid."""
    if manifest['config'].get('storage', {}).get('keep_intermediates', False):
        return
    from cleanup_work import cleanup, owned_path
    path = Path(path)
    work = workspace(path, manifest, item)
    relative = work.relative_to(path) / 'output/transcriptome_assembly'
    targets = []
    # The upstream stages share tmp/. Never clear another failed stage's work
    # just because a previously completed stage was requested again.
    errors = []
    try:
        other_failed = False
        for other in STAGES:
            receipt = path / 'jobs/status' / f"{item['species']}.{other}.json"
            if other != stage and receipt.exists():
                other_failed |= json.loads(receipt.read_text()).get('state') in {'running', 'failed'}
        if not other_failed:
            targets.append(relative / 'tmp')
        # Assembly FASTQs are also used for quantification. Require both products,
        # irrespective of which stage is being retried. Keep caller-owned input/reads.
        if products.get('reference') and products.get('quant'):
            reads = relative / 'amalgkit_getfastq' / item['row']['species_id']
            # Preserve getfastq logs/QC and never descend into a linked input tree.
            directory = owned_path(path, reads)
            if not directory.is_symlink():
                for root, _, names in os.walk(directory, followlinks=False):
                    for name in names:
                        if name.endswith(('.fastq', '.fastq.gz', '.fq', '.fq.gz', '.sra')):
                            targets.append((Path(root) / name).relative_to(path))
    except (OSError, ValueError) as error:
        targets = []
        errors.append({'path': str(relative), 'error': str(error)})
    if not apply:
        from cleanup_work import scratch_usage
        sizes = [scratch_usage(owned_path(path, target)) for target in targets]
        return {'targets': [str(t) for t in targets], 'errors': errors,
                **{key: sum(size[key] for size in sizes)
                   for key in ('files', 'allocated_bytes', 'reclaimable_bytes')}}
    return cleanup(path, targets, path / 'jobs/cleanup' / f"{item['species']}.{stage}.json", errors=errors)


def worker(path, stage, task_id):
    path = Path(path).resolve()
    manifest = load(path, check_code=True)
    if stage not in STAGES or not 1 <= task_id <= len(manifest["items"]): raise ValueError("invalid stage/task index")
    item = manifest["items"][task_id - 1]
    species = item["species"]
    native = item["row"]["species_id"]
    from relabel_sample import relabel
    receipt_dir = path / "jobs/status"
    receipt_path = receipt_dir / f"{species}.{stage}.json"
    # Serialize all work on a species, including jobs accidentally submitted twice.
    with locked(Path(manifest["config"]["store"]) / species / ".worker.lock"):
        products = item_products(manifest, item)
        key = dict(zip(STAGES, ("reference", "busco", "quant")))[stage]
        if products[key]:
            write_json(receipt_path, {"state": "reused", "at": now()})
            cleanup_genegalleon(path, manifest, item, stage, products)
            return
        work = workspace(path, manifest, item)
        if not (work / "input/amalgkit_metadata" / f"{native}_metadata.tsv").is_file():
            raise ValueError("prepare job workspace with submit --dry-run or submit before running workers")
        expected_row = workspace_metadata(work, item, manifest["raw_inputs"].get(species, {}))
        if read_tsv(work / "input/amalgkit_metadata" / f"{native}_metadata.tsv") != [expected_row]:
            raise ValueError("staged GeneGalleon metadata changed")
        ref = products["reference"]
        if stage != "assembly" and not ref: raise ValueError(f"assembly prerequisite incomplete: {species}")
        out = work / "output/transcriptome_assembly"
        if ref:
            # Native tools retain biological names inside this isolated workspace.
            cds_native = out / "longest_cds" / f"{native}_longestCDS.fa.gz"
            relabel(verify(ref["cds"]), cds_native, species, native)
            link_file(cds_native, work / "input/species_cds" / f"{native}_longestCDS.fa.gz")
        # Incomplete native outputs must not satisfy GeneGalleon's existence checks.
        # Preserve them for diagnosis, and rerun only this unfinished stage.
        patterns = {"assembly": [("assembled_transcripts_with_isoforms", f"{native}_isoform.fa.gz"), ("corset_clusters", f"{native}_corset.clusters.tsv"), ("corset_counts", f"{native}_corset.counts.tsv"), ("assembly_stat", f"{native}_assembly_stat.tsv"), ("longest_cds", f"{native}_longestCDS.fa.gz"), ("longest_cds_transcript", f"{native}_longestCDS.transcript.fa.gz")],
                    "busco": [("busco_full_longest_cds", f"{native}_busco.full.tsv"), ("busco_short_longest_cds", f"{native}_busco.short.txt")],
                    "quant": [("amalgkit_quant", native), ("amalgkit_merge", native)]}
        previous = receipt_path.exists()
        if previous and json.loads(receipt_path.read_text())["state"] in {"running", "failed"}:
            quarantine = path / "jobs/incomplete" / species / stage / str(len(list((path / "jobs/incomplete" / species / stage).glob("*"))))
            for directory, pattern in patterns[stage]:
                for source in (out / directory).glob(pattern):
                    quarantine.mkdir(parents=True, exist_ok=True)
                    source.rename(quarantine / (directory + "-" + source.name))
            # A failed registration can leave a relabelled file from this attempt.
            published_name = {"assembly": f"{species}_longestCDS.fa.gz",
                              "busco": f"{species}.busco.full.tsv",
                              "quant": f"{item['row']['run']}_abundance.tsv"}[stage]
            published_file = work / "products" / published_name
            if published_file.exists():
                quarantine.mkdir(parents=True, exist_ok=True)
                published_file.rename(quarantine / ("products-" + published_name))
        write_json(receipt_path, {"state": "running", "started_at": now(), "species": species, "run": item["row"]["run"],
                                      "stage": stage, "job_id": os.environ.get("SLURM_JOB_ID")})
        try:
            actual = sorted(p.name for p in (work / "input/amalgkit_metadata").iterdir()
                            if p.is_file() and not p.name.startswith("."))
            if actual != [native + "_metadata.tsv"]:
                raise ValueError("staged GeneGalleon metadata file set changed")
            env = gg_environment(manifest, item, products, stage, work, 1)
            repository = Path(manifest["config"]["genegalleon"]["repository"])
            command = ["bash", str(repository / "workflow/gg_transcriptome_generation_entrypoint.sh")]
            subprocess.run(command, cwd=repository, env=env, check=True)
            provenance = {"source": "genegalleon", "dataset": str(path), "stage": stage,
                          "settings": {k: v for k, v in env.items() if k.startswith(("GG_TRANSCRIPTOME_", "GG_COMMON_"))},
                          "software": manifest["genegalleon"], "run": item["row"]["run"],
                          "raw_inputs": manifest["raw_inputs"].get(species, {}),
                          "condition": manifest["config"]["conditions"][stage]}
            store = manifest["config"]["store"]
            published = work / "products"
            if stage == "assembly":
                cds = published / f"{species}_longestCDS.fa.gz"
                provenance["relabel"] = relabel(out / "longest_cds" / f"{native}_longestCDS.fa.gz", cds, native, species)
                ref = register_reference(store, item, cds, provenance)
            elif stage == "busco":
                full = published / f"{species}.busco.full.tsv"
                provenance["relabel"] = relabel(out / "busco_full_longest_cds" / f"{native}_busco.full.tsv", full, native, species, 'busco')
                register_busco(store, ref, full=full,
                               short=out / "busco_short_longest_cds" / f"{native}_busco.short.txt",
                               lineage=manifest["analysis"]["phylogeny"]["lineage"], provenance=provenance)
            else:
                run = item["row"]["run"]
                abundance = published / f"{run}_abundance.tsv"
                for suffix in ("eff_length", "est_counts", "tpm", "metadata"):
                    p = out / "amalgkit_merge" / native / f"{native}_{suffix}.tsv"
                    if not p.is_file() or not p.stat().st_size: raise ValueError(f"quant/merge did not finish: {p}")
                provenance["relabel"] = relabel(out / "amalgkit_quant" / native / run / f"{run}_abundance.tsv",
                                               abundance, native, species, 'quant')
                register_quant(store, ref, item, abundance, provenance)
            write_json(receipt_path, {"state": "complete", "finished_at": now(), "reference_id": ref["reference_id"],
                                      "job_id": os.environ.get("SLURM_JOB_ID")})
        except BaseException as error:
            write_json(receipt_path, {"state": "failed", "at": now(), "error": str(error), "species": species,
                                      "run": item["row"]["run"], "stage": stage, "job_id": os.environ.get("SLURM_JOB_ID")})
            raise
        cleanup_genegalleon(path, manifest, item, stage, item_products(manifest, item))


def materialize(path):
    path = Path(path).resolve()
    manifest = load(path, check_code=True)
    target = Path(manifest["analysis"]["input_root"])
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
        target.parent.mkdir(parents=True, exist_ok=True)
        staging = Path(tempfile.mkdtemp(prefix=".input-", dir=target.parent))
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
                summaries.append({"Species": species, **p["busco"]["counts"]})
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
    from phase_config import validate_slurm, workflow_jobs, write_profile
    path = Path(path).resolve()
    manifest = load(path, check_code=True)
    materialize(path)
    root = Path(manifest["root"])
    slurm = load_execution(execution) if execution else manifest["config"]["slurm"]
    slurm = validate_slurm(slurm)
    profile = write_profile(path / "jobs" / (Path(execution).stem if execution else "local") / "profile", slurm)
    command = [str(root / "run_pipeline.sh"), "--slurm", "--profile", str(profile),
               "--jobs", str(workflow_jobs(slurm)), "--configfile", str(path / "pipeline.yaml")]
    subprocess.run([*command, "--", "database" if manifest["analysis"].get("output_root") else "mapping"], cwd=root, check=True)
    complete(path)


def load_execution(path):
    wrapper = json.loads(Path(path).read_text())
    if digest(wrapper["slurm"]) != wrapper["sha256"]: raise ValueError("submission resources changed")
    return wrapper["slurm"]


def submit(path, until="database", species=None, dry_run=False, resources=None):
    path = Path(path).resolve()
    manifest = load(path, check_code=True)
    wanted = set(Path(species).read_text().splitlines()) if species else None
    if wanted is not None:
        if not wanted: raise ValueError("pilot species must be present in frozen metadata")
        wanted = select_samples([dict(i["row"], species=i["species"]) for i in manifest["items"]], wanted)
    report = status(path)
    if (path / "completed.json").exists():
        print("Build already complete; all recorded products verified.")
        return []
    if any(r["assembly"] == "conflict" for r in report): raise ValueError("resolve reported input conflicts before submitting")
    stop = STAGES.index(until) if until in STAGES else len(STAGES) - 1
    pending = {stage: [index for index, state in enumerate(report, 1)
                       if state[stage] == "pending" and (wanted is None or state["species"] in wanted)]
               for stage in STAGES[:stop + 1]}
    from phase_config import array_concurrency, execution_settings, resolve_array_size, worker_budgets
    slurm = execution_settings(manifest["config"]["slurm"], resources)
    caps = {stage: array_concurrency(slurm, stage) for stage, indices in pending.items() if indices}
    if until in {"mapping", "database"} and wanted is None: worker_budgets(slurm)
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
        array_size = None
        if any(pending.values()):
            array_size = resolve_array_size()
            print(f"Using Slurm array limit: {array_size} sample index slots per batch", flush=True)
            stage_workspace(path, manifest)
        batch_number = 1 + max([int(p.name.split("_")[1].split(".")[0]) for p in jobs.glob("submission_*.resources.json")] + [0])
        execution = jobs / f"submission_{batch_number:04d}.resources.json"
        write_json(execution, {"slurm": slurm, "sha256": digest(slurm),
                               "resolved_array_size": array_size})
        batch = {"created_at": now(), "until": until, "pilot_species": sorted(wanted) if wanted else None, "resources": record(execution), "jobs": []}
        receipt = jobs / f"submission_{batch_number:04d}.json"
        previous = None
        commands = []
        scheduled = []
        for stage, indices in pending.items():
            for offset in sorted({((i - 1) // array_size) * array_size for i in indices}):
                scheduled.append((stage, [i for i in indices if offset < i <= offset + array_size], offset))
        if until in {"mapping", "database"} and wanted is None:
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
                   f"--mem={job_resources['mem_gb'] * 1000}M", f"--time={job_resources['time']}", "--chdir=" + manifest["root"],
                   f"--job-name={manifest['name']}_{label}", f"--output={jobs}/logs/{batch_number:04d}_{label}_%A_%a.out",
                   f"--error={jobs}/logs/{batch_number:04d}_{label}_%A_%a.err"]
            if slurm.get("partition"): cmd.append(f"--partition={slurm['partition']}")
            if indices:
                array = ",".join(str(i - offset) for i in indices)
                if caps[stage] is not None: array += "%" + str(caps[stage])
                cmd.append("--array=" + array)
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


def submit_named(root, name=None, config=None, metadata=None, until="database", species=None, dry_run=False, resources=None):
    """Freeze a new build, or submit an existing name with its saved conditions."""
    root = Path(root).resolve()
    config = absolute(root, config or "config/build.yaml")
    target = build_directory(root, name, config)
    with locked(target.parent / f".{target.name}.prepare.lock"):
        if target.exists():
            if not (target / "build.json").is_file():
                raise ValueError(f"existing path is not a prepared build: {target}; choose a new name")
            print(f"Using saved build: {target}. Source edits require a new name.", flush=True)
        else:
            target = prepare(root, target.name, config, metadata)
            print(f"Prepared build: {target}", flush=True)
    return submit(target, until, absolute(root, species) if species else None, dry_run,
                  absolute(root, resources) if resources else None)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("plan", "prepare", "fetch-software"):
        command = sub.add_parser(name)
        command.add_argument("--root", default=".")
        command.add_argument("--config", default="config/build.yaml")
        if name != "fetch-software": command.add_argument("--metadata")
        if name == "prepare": command.add_argument("--name", help="Override name from the build config")
    for name in ("status", "cleanup", "submit", "materialize", "worker", "mapping", "complete"):
        command = sub.add_parser(name)
        command.add_argument("--build", "--dataset", dest="dataset", required=name != "submit",
                             help="Path to an already prepared build")
        if name == "submit": command.add_argument("--until", choices=UNTIL, default="database")
        if name == "mapping": command.add_argument("--execution")
        if name == "status": command.add_argument("--storage", action="store_true", help="Include intermediate cleanup status")
        if name == "cleanup": command.add_argument("--apply", action="store_true", help="Delete verified successful scratch; default is preview")
        if name == "submit":
            command.add_argument("--root", default=".")
            command.add_argument("--config", help="Settings for a new build; default: config/build.yaml")
            command.add_argument("--name", help="Build name; defaults to name in the build config")
            command.add_argument("--metadata", help="Metadata for a new build")
            command.add_argument("--resources")
            command.add_argument("--species-list")
            command.add_argument("--dry-run", action="store_true")
        if name == "worker":
            command.add_argument("--stage", choices=STAGES, required=True)
            command.add_argument("--task-id", type=int, required=True)
    args = parser.parse_args()
    if args.command in {"plan", "prepare", "fetch-software"}:
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
    elif args.command == "status":
        report = status(args.dataset)
        if args.storage:
            from storage_management import run_storage
            print(json.dumps({'samples': report, 'storage': run_storage(args.dataset, 'build')}, indent=2))
        else:
            print(json.dumps(report, indent=2))
        return int(any(r["assembly"] == "conflict" for r in report))
    elif args.command == "cleanup":
        from storage_management import run_storage
        report = run_storage(args.dataset, 'build', inspect=True, apply=args.apply)
        print(json.dumps(report, indent=2))
        return int(report['state'] == 'pending')
    elif args.command == "submit":
        if args.dataset:
            if any(value is not None for value in (args.name, args.config, args.metadata)):
                parser.error("--build cannot be combined with --name, --config, or --metadata; use --resources for retry resources")
            submit(absolute(args.root, args.dataset), args.until,
                   absolute(args.root, args.species_list) if args.species_list else None, args.dry_run,
                   absolute(args.root, args.resources) if args.resources else None)
        else:
            submit_named(args.root, args.name, args.config, args.metadata, args.until,
                         args.species_list, args.dry_run, args.resources)
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
