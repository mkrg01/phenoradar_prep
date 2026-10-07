"""Exercise incremental datasets with synthetic RNA-seq artifacts and scheduler doubles."""
import copy
import gzip
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from common import file_record, read_tsv, write_json, write_tsv
from dataset import (gg_environment, load, materialize, plan, prepare, status, submit, worker)
from dataset_assets import (COUNTS, identities, register_busco, register_quant,
                            register_reference, resolve)

from phase_config import resolve_array_size

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def dataset_project(tmp_path, tiny_inputs, monkeypatch):
    # General workflow tests must not depend on a real Slurm controller.
    # Query integration tests restore the real helper with a fake scontrol.
    monkeypatch.setattr('phase_config.resolve_array_size', lambda: 1000)
    root = tmp_path / "project"
    root.mkdir()
    shutil.copytree(ROOT / "config", root / "config")
    # Tests start without any external database.
    build_config = root / "config/build.yaml"
    cfg = yaml.safe_load(build_config.read_text())
    cfg["reuse_from"] = None
    cfg["metadata"] = "input/metadata.tsv"
    build_config.write_text(yaml.safe_dump(cfg))
    shutil.copytree(ROOT / "workflow", root / "workflow", ignore=shutil.ignore_patterns("__pycache__"))
    shutil.copy2(ROOT / "run_pipeline.sh", root / "run_pipeline.sh")
    shutil.copytree(Path(tiny_inputs["metadata"]).parent, root / "input")
    metadata = root / "input/metadata.tsv"
    rows = [r for r in read_tsv(metadata) if r["run"] != "A2"]
    write_tsv(metadata, list(rows[0]), rows)
    full = root / "input/busco/full"
    full.mkdir(parents=True, exist_ok=True)
    for row in read_tsv(root / "input/busco/summary.tsv"):
        species = row["Species"].replace(" ", "_")
        statuses = []
        for key, state in zip(COUNTS[:-1], ("Complete", "Duplicated", "Fragmented", "Missing")):
            statuses.extend([state] * int(row[key]))
        (full / f"{species}.busco.full.tsv").write_text("# The lineage dataset is: embryophyta_odb12\n" + "".join(
            f"M{i}\t{state}\t{species}_g1\t100\t3\n" for i,state in enumerate(statuses)))
    # New builds consume sample-prefixed products; low-level legacy fixtures remain separate.
    from relabel_sample import relabel
    from sample_identity import annotate
    original_rows = read_tsv(metadata)
    summaries = {r["Species"]: r for r in read_tsv(root / "input/busco/summary.tsv")}
    upgraded = []
    for raw in original_rows:
        row = annotate(raw)
        old, new, run = row["species_id"], row["analysis_sample_id"], row["run"]
        relabel(root / "input/cds" / f"{old}_longestCDS.fa.gz", root / "input/cds" / f"{new}_longestCDS.fa.gz", old, new)
        relabel(full / f"{old}.busco.full.tsv", full / f"{new}.busco.full.tsv", old, new, "busco")
        relabel(root / "input/quant" / old / run / f"{run}_abundance.tsv", root / "input/quant" / new / run / f"{run}_abundance.tsv", old, new, "quant")
        upgraded.append(dict(summaries[raw["scientific_name"]], Species=new))
    write_tsv(root / "input/busco/summary.tsv", list(upgraded[0]), upgraded)
    return root


@pytest.mark.parametrize("ncbi_tax_id", [3193, 33090])
def test_build_preserves_odb_ncbi_tax_id_in_frozen_configs(dataset_project, ncbi_tax_id):
    root = dataset_project
    fake_genegalleon(root)
    config = root / "config/build.yaml"
    cfg = yaml.safe_load(config.read_text())
    cfg["odb"]["ncbi_tax_id"] = ncbi_tax_id
    config.write_text(yaml.safe_dump(cfg))
    build = prepare(root, "taxid_test", config)
    manifest = load(build)
    pipeline = yaml.safe_load((build / "pipeline.yaml").read_text())
    for odb in (manifest["config"]["odb"], manifest["analysis"]["odb"], pipeline["odb"]):
        assert odb["ncbi_tax_id"] == ncbi_tax_id
        assert "node" not in odb


@pytest.mark.parametrize("policy", ["error", "drop", "split"])
def test_tpm_policy_is_not_a_build_setting(dataset_project, policy):
    config = dataset_project / "config/build.yaml"
    cfg = yaml.safe_load(config.read_text())
    cfg["tpm"] = {"multimap": policy}
    config.write_text(yaml.safe_dump(cfg))
    with pytest.raises(ValueError, match="unknown build settings:.*tpm"):
        plan(dataset_project, config)


def native_events(build):
    return [json.loads(line) for p in sorted((build / "work/genegalleon").glob("*/events.jsonl"))
            for line in p.read_text().splitlines()]


def imported(root):
    store = root / "resources/dataset_assets"
    # Seed completed worker products through the normal stage registry API.
    _, items = identities(root / "input/metadata.tsv")
    summaries = {r["Species"]: r for r in read_tsv(root / "input/busco/summary.tsv")}
    for item in items:
        name, run = item["species"], item["row"]["run"]
        cds = root / "input/cds" / f"{name}_longestCDS.fa.gz"
        if not cds.exists():
            continue
        ref = register_reference(store, item, cds)
        full = root / "input/busco/full" / f"{name}.busco.full.tsv"
        register_busco(store, ref, summaries[name], full=full if full.exists() else None)
        abundance = root / "input/quant" / name / run / f"{run}_abundance.tsv"
        if abundance.exists():
            register_quant(store, ref, item, abundance)
    from database_fixtures import database_from_stages
    database = database_from_stages(root, store, items)
    config = root / 'config/build.yaml'
    cfg = yaml.safe_load(config.read_text()); cfg['reuse_from'] = str(database)
    config.write_text(yaml.safe_dump(cfg))
    return database


@pytest.mark.parametrize("value", [None, "resources/odb_existing/tlight"])
def test_build_rejects_external_odb_import_setting(dataset_project, value):
    root = dataset_project
    config = root / "config/build.yaml"
    cfg = yaml.safe_load(config.read_text())
    cfg["odb"]["existing_results"] = value
    config.write_text(yaml.safe_dump(cfg))
    with pytest.raises(ValueError, match="unknown build.odb settings"):
        plan(root, config)


@pytest.mark.parametrize("override", [None, "pilot_20260925"])
def test_prepare_cli_uses_config_name_with_optional_override(dataset_project, override):
    root = dataset_project
    imported(root)
    cfg = yaml.safe_load((root / "config/build.yaml").read_text())
    cfg["name"] = "angiosperm_leaf_20260925"
    cfg["name_mode"] = "fixed"
    config = root / "config/named.yaml"
    config.write_text(yaml.safe_dump(cfg))
    command = [sys.executable, str(root / "workflow/scripts/dataset.py"),
               "prepare", "--config", "config/named.yaml"]
    if override is not None:
        command.extend(["--name", override])
    result = subprocess.run(command, cwd=root, text=True, capture_output=True)
    assert result.returncode == 0, result.stderr
    name = override or cfg["name"]
    path = root / "results" / name
    assert result.stdout.strip() == str(path)
    frozen = load(path, check_code=True)
    assert frozen["name"] == frozen["config"]["name"] == name
    assert frozen["analysis"]["run_name"] == name
    assert frozen["analysis"]["output_root"] == f"results/{name}/work/database"
    assert {row["assembly"] for row in status(path)} == {"reuse"}
    # Renaming in the source config cannot rename an existing preparation.
    cfg["name"] = "angiosperm_leaf_20260926"
    config.write_text(yaml.safe_dump(cfg))
    assert load(path, check_code=True)["name"] == name
    with pytest.raises(ValueError, match="name already exists"):
        prepare(root, name, config)


def test_prepare_legacy_config_requires_explicit_name(dataset_project):
    root = dataset_project
    imported(root)
    config = root / "config/build.yaml"
    cfg = yaml.safe_load(config.read_text())
    cfg.pop("name", None)
    cfg.pop("name_mode", None)
    config.write_text(yaml.safe_dump(cfg))
    with pytest.raises(ValueError, match="set name in build config or pass --name"):
        prepare(root, None, config)
    path = prepare(root, "explicit_20260925", config)
    assert load(path)["name"] == "explicit_20260925"


@pytest.mark.parametrize("name", ["", "../escape", "leaf/20260925", "two words", 20260925, True])
def test_prepare_rejects_invalid_config_name(dataset_project, name):
    root = dataset_project
    config = root / "config/build.yaml"
    cfg = yaml.safe_load(config.read_text())
    cfg["name"] = name
    config.write_text(yaml.safe_dump(cfg))
    with pytest.raises(ValueError, match="build name must be a simple directory name"):
        prepare(root, None, config)
    assert not (root / "results").exists()


def test_manual_metadata_allows_multiple_runs_with_unique_sample_identities(tmp_path):
    path = tmp_path / "metadata.tsv"
    fields = ["scientific_name", "run", "taxid"]
    write_tsv(path, fields, [dict(zip(fields, ["Alpha plant", f"R{i}", "42"])) for i in range(2)])
    assert len(identities(path)[1]) == 2
    write_tsv(path, fields, [dict(zip(fields, ["Alpha plant", "R1", "42"]))] * 2)
    with pytest.raises(ValueError, match="unique runs"):
        identities(path)


def test_database_and_changed_metadata_only_reuses_products(dataset_project):
    root = dataset_project
    imported(root)
    cfg = root / "config/build.yaml"
    report = plan(root, cfg)[-1]
    assert [(r["assembly"], r["busco"], r["quant"]) for r in report] == [
        ("reuse", "reuse", "reuse"), ("reuse", "reuse", "reuse"), ("reuse", "reuse", "reuse")]
    path = prepare(root, "base", cfg)
    assert submit(path, until="quant", dry_run=True) == []
    assert not (path / "work/genegalleon").exists()
    assert load(path)["analysis"]["odb"]["incremental"] is True
    assert load(path)["software_lock"] is None  # Reuse-only preparation does not acquire software.


def test_status_resolves_each_run_once_and_submit_loads_once(dataset_project, monkeypatch):
    import dataset
    root = dataset_project
    imported(root)
    path = prepare(root, "status_checks", root / "config/build.yaml")
    expected = status(path)
    calls = []
    resolve = dataset.resolve
    def counted(store, item, *args, **kwargs):
        calls.append(item["species"])
        return resolve(store, item, *args, **kwargs)
    monkeypatch.setattr(dataset, "resolve", counted)
    assert status(path) == expected
    assert len(calls) == len(load(path)["items"])
    loads = []
    load_build = dataset.load
    def counted_load(*args, **kwargs):
        loads.append(args[0])
        return load_build(*args, **kwargs)
    monkeypatch.setattr(dataset, "load", counted_load)
    assert submit(path, until="quant", dry_run=True) == []
    assert loads == [path.resolve()]


def test_removal_readdition_and_frozen_membership(dataset_project):
    root = dataset_project
    store = imported(root)
    cfg = root / "config/build.yaml"
    first = prepare(root, "base", cfg)
    first_input = materialize(first)
    assert {r["scientific_name"] for r in read_tsv(first_input / "metadata.tsv")} == {"Alpha plant", "Beta sp-X", "Gamma plant"}
    metadata = root / "input/metadata.tsv"
    original = read_tsv(metadata)
    write_tsv(metadata, list(original[0]), [original[0]])
    # Existing batches are unaffected by later edits to the source metadata.
    assert len(load(first)["items"]) == 3
    second = prepare(root, "removed", cfg)
    second_input = materialize(second)
    assert [r["scientific_name"] for r in read_tsv(second_input / "metadata.tsv")] == ["Alpha plant"]
    assert not (second_input / "cds/Beta_sp-X_B1_longestCDS.fa.gz").exists()
    assert [p.name for p in (second_input / "quant").iterdir()] == ["Alpha_plant_A1"]
    assert (store / "cds/Beta_sp-X_B1_longestCDS.fa.gz").exists()
    write_tsv(metadata, list(original[0]), original)
    third = prepare(root, "restored", cfg)
    assert submit(third, until="quant", dry_run=True) == []
    assert len(read_tsv(materialize(third) / "metadata.tsv")) == 3
    assert len(read_tsv(first_input / "metadata.tsv")) == 3


def test_changed_run_requires_independent_assembly_and_missing_is_not_silently_dropped(dataset_project):
    root = dataset_project
    imported(root)
    metadata = root / "input/metadata.tsv"
    rows = read_tsv(metadata)
    rows[0]["run"] = "Anew"
    write_tsv(metadata, list(rows[0]), rows)
    report = plan(root, root / "config/build.yaml")[-1]
    assert [report[0][s] for s in ("assembly", "busco", "quant")] == ["pending", "pending", "pending"]
    fake_genegalleon(root)
    path = prepare(root, "newrun", root / "config/build.yaml")
    with pytest.raises(ValueError, match="dataset incomplete: Alpha_plant_Anew: reference, busco, quant"):
        materialize(path)
    assert not (path / "work/input").exists()


def test_modified_registered_cds_is_a_conflict(dataset_project):
    root = dataset_project
    imported(root)
    cds = root / "input/cds/Alpha_plant_A1_longestCDS.fa.gz"
    with gzip.open(cds, "wt") as handle: handle.write(">Alpha_plant_A1_g1\nATGCCC\n")
    report = plan(root, root / "config/build.yaml")[-1]
    assert report[0]["assembly"] == "conflict"
    assert "registered file changed" in report[0]["reason"]


def test_frozen_metadata_and_configuration_are_verified(dataset_project):
    root = dataset_project
    imported(root)
    path = prepare(root, "immutable", root / "config/build.yaml")
    (path / "metadata.tsv").write_text("modified")
    with pytest.raises(ValueError, match="registered file changed"):
        status(path)


def fake_genegalleon(root):
    repo = root / "fake_gg"
    (repo / "workflow").mkdir(parents=True)
    (repo / "genegalleon.sif").write_text("test container identifier")
    implementation = r'''

import csv, gzip, hashlib, json, os, sys, time
from pathlib import Path
work = Path(os.environ['gg_workspace_dir'])
files = sorted((work / 'input/amalgkit_metadata').glob('*.tsv'))
metadata = files[int(os.environ['GG_ARRAY_TASK_ID']) - 1]
row = next(csv.DictReader(metadata.open(), delimiter='\t'))
species = row['scientific_name'].replace(' ', '_')
run = row['run']
prefix = 'GG_TRANSCRIPTOME_'
out = work / 'output/transcriptome_assembly'
out.mkdir(parents=True, exist_ok=True)
with (work / 'events.jsonl').open('a') as handle:
    handle.write(json.dumps({'species':species, 'env':{k:v for k,v in os.environ.items() if k.startswith(prefix)},
                            'budgets':{k:os.environ[k] for k in ('GG_TASK_CPUS', 'GG_MEM_TOTAL_GB', 'GG_MEM_TOOL_GB')}})+'\n')
assert os.environ['LC_ALL'] == os.environ['SINGULARITYENV_LC_ALL'] == os.environ['APPTAINERENV_LC_ALL'] == 'C'
assert os.environ[prefix+'RUN_MULTISPECIES_SUMMARY'] == '0'
assert os.environ['GG_OBSERVABILITY'] == '1'
assert os.environ['GG_COMMON_TMP_ROOT'] == 'workspace'
scratch = out / 'tmp/1_native/large'
scratch.parent.mkdir(parents=True, exist_ok=True)
scratch.write_bytes(b'scratch')
attempt = work / 'output/observations' / str(time.time_ns())
attempt.mkdir(parents=True)
record = {'schema':'genegalleon-observation-v1', 'attempt_id':attempt.name,
          'workflow':'gg_transcriptome_generation', 'started_at_ns':time.time_ns(), 'execution_state':'started'}
def save_run(): (attempt / 'run.json').write_text(json.dumps(record))
save_run()
def finish_run(code):
    record.update(execution_state='exited', exit_code=code)
    save_run()
    sys.exit(code)
def emit(step, operation, code, proof):
    data={'schema':'genegalleon-contract-observation-v1','attempt_id':attempt.name,
          'observed_at_ns':time.time_ns(),'step':step,'operation':operation,'exit_code':code,
          'manifest_sha256':'proof' if proof.exists() else None, 'proof':str(proof)}
    (attempt / ('contract-'+step+'.json')).write_text(json.dumps(data))
def begin(step):
    proof=out / 'artifact_provenance' / (step+'.json')
    if proof.exists():
        data=json.loads(proof.read_text())
        if not all(Path(p).is_file() and hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in data.items()):
            emit(step, 'needs-run', 3, proof)
            finish_run(17)
        emit(step, 'needs-run', 1, proof)
        return False
    emit(step, 'needs-run', 0, proof)
    with (work/'steps.jsonl').open('a') as log: log.write(json.dumps(step)+'\n')
    return True
def finish(step, paths):
    proof=out / 'artifact_provenance' / (step+'.json')
    proof.parent.mkdir(exist_ok=True)
    proof.write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [metadata,*paths]}))
    emit(step, 'record', 0, proof)
if os.environ[prefix+'RUN_AMALGKIT_GETFASTQ']=='1' and os.environ.get('FAKE_GG_FAIL_DOWNLOAD'):
    finish_run(29)
if os.environ[prefix+'RUN_ASSEMBLY']=='1' and begin('transcriptome_longest_cds'):
    if os.environ.get('FAKE_GG_FAIL_ASSEMBLY'): finish_run(23)
    cds=out/'longest_cds'/(species+'_longestCDS.fa.gz')
    cds.parent.mkdir(parents=True,exist_ok=True)
    with gzip.open(cds,'wt') as handle:
        handle.write('>'+species+'_g1\nATGAAATAA\n>'+species+'_g2\nATGCCCTAA\n')
    finish('transcriptome_longest_cds',[cds])
if os.environ[prefix+'RUN_BUSCO_LONGEST_CDS']=='1' and begin('transcriptome_busco_longest_cds'):
    if os.environ.get('FAKE_GG_FAIL_BUSCO'): finish_run(31)
    outputs=[]
    for directory,suffix,content in [('busco_full_longest_cds','full.tsv','# The lineage dataset is: embryophyta_odb12 (test)\nB1\tComplete\t'+species+'_g1:0-9\t100\t3\nB2\tComplete\t'+species+'_g2:0-9\t100\t3\n'),('busco_short_longest_cds','short.txt','C:100%[S:100%,D:0%],F:0%,M:0%,n:2\n')]:
        path=out/directory/(species+'_busco.'+suffix)
        if os.environ.get('FAKE_GG_LOW_BUSCO'):
            content=('# The lineage dataset is: embryophyta_odb12 (test)\nB1\tMissing\nB2\tMissing\n' if suffix=='full.tsv' else 'C:0%[S:0%,D:0%],F:0%,M:0%,n:2\n')
        elif os.environ.get('FAKE_GG_HALF_BUSCO'):
            content=('# The lineage dataset is: embryophyta_odb12 (test)\nB1\tDuplicated\t'+species+'_g1:0-9\t100\t3\nB1\tDuplicated\t'+species+'_g2:0-9\t100\t3\nB2\tMissing\n' if suffix=='full.tsv' else 'C:50%[S:0%,D:50%],F:0%,M:0%,n:2\n')
        path.parent.mkdir(parents=True,exist_ok=True); path.write_text(content); outputs.append(path)
    finish('transcriptome_busco_longest_cds',outputs)
if os.environ[prefix+'RUN_AMALGKIT_QUANT']=='1' and begin('transcriptome_quant'):
    if os.environ.get('FAKE_GG_FAIL_QUANT'): finish_run(32)
    abundance=out/'amalgkit_quant'/species/run/(run+'_abundance.tsv')
    abundance.parent.mkdir(parents=True,exist_ok=True)
    abundance.write_text('target_id\ttpm\n'+species+'_g1\t500000\n'+species+'_g2\t500000\n')
    finish('transcriptome_quant',[abundance])
if os.environ[prefix+'RUN_AMALGKIT_MERGE']=='1' and begin('transcriptome_merge'):
    if os.environ.get('FAKE_GG_INCOMPLETE_MERGE'): finish_run(33)
    outputs=[]
    for suffix in ('eff_length','est_counts','tpm','metadata'):
        path=out/'amalgkit_merge'/species/(species+'_'+suffix+'.tsv')
        path.parent.mkdir(parents=True,exist_ok=True); path.write_text('synthetic\n'); outputs.append(path)
    finish('transcriptome_merge',outputs)
if os.environ[prefix+'REMOVE_AMALGKIT_FASTQ_AFTER_COMPLETION']=='1':
    for p in (out/'amalgkit_getfastq').rglob('*'):
        if p.is_file() and p.name.endswith(('.fastq','.fastq.gz','.fq','.fq.gz','.sra')): p.unlink()
finish_run(0)
'''
    script = repo / "workflow/gg_transcriptome_generation_entrypoint.sh"
    script.write_text("#!/usr/bin/env bash\nexec " + sys.executable + " - <<'PY'\n" + implementation + "\nPY\n")
    api = repo / 'workflow/support/workflow_api.py'
    api.parent.mkdir(parents=True)
    api.write_text("import hashlib, json, sys\nfrom pathlib import Path\nattempt=Path(sys.argv[sys.argv.index('--attempt')+1])\ncontracts=[]\nfor p in attempt.glob('contract-*.json'):\n    r=json.loads(p.read_text()); proof=Path(r['proof'])\n    current=False\n    if proof.is_file():\n        data=json.loads(proof.read_text())\n        current=all(Path(p).is_file() and hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in data.items())\n    contracts.append({'step':r['step'],'state':'verified_current' if current else 'needs_run'})\nprint(json.dumps({'schema':'genegalleon-api-v1','command':'preflight','contracts':contracts}))\n")
    cfg = yaml.safe_load((root / "config/build.yaml").read_text())
    cfg["genegalleon"]["repository"] = str(repo)
    (root / "config/build.yaml").write_text(yaml.safe_dump(cfg))
    return repo


def new_dataset(root, names=("New plant",)):
    fake_genegalleon(root)
    fields = ["scientific_name", "run", "taxid"]
    write_tsv(root / "input/new.tsv", fields, [dict(zip(fields, [name, f"SRR{i+1}", "42"])) for i,name in enumerate(names)])
    return prepare(root, "addition", root / "config/build.yaml", "input/new.tsv")


@pytest.mark.parametrize("key", ["amalgkit_rrna_filter", "amalgkit_contam_filter"])
@pytest.mark.parametrize("value", ["yes", "no"])
def test_read_filter_settings_are_not_configurable(dataset_project, key, value):
    config = dataset_project / "config/build.yaml"
    cfg = yaml.safe_load(config.read_text())
    cfg["genegalleon"]["settings"][key] = value
    config.write_text(yaml.safe_dump(cfg))
    with pytest.raises(ValueError, match=f"managed/invalid GeneGalleon setting: {key}"):
        plan(dataset_project, config)


def test_sample_workers_reuse_and_native_array_filename_order(dataset_project, monkeypatch):
    root = dataset_project
    for key in ('AMALGKIT_RRNA_FILTER', 'AMALGKIT_CONTAM_FILTER'):
        monkeypatch.setenv('GG_TRANSCRIPTOME_' + key, 'yes')
    path = new_dataset(root, ('New plant', 'New plant alba'))
    commands = submit(path, until='quant', dry_run=True)
    assert len(commands) == 1 and '--array=1,2' in commands[0]
    assert '--stage' not in Path(commands[0][-1]).read_text()
    for index in (1, 2): worker(path, index)
    events = native_events(path)
    assert sorted(e['species'] for e in events) == ['New_plant', 'New_plant_alba']
    for event in events:
        assert event['env']['GG_TRANSCRIPTOME_AMALGKIT_RRNA_FILTER'] == 'no'
        assert event['env']['GG_TRANSCRIPTOME_AMALGKIT_CONTAM_FILTER'] == 'no'
        assert event['env']['GG_TRANSCRIPTOME_RUN_ASSEMBLY'] == '1'
        assert event['env']['GG_TRANSCRIPTOME_RUN_BUSCO_LONGEST_CDS'] == '1'
        assert event['env']['GG_TRANSCRIPTOME_RUN_AMALGKIT_QUANT'] == '1'
        assert event['env']['GG_TRANSCRIPTOME_DELETE_TMP_DIR'] == '1'
    assert all(r['quant'] == 'reuse' for r in status(path))
    worker(path, 1)
    assert len(native_events(path)) == 2
    assert submit(path, until='quant', dry_run=True) == []
    assert len(read_tsv(materialize(path) / 'metadata.tsv')) == 2


def test_failed_stage_is_not_registered_and_retry_is_limited(dataset_project, monkeypatch):
    path = new_dataset(dataset_project)
    submit(path, until='quant', dry_run=True)
    monkeypatch.setenv('FAKE_GG_FAIL_ASSEMBLY', '1')
    with pytest.raises(subprocess.CalledProcessError): worker(path, 1)
    first = status(path)[0]
    assert first['assembly'] == 'pending' and first['stopped_at'] == 'assembly'
    assert first['jobs']['sample']['state'] == 'failed'
    assert not (path / 'work/genegalleon/New_plant_SRR1/output/transcriptome_assembly/tmp').exists()
    monkeypatch.delenv('FAKE_GG_FAIL_ASSEMBLY')
    monkeypatch.setenv('FAKE_GG_FAIL_QUANT', '1')
    with pytest.raises(subprocess.CalledProcessError): worker(path, 1)
    partial = status(path)[0]
    assert [partial[s] for s in ('assembly', 'busco', 'quant')] == ['reuse', 'reuse', 'pending']
    assert partial['stopped_at'] == 'quant'
    monkeypatch.delenv('FAKE_GG_FAIL_QUANT')
    worker(path, 1)
    assert status(path)[0]['quant'] == 'reuse'
    event = native_events(path)[-1]['env']
    assert event['GG_TRANSCRIPTOME_RUN_ASSEMBLY'] == '0'
    assert event['GG_TRANSCRIPTOME_RUN_BUSCO_LONGEST_CDS'] == '0'
    assert event['GG_TRANSCRIPTOME_RUN_AMALGKIT_QUANT'] == '1'
    assert not (path / 'jobs/incomplete').exists()


def test_completed_new_species_reused_by_next_dataset(dataset_project):
    root = dataset_project
    path = new_dataset(root)
    submit(path, until="quant", dry_run=True)
    worker(path, 1)
    from database_fixtures import database_from_stages
    frozen = load(path)
    source = database_from_stages(root, frozen['config']['store'], frozen['items'], 'completed_source')
    config = root/'config/build.yaml'
    cfg = yaml.safe_load(config.read_text()); cfg['reuse_from'] = str(source)
    config.write_text(yaml.safe_dump(cfg))
    second = prepare(root, "next", config, "input/new.tsv")
    assert submit(second, until="quant", dry_run=True) == []
    assert read_tsv(materialize(second) / "metadata.tsv")[0]["run"] == "SRR1"


def test_partial_pilot_never_exports_incomplete_dataset(dataset_project):
    path = new_dataset(dataset_project, ("New plant", "Other plant"))
    subset = dataset_project / "pilot.txt"
    subset.write_text("Other_plant\n")
    commands = submit(path, until="mapping", species=subset, dry_run=True)
    assert len(commands) == 1
    assert all("--array=2" in cmd for cmd in commands)
    assert not any("mapping" in cmd[-1] for cmd in commands)
    worker(path, 2)
    with pytest.raises(ValueError, match="dataset incomplete: New_plant"):
        materialize(path)
    assert "--array=1" in submit(path, until='quant', dry_run=True)[0]


def test_slurm_submission_dependencies_and_duplicate_submission_guard(dataset_project, monkeypatch):
    path = new_dataset(dataset_project)
    calls = []
    active = False
    def scheduler(command, **kwargs):
        nonlocal active
        calls.append(command)
        if command[0] == "squeue": return "1001\n" if active else ""
        return str(1000 + len([c for c in calls if c[0] == "sbatch"])) + ";cluster\n"
    monkeypatch.setattr(subprocess, "check_output", scheduler)
    submit(path, until="mapping")
    batch = json.loads((path / "jobs/submission_0001.json").read_text())
    assert len(batch["jobs"]) == 2
    assert "--dependency=afterok:1001" in batch["jobs"][1]["command"]
    active = True
    with pytest.raises(ValueError, match="queued/running"):
        submit(path, until="mapping")
    assert len([c for c in calls if c[0] == "sbatch"]) == 2


def test_private_relative_reads_are_frozen_and_reuse_detects_changed_bytes(dataset_project):
    root = dataset_project
    fake_genegalleon(root)
    reads = root / "input/local.fastq"
    reads.write_text("@r1\nATGC\n+\nIIII\n")
    fields = ["scientific_name", "run", "taxid", "private_file", "lib_layout", "read1_path"]
    write_tsv(root / "input/private.tsv", fields,
              [dict(zip(fields, ["Private plant", "LOCAL1", "42", "yes", "single", "local.fastq"]))])
    cfg = root / "config/build.yaml"
    path = prepare(root, "private", cfg, "input/private.tsv")
    submit(path, until="quant", dry_run=True)
    staged = path / "work/genegalleon/Private_plant_LOCAL1/input/reads/Private_plant_LOCAL1/read1_path.fastq"
    assert staged.read_bytes() == reads.read_bytes()
    assert read_tsv(path / "work/genegalleon/Private_plant_LOCAL1/input/amalgkit_metadata/Private_plant_metadata.tsv")[0]["read1_path"] == "/workspace/input/reads/Private_plant_LOCAL1/read1_path.fastq"
    original = reads.read_bytes()
    reads.write_text("changed after submission")
    with pytest.raises(ValueError, match="registered file changed"):
        worker(path, 1)
    reads.write_bytes(original)
    worker(path, 1)
    from database_fixtures import database_from_stages
    frozen = load(path)
    source = database_from_stages(root, frozen['config']['store'], frozen['items'], 'private_source')
    config = yaml.safe_load(cfg.read_text()); config['reuse_from'] = str(source)
    cfg.write_text(yaml.safe_dump(config))
    assert plan(root, cfg, "input/private.tsv")[-1][0]["quant"] == "reuse"
    second = prepare(root, "private_reused", cfg, "input/private.tsv")
    assert submit(second, until="quant", dry_run=True) == []
    reads.write_text("different reads under same run")
    report = plan(root, cfg, "input/private.tsv")[-1][0]
    assert report["quant"] == "conflict"
    assert "registered file changed" in report["reason"]
    reads.unlink()
    assert plan(root, cfg, "input/private.tsv")[-1][0]["quant"] == "reuse"


def test_retry_quarantine_does_not_touch_another_species_with_same_prefix(dataset_project, monkeypatch):
    path = new_dataset(dataset_project, ("New plant", "New plant alba"))
    submit(path, until='quant', dry_run=True)
    monkeypatch.setenv("FAKE_GG_FAIL_ASSEMBLY", "1")
    with pytest.raises(subprocess.CalledProcessError): worker(path, 1)
    monkeypatch.delenv("FAKE_GG_FAIL_ASSEMBLY")
    worker(path, 2)
    other = path / "work/genegalleon/New_plant_alba_SRR2/output/transcriptome_assembly/assembled_transcripts_with_isoforms/New_plant_alba_isoform.fa.gz"
    other.parent.mkdir(parents=True, exist_ok=True)
    other.write_bytes(b"completed output of another species")
    worker(path, 1)
    assert other.read_bytes() == b"completed output of another species"
    assert all(r["assembly"] == "reuse" for r in status(path))


def test_split_slurm_arrays_preserve_species_identity_and_bound_concurrency(dataset_project, monkeypatch):
    path = new_dataset(dataset_project, ("Alpha new", "Beta new", "Gamma new"))
    monkeypatch.setattr('phase_config.resolve_array_size', lambda: 2)
    commands = submit(path, until='quant', dry_run=True)
    assert len(commands) == 2
    assert "--array=1,2" in commands[0]
    assert "--array=1" in commands[1]
    assert "--dependency=afterany:JOB_ID_sample" in commands[1]
    # Execute the generated high-index batch locally with the GeneGalleon double.
    # Local Slurm index 1 must select logical species 3, not species 1.
    subprocess.run(["bash", commands[1][-1]], cwd=dataset_project,
                   env=dict(os.environ, SLURM_ARRAY_TASK_ID="1"), check=True)
    assert [r["assembly"] for r in status(path)] == ["pending", "pending", "reuse"]
    calls = []
    def scheduler(command, **kwargs):
        calls.append(command)
        assert command[0] == "sbatch"
        return str(1000 + len(calls)) + "\n"
    monkeypatch.setattr(subprocess, "check_output", scheduler)
    submit(path, until='quant')
    batch = json.loads((path / "jobs/submission_0002.json").read_text())
    assert [j["array_offset"] for j in batch["jobs"]] == [0]
    assert [j["indices"] for j in batch["jobs"]] == [[1, 2]]
    assert not any(arg.startswith("--dependency=") for arg in calls[0])


@pytest.mark.parametrize('limits,suffix', [({}, ''), ({'jobs': 2}, '%2'), ({'cpus': 8}, '%2')])
def test_auto_array_size_preserves_stage_dependencies_and_sparse_retry(dataset_project, monkeypatch, limits, suffix):
    config = dataset_project / 'config/build.yaml'
    cfg = yaml.safe_load(config.read_text())
    cfg['slurm']['total_limits'].update(limits)
    config.write_text(yaml.safe_dump(cfg))
    path = new_dataset(dataset_project, ('Alpha new', 'Beta new', 'Gamma new'))
    frozen = (path / 'build.json').read_bytes()
    queries = []
    def detect():
        queries.append(True)
        return 3
    monkeypatch.setattr('phase_config.resolve_array_size', detect)
    commands = submit(path, dry_run=True)
    assert len(commands) == 2  # One sample array, followed by the mapping controller.
    assert '--array=1,2,3' + suffix in commands[0]
    assert not any(arg.startswith('--dependency=') for arg in commands[0])
    assert '--dependency=afterok:JOB_ID_sample' in commands[1]
    assert '--kill-on-invalid-dep=yes' in commands[1]
    # An array task failure must not change the identities of the remaining tasks.
    for task_id in (1, 2, 3):
        env = dict(os.environ, SLURM_ARRAY_TASK_ID=str(task_id))
        if task_id == 2:
            env['FAKE_GG_FAIL_ASSEMBLY'] = '1'
        result = subprocess.run(['bash', commands[0][-1]], cwd=dataset_project, env=env,
                                text=True, capture_output=True)
        assert bool(result.returncode) == (task_id == 2), result.stderr
    assert [row['assembly'] for row in status(path)] == ['reuse', 'pending', 'reuse']
    calls = []
    def scheduler(command, **kwargs):
        assert command[0] == 'sbatch'
        calls.append(command)
        return str(1000 + len(calls)) + '\n'
    monkeypatch.setattr(subprocess, 'check_output', scheduler)
    retry = submit(path)
    assert '--array=2' + suffix in retry[0]
    assert len(retry) == 2
    batch = json.loads((path / 'jobs/submission_0002.json').read_text())
    assert [job['array_offset'] for job in batch['jobs']] == [0, 0]
    assert [job['indices'] for job in batch['jobs']] == [[2], None]
    resources = json.loads((path / 'jobs/submission_0002.resources.json').read_text())
    assert 'array_size' not in resources['slurm']
    assert resources['resolved_array_size'] == 3
    assert (path / 'build.json').read_bytes() == frozen
    subprocess.run(['bash', retry[0][-1]], cwd=dataset_project,
                   env=dict(os.environ, SLURM_ARRAY_TASK_ID='2'), check=True)
    assert all(row['assembly'] == 'reuse' for row in status(path))
    assert submit(path, until='quant', dry_run=True) == []
    assert len(queries) == 2  # Once per submission; no query when no array is needed.


def test_legacy_array_sizes_in_saved_build_and_override_are_ignored(dataset_project, monkeypatch):
    path = new_dataset(dataset_project, ('Alpha new', 'Beta new', 'Gamma new'))
    # Model a frozen build prepared before array_size was retired.
    manifest = load(path)
    manifest['config']['slurm']['array_size'] = 2
    write_json(path / 'build.json', manifest)
    hashes = json.loads((path / 'checksums.json').read_text())
    hashes['build.json'] = file_record(path / 'build.json')['sha256']
    write_json(path / 'checksums.json', hashes)
    original = (path / 'build.json').read_bytes()
    resources = dataset_project / 'retry.yaml'
    resources.write_text('slurm:\n  array_size: 1\n')
    monkeypatch.setattr('phase_config.resolve_array_size', lambda: 3)
    commands = submit(path, until='quant', dry_run=True, resources=resources)
    assert len(commands) == 1
    assert '--array=1,2,3' in commands[0]
    saved = json.loads((path / 'jobs/submission_0001.resources.json').read_text())
    assert 'array_size' not in saved['slurm']
    assert saved['resolved_array_size'] == 3
    assert (path / 'build.json').read_bytes() == original
    assert load(path)['config']['slurm']['array_size'] == 2


def test_auto_array_size_preserves_sample_ids_when_cluster_limit_changes(dataset_project, monkeypatch):
    path = new_dataset(dataset_project, ('Alpha new', 'Beta new', 'Gamma new'))
    monkeypatch.setattr('phase_config.resolve_array_size', resolve_array_size)
    outputs = iter(['MaxArraySize = 4\nSchedulerParameters = (null)\n',
                    'MaxArraySize = 1001\nSchedulerParameters = max_array_tasks=2\n'])
    calls = []
    def query(command, **kwargs):
        assert command == ['scontrol', 'show', 'config']
        calls.append(command)
        return next(outputs)
    monkeypatch.setattr(subprocess, 'check_output', query)
    commands = submit(path, until='quant', dry_run=True)
    assert len(commands) == 1
    assert '--array=1,2,3' in commands[0]
    for task_id in (1, 2):
        subprocess.run(['bash', commands[0][-1]], cwd=dataset_project,
                       env=dict(os.environ, SLURM_ARRAY_TASK_ID=str(task_id)), check=True)
    # A lower detected limit changes local array indices, never sample identity.
    retry = submit(path, until='quant', dry_run=True)
    assert len(retry) == 1
    assert '--array=1' in retry[0]
    saved = json.loads((path / 'jobs/submission_0002.resources.json').read_text())
    assert 'array_size' not in saved['slurm']
    assert saved['resolved_array_size'] == 2
    subprocess.run(['bash', retry[0][-1]], cwd=dataset_project,
                   env=dict(os.environ, SLURM_ARRAY_TASK_ID='1'), check=True)
    assert all(row['assembly'] == 'reuse' for row in status(path))
    assert submit(path, until='quant', dry_run=True) == []
    assert len(calls) == 2


def test_auto_array_query_failure_prevents_submission(dataset_project, monkeypatch):
    path = new_dataset(dataset_project)
    monkeypatch.setattr('phase_config.resolve_array_size', resolve_array_size)
    calls = []
    def fail(command, **kwargs):
        calls.append(command)
        raise subprocess.CalledProcessError(1, command)
    monkeypatch.setattr(subprocess, 'check_output', fail)
    with pytest.raises(ValueError, match='Cannot query Slurm array limits'):
        submit(path)
    assert calls == [['scontrol', 'show', 'config']]
    assert not list((path / 'jobs').glob('submission_*.json'))


def test_extra_native_metadata_file_cannot_shift_species_array_index(dataset_project):
    path = new_dataset(dataset_project)
    submit(path, until='quant', dry_run=True)
    (path / "work/genegalleon/New_plant_SRR1/input/amalgkit_metadata/Aardvark_backup.tsv").write_text("unexpected metadata")
    with pytest.raises(ValueError, match="metadata file set changed"):
        worker(path, 1)
    assert not native_events(path)


def test_prepare_resolves_automatic_dependencies_once_and_freezes_the_lock(dataset_project, monkeypatch):
    import dataset_software
    root = dataset_project
    repository = fake_genegalleon(root)
    config_path = root / "config/build.yaml"
    cfg = yaml.safe_load(config_path.read_text())
    cfg["genegalleon"].pop("repository", None)
    cfg["genegalleon"].pop("image", None)
    config_path.write_text(yaml.safe_dump(cfg))
    fields = ["scientific_name", "run", "taxid"]
    write_tsv(root / "input/automatic.tsv", fields, [dict(zip(fields, ["New plant", "SRR1", "42"]))])
    calls = []
    def source(config):
        calls.append("source")
        return repository, {"kind": "downloaded_source", "identity": {"revision": config["revision"]},
                            "files": dataset_software.source_records(repository)}
    def image(config):
        calls.append("image")
        path = repository / "genegalleon.sif"
        return path, {"kind": "downloaded_image", "identity": {"uri": config["image_uri"]},
                      "files": [dataset_software.record(path)]}
    monkeypatch.setattr(dataset_software, "fetch_source", source)
    monkeypatch.setattr(dataset_software, "fetch_image", image)
    # Planning must not fetch dependencies even when every upstream stage is missing.
    assert plan(root, config_path, "input/automatic.tsv")[-1][0]["assembly"] == "pending"
    assert calls == []
    path = prepare(root, "automatic", config_path, "input/automatic.tsv")
    manifest = load(path, check_code=True)
    assert calls == ["source", "image"]
    assert manifest["config"]["genegalleon"]["repository"] == str(repository)
    assert manifest["software_lock"]["source"]["identity"]["revision"] == cfg["genegalleon"]["revision"]
    cfg["genegalleon"]["revision"] = "f" * 40
    config_path.write_text(yaml.safe_dump(cfg))
    submit(path, until='quant', dry_run=True)
    worker(path, 1)
    assert calls == ["source", "image"]
    assert load(path)["software_lock"] == manifest["software_lock"]


def test_database_controller_finishes_expression_before_publication(dataset_project, monkeypatch):
    import build_products
    import dataset
    root = dataset_project
    imported(root)
    build = prepare(root, 'controller', root / 'config/build.yaml')
    calls = []
    monkeypatch.setattr(dataset.subprocess, 'run', lambda command, **kwargs: calls.append(command))
    monkeypatch.setattr(build_products, 'complete', lambda path: calls.append(path))
    dataset.run_mapping(build)
    assert calls[0][-2:] == ['--', 'database']
    assert calls[1] == build
