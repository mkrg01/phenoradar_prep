"""Inspect real Snakemake metadata and grouping without altering a saved run."""
import base64
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from analysis_inspection import inspect_workflow
from phase_config import merge_slurm

ROOT = Path(__file__).resolve().parents[1]


def snapshot(root):
    return {str(p.relative_to(root)): (p.read_bytes(), p.stat().st_mtime_ns, p.stat().st_ctime_ns)
            for p in root.rglob('*') if p.is_file()}


@pytest.fixture
def workflow(tmp_path):
    pytest.importorskip('snakemake')
    root = tmp_path
    (root / 'workflow').mkdir()
    (root / 'results/test').mkdir(parents=True)
    (root / 'results/test/pipeline.yaml').write_text('{}')
    snakefile = root / 'workflow/Snakefile'
    snakefile.write_text('''
rule all:
    input: "out", "local"
rule first:
    input: "input"
    output: "middle"
    group: "marker"
    params: mode="initial"
    shell: "cp {input} {output}"
rule second:
    input: "middle"
    output: "out"
    group: "marker"
    shell: "cp {input} {output}"
localrules: publish
rule publish:
    output: "local"
    shell: "echo ok > {output}"
''')
    (root / 'input').write_text('sequence\n')
    public = yaml.safe_load((ROOT / 'config/analysis.yaml').read_text())['slurm']
    return root, {'root': str(root), 'pipeline': {'output_root': 'results/test'}}, merge_slurm(public, {}, build=False)


def inspect(workflow):
    root, manifest, slurm = workflow
    before = snapshot(root)
    report = inspect_workflow(manifest, slurm, 'all')
    assert snapshot(root) == before
    return report


def execute(root):
    result = subprocess.run([sys.executable, '-m', 'snakemake', '--snakefile', 'workflow/Snakefile',
                             '--cores', '2', '--shared-fs-usage', 'persistence', 'input-output',
                             'software-deployment', 'sources', '--', 'all'],
                            cwd=root, capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr


def test_readonly_plan_counts_grouped_workers_separately(workflow):
    report = inspect(workflow)
    assert report['workflow_steps'] == 3
    assert report['worker_jobs'] == 1
    assert report['local_steps'] == 1
    assert report['pending'] == {'first': 1, 'second': 1, 'publish': 1}
    assert report['reasons'] == {'missing_output': 3}
    assert not report['partial']


def test_container_declarations_are_inspected_without_preparing_images(workflow):
    root, _, _ = workflow
    (root / 'VERSION').write_text('0.4.20\n')
    (root / 'workflow/envs').mkdir()
    (root / 'workflow/envs/test.yaml').write_text('dependencies: [python]\n')
    snakefile = root / 'workflow/Snakefile'
    snakefile.write_text('containerized: "docker://example.invalid/test:fixed"\n' +
                        snakefile.read_text().replace('params: mode="initial"',
                        'conda: "envs/test.yaml"\n    params: mode="initial"'))
    assert inspect(workflow)['worker_jobs'] == 1


@pytest.mark.parametrize('containerized', [False, True])
def test_completed_environment_metadata_uses_real_readonly_cache_paths(workflow, monkeypatch, containerized):
    root, _, _ = workflow
    execute(root)
    (root / 'VERSION').write_text('0.5.5\n')
    env = root / 'workflow/envs/test.yaml'
    env.parent.mkdir(); env.write_text('dependencies: [python]\n')
    uri = 'docker://example.invalid/test:fixed'
    snakefile = root / 'workflow/Snakefile'
    source = snakefile.read_text().replace('params: mode="initial"',
                                         'conda: "envs/test.yaml"\n    params: mode="initial"')
    if containerized:
        source = f'containerized: "{uri}"\n' + source
        cache = root / '.snakemake/singularity'
        cache.mkdir(exist_ok=True)
        (cache / (hashlib.md5(uri.encode()).hexdigest() + '.simg')).write_bytes(b'cached test image')
    snakefile.write_text(source)
    fake_bin = root / 'test_bin'; fake_bin.mkdir()
    tool = f'''#!{sys.executable}
import json,sys
from pathlib import Path
if Path(sys.argv[0]).name == 'singularity':
    if sys.argv[1:] == ['--version']:
        print('singularity version 3.11.0'); raise SystemExit(0)
    images = [Path(x) for x in sys.argv if x.endswith('.simg')]
    assert len(images) == 1 and images[0].is_file(), 'inspection missed the existing SIF cache'
    command = sys.argv[-1]
else:
    command = 'conda ' + ' '.join(sys.argv[1:])
if command == 'conda info --json':
    print(json.dumps({{'conda_prefix':'/opt/conda','platform':'linux-64'}}))
elif command == 'conda --version': print('conda 26.1.1')
elif command == 'type conda': print('conda is /opt/conda/bin/conda')
elif command == 'conda config --get channel_priority --json': print('{{"get":{{"channel_priority":"strict"}}}}')
else: raise SystemExit('inspection attempted an unexpected environment operation: '+command)
'''
    for name in ('singularity', 'conda'):
        p = fake_bin / name; p.write_text(tool); p.chmod(0o755)
    monkeypatch.setenv('PATH', str(fake_bin) + os.pathsep + os.environ['PATH'])
    env_content = env.read_bytes()
    prefix = b'' if containerized else str((root / '.snakemake/conda').resolve()).encode()
    env_hash = hashlib.md5(prefix + env_content).hexdigest()
    for output in ('middle', 'out', 'local'):
        path = root / '.snakemake/metadata' / base64.urlsafe_b64encode(output.encode()).decode()
        record = json.loads(path.read_text())
        stack = (env_hash.encode() if output == 'middle' else b'') + (uri.encode() if containerized else b'')
        record['software_stack_hash'] = hashlib.md5(stack).hexdigest()
        path.write_text(json.dumps(record))
    report = inspect(workflow)
    assert report['pending'] == {}
    assert sum(report['reused'].values()) == 3


def test_pending_rule_resources_are_checked_against_worker_budget(workflow):
    root, manifest, slurm = workflow
    scripts = ROOT / 'workflow/scripts'
    snakefile = root / 'workflow/Snakefile'
    snakefile.write_text(f'import sys\nsys.path.insert(0, {str(scripts)!r})\n' + snakefile.read_text() +
                        '\nfrom resource_limits import apply_workflow_limits\napply_workflow_limits(workflow)\n')
    limits = merge_slurm(slurm, {'total_limits': {'cpus': 5},
                        'per_job_resources': {'first': {'cpus': 8}}}, build=False)
    with pytest.raises(ValueError, match='remaining.*total_limits.cpus budget'):
        inspect_workflow(manifest, limits, 'all')


def test_resume_reads_real_metadata_and_explains_changed_parameters(workflow):
    root, _, _ = workflow
    execute(root)
    report = inspect(workflow)
    assert report['workflow_steps'] == report['worker_jobs'] == 0
    assert sum(report['reused'].values()) == 3
    snakefile = root / 'workflow/Snakefile'
    snakefile.write_text(snakefile.read_text().replace('mode="initial"', 'mode="changed"'))
    report = inspect(workflow)
    assert report['pending'] == {'first': 1, 'second': 1}
    assert report['worker_jobs'] == 1
    assert report['reasons']['conditions_or_implementation_changed'] == 1


def test_incomplete_and_missing_stages_resume_without_touching_metadata(workflow):
    root, _, _ = workflow
    execute(root)
    (root / 'out').unlink()
    report = inspect(workflow)
    assert report['pending'] == {'second': 1}
    (root / 'out').write_text('unfinished\n')
    incomplete = root / '.snakemake/incomplete' / base64.urlsafe_b64encode(b'out').decode()
    incomplete.write_text(json.dumps({'external_jobid': '123'}))
    report = inspect(workflow)
    assert report['pending'] == {'second': 1}
    assert report['reasons']['incomplete_output'] == 1
    assert incomplete.is_file()


def test_checkpoint_counts_are_explicitly_partial(workflow):
    root, _, _ = workflow
    (root / 'workflow/Snakefile').write_text('''
def dynamic(wc):
    plan = checkpoints.discover.get().output[0]
    return [line.strip() + ".out" for line in open(plan)]
rule all:
    input: dynamic
checkpoint discover:
    output: "plan"
    shell: "echo sample > {output}"
rule future:
    output: "{sample}.out"
    shell: "echo result > {output}"
''')
    report = inspect(workflow)
    assert report['partial']
    assert report['pending_checkpoints'] == ['discover']
    assert report['workflow_steps'] == report['worker_jobs'] == 1
    assert 'future' not in report['pending']


def test_invalid_batch_receipt_is_planned_for_repair_without_deleting_it(workflow):
    root, _, _ = workflow
    scripts = ROOT / 'workflow/scripts'
    (root / 'plan').mkdir()
    (root / 'alignments/batches').mkdir(parents=True)
    (root / 'plan/index.json').write_text(json.dumps({'batches': ['0001']}))
    receipt = root / 'alignments/batches/0001.json'
    receipt.write_text('{"members": {}}')
    (root / 'workflow/Snakefile').write_text(f'''
import sys
sys.path.insert(0, {str(scripts)!r})
from alignment_batches import batch_targets
rule all:
    input: lambda wc: batch_targets("plan", "alignments", "reports")
rule batch:
    output: "alignments/batches/{{batch}}.json"
    shell: "echo result > {{output}}"
''')
    report = inspect(workflow)
    assert report['pending'] == {'batch': 1}
    assert report['invalid_batch_receipts'] == ['alignments/batches/0001.json']
    assert report['reasons']['invalid_batch_receipt'] == 1


def test_runtime_observer_uses_real_snakemake_progress_events(workflow, monkeypatch):
    from analysis_runtime import RunState, read_state
    root, _, _ = workflow
    state = RunState(root / 'results/test', 'all')
    state.phase('computing')
    monkeypatch.setenv('PHENORADAR_RUN_STATE', str(state.path))
    monkeypatch.setenv('PHENORADAR_RUN_TOKEN', state.token)
    snakefile = root / 'workflow/Snakefile'
    snakefile.write_text(f'import sys\nsys.path.insert(0, {str(ROOT / "workflow/scripts")!r})\n'
                        'from analysis_runtime import install_progress\ninstall_progress(workflow)\n' + snakefile.read_text())
    execute(root)
    progress = read_state(state.path)['progress']
    assert progress['done'] == progress['total'] == 4
