"""Exercise numerical-library initialization in fresh workflow/tool processes."""
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

import pytest

import angiosperms353 as qc
from dataset_assets import record
import taxonomy_screen_jobs as jobs


LIBRARY_THREADS = ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
                   'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')


@pytest.mark.parametrize('entry', ['metadata', 'screen', 'batch', 'finalize'])
def test_fresh_taxonomy_entrypoints_bound_library_pools(tmp_path, entry):
    scripts = Path(jobs.__file__).parent
    # An inherited setting must not expand to the CPUs of the whole node.
    env = dict(os.environ, **dict.fromkeys(LIBRARY_THREADS, '64'))
    env['OMP_THREAD_LIMIT'] = '128'
    program = r'''
import json, os, pathlib, sys
sys.path.insert(0, sys.argv[1])
import metadata_catalog as catalog
import taxonomy_screen_jobs as jobs
assert 'numpy' not in sys.modules, 'numerical initialization precedes the runtime guard'
root = pathlib.Path(sys.argv[2])
entry = sys.argv[3]
def probe(*args, **kwargs):
    import numpy
    return dict(libraries={name: os.environ.get(name) for name in (
        'OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
        'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')},
        thread_limit=os.environ.get('OMP_THREAD_LIMIT'),
        native_threads=len(list(pathlib.Path('/proc/self/task').iterdir()))
            if pathlib.Path('/proc/self/task').is_dir() else None)
plan = dict(settings={'batch_size': 24}, samples=[{'run': 'ERR1'}])
if entry == 'metadata':
    import taxonomy_screen_filter
    catalog.configuration = lambda *args, **kwargs: {'taxonomy_screen': {'enabled': True}}
    taxonomy_screen_filter.filter_database = probe
    sys.argv = ['metadata', 'filter-taxonomy', '--root', str(root), '--database', 'db',
                '--screen', 'screen.json', '--output', 'out']
    catalog.main()
elif entry == 'screen':
    jobs.prepare_plan = lambda *args, **kwargs: plan
    jobs.batch_result = lambda *args: None
    def submit(output, plan, missing, execution, job, ram, storage):
        return dict(probe(), cpus=job['cpus'], ram=ram)
    jobs.submit = submit
    result = jobs.screen([{'run': 'ERR1'}], cache=root, taxonomy=root / 'tax.sqlite',
        output=root / 'screen', cluster=True, threads=8, ram=24,
        config={'execution': {'per_job_resources': {'taxonomy_screen': {'cpus': 8, 'mem_gb': 32}}}})
    print(json.dumps(result))
else:
    (root / 'plan.json').write_text(json.dumps(plan))
    if entry == 'batch':
        spec = root / 'task.json'
        spec.write_text(json.dumps(dict(indices=[1], threads=8, ram=24, storage={})))
        os.environ['SGE_TASK_ID'] = '1'
        os.environ.pop('SLURM_ARRAY_TASK_ID', None)
        def run_batch(output, plan, index, samples, threads, ram, storage):
            print(json.dumps(dict(probe(), cpus=threads, ram=ram)))
        jobs.run_batch = run_batch
        sys.argv = ['worker', str(root), 'batch', str(spec)]
    else:
        jobs.finalize = probe
        sys.argv = ['worker', str(root), 'finalize']
    jobs.main()
'''
    completed = subprocess.run([sys.executable, '-c', program, str(scripts), str(tmp_path), entry],
                               env=env, capture_output=True, text=True, check=True)
    result = json.loads(completed.stdout)
    assert result['libraries'] == dict.fromkeys(LIBRARY_THREADS, '1')
    # Do not set OMP_THREAD_LIMIT=1, which would override explicit tool requests.
    assert result['thread_limit'] == '128'
    if result['native_threads'] is not None:
        # Matplotlib may briefly retain its font-cache warning timer thread.
        assert result['native_threads'] <= 2
    if entry in {'screen', 'batch'}:
        assert (result['cpus'], result['ram']) == (8, 24)


def test_captus_child_is_bounded_without_changing_caller_or_tool_resources(tmp_path, monkeypatch):
    for name in LIBRARY_THREADS:
        monkeypatch.setenv(name, '64')
    monkeypatch.setenv('OMP_THREAD_LIMIT', '128')
    monkeypatch.setenv('TAXONOMY_TEST_VALUE', 'preserved')
    executable_dir = tmp_path / 'bin'
    executable_dir.mkdir()
    version_env = tmp_path / 'version.json'
    program = tmp_path / 'captus.py'
    program.write_text('''import json, os, pathlib, sys
import numpy
payload = dict(environment={name: os.environ.get(name) for name in ''' + repr(LIBRARY_THREADS) + '''},
               preserved=os.environ.get('TAXONOMY_TEST_VALUE'), tmpdir=os.environ.get('TMPDIR'))
if '--version' in sys.argv:
    pathlib.Path(''' + repr(str(version_env)) + ''').write_text(json.dumps(payload))
    print('Captus test')
else:
    args=sys.argv
    payload.update(threads=args[args.index('--threads')+1], ram=args[args.index('--ram')+1])
    target=pathlib.Path(args[args.index('-o')+1])/'sample0__captus-ext/01_coding_NUC'
    target.mkdir(parents=True)
    (target/'NUC_coding_NT.fna').write_text('>sample__1\\nATGAAATAA\\n')
    print(json.dumps(payload))
''')
    captus = executable_dir / 'captus'
    captus.write_text('#!/bin/sh\nexec ' + shlex.join([sys.executable, str(program)]) + ' "$@"\n')
    captus.chmod(0o755)
    mmseqs = executable_dir / 'mmseqs'
    mmseqs.write_text('#!/bin/sh\nprintf "MMseqs test\\n"\n')
    mmseqs.chmod(0o755)
    monkeypatch.setenv('PATH', str(executable_dir) + os.pathsep + os.environ['PATH'])
    assert jobs.tool_versions() == {'captus': 'Captus test', 'mmseqs': 'MMseqs test'}
    assert json.loads(version_env.read_text())['environment'] == dict.fromkeys(LIBRARY_THREADS, '1')
    cds = tmp_path / 'input.fna'
    cds.write_text('>gene\nATGAAATAA\n')
    work = tmp_path / 'work'
    work.mkdir()
    paths = qc.extract([dict(run='ERR1', cds=record(cds))], work,
                       threads=8, ram=24, blat='blat', log_dir=tmp_path)
    assert len(paths) == 1 and paths[0].is_file()
    result = json.loads((tmp_path / 'captus.log').read_text())
    assert result['environment'] == dict.fromkeys(LIBRARY_THREADS, '1')
    assert (result['threads'], result['ram'], result['preserved']) == ('8', '24', 'preserved')
    assert result['tmpdir'] == str(work)
    assert {name: os.environ[name] for name in LIBRARY_THREADS} == dict.fromkeys(LIBRARY_THREADS, '64')
