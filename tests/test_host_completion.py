"""Exercise final tree validation in the documented host environment."""
import json
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

import analysis
from analysis_runtime import read_state
from analysis_validation import check_dependencies, validate_tree
from common import write_json
from dataset_assets import verify
from test_analysis_validation import completed_tree, snapshot

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def saved_tree(snapshot, monkeypatch):
    pipeline = yaml.safe_load((ROOT / 'workflow/pipeline_defaults.yaml').read_text())
    pipeline.update(run_name='completed', output_root='completed', trait='C4')
    pipeline['phylogeny']['trees'] = ['all']
    for option in ('contrast_pairs', 'dating', 'taxonomy_check'):
        pipeline['phylogeny'][option]['enabled'] = False
    pipeline['kegg']['enabled'] = False
    slurm = yaml.safe_load((ROOT / 'config/analysis.yaml').read_text())['execution']
    manifest = {'root': str(snapshot.parent), 'pipeline': pipeline, 'config': {'slurm': slurm},
                'build': {'path': 'database/manifest.json'}}
    write_json(snapshot / 'build_inputs/index.json', {})
    write_json(snapshot / 'run.json', {})
    write_json(snapshot / 'metadata/selection.json', {})
    completed_tree(snapshot, '((Plant_A:0.1,Plant_B-x:0.2):0.1,Plant_C:0.3);\n')
    monkeypatch.setattr(analysis, 'load', lambda *args, **kwargs: manifest)
    return snapshot, manifest


@pytest.mark.parametrize('invalid_tree', [False, True])
def test_host_controller_records_completion_only_after_real_tree_validation(saved_tree, monkeypatch, invalid_tree):
    path, manifest = saved_tree
    tree = path / 'phylogeny/all/species_tree.nwk'
    if invalid_tree:
        tree.write_text('(Plant_A:1,Unknown_plant:1);\n')
    protected = {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in path.rglob('*') if p.is_file()}
    calls = []

    def reused_workflow(command, **kwargs):
        assert command[0] == str(path.parent / 'run_pipeline.sh')
        assert command[-1] == 'all'
        assert not (path / 'completed.json').exists()
        calls.append(command)
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr(subprocess, 'run', reused_workflow)
    monkeypatch.setenv('SLURM_JOB_ID', '987654')
    if invalid_tree:
        with pytest.raises(ValueError, match='tree species differ'):
            analysis.run(path, local=True)
        assert not (path / 'completed.json').exists()
        assert read_state(path / 'jobs/run_state.json')['phase'] == 'failed'
    else:
        analysis.run(path, local=True)
        done = json.loads((path / 'completed.json').read_text())
        runtime = read_state(path / 'jobs/run_state.json')
        assert runtime['phase'] == 'complete'
        assert done['controller_id'] == '987654'
        assert done['run_token'] == runtime['run_token']
        assert str(tree) in {entry['path'] for entry in done['files']}
        for entry in done['files']:
            verify(entry)
    assert len(calls) == 1
    assert all((p.read_bytes(), p.stat().st_mtime_ns) == expected for p, expected in protected.items())


@pytest.mark.parametrize('entry', ['preflight', 'run'])
def test_missing_host_ete4_fails_before_workflow_or_submission(saved_tree, monkeypatch, entry):
    path, _ = saved_tree
    monkeypatch.setitem(sys.modules, 'ete4', None)
    monkeypatch.setattr(subprocess, 'run', lambda *args, **kwargs: pytest.fail('started workflow'))
    monkeypatch.setattr(analysis, 'submission_guard', lambda *args: pytest.fail('queried Slurm'))
    with pytest.raises(ValueError, match='ete4 in the host Python environment'):
        if entry == 'preflight':
            analysis.preflight(path, 'all', None)
        else:
            analysis.run(path, local=True)
    assert not (path / 'jobs').exists()
    assert not (path / 'completed.json').exists()


@pytest.mark.parametrize('pipeline,target', [({'phylogeny': {'trees': []}}, 'all'),
                                           ({'phylogeny': {'trees': ['all']}}, 'alignments')])
def test_no_tree_completion_check_does_not_require_ete4(monkeypatch, pipeline, target):
    monkeypatch.setitem(sys.modules, 'ete4', None)
    check_dependencies(pipeline, target)


def test_explicit_tree_verification_explains_missing_host_dependency(tmp_path, monkeypatch):
    monkeypatch.setitem(sys.modules, 'ete4', None)
    with pytest.raises(ValueError, match='update the active environment from environment.yaml'):
        validate_tree(tmp_path / 'tree.nwk', {'A'})


def test_documented_host_environment_includes_tree_validation_dependency():
    dependencies = yaml.safe_load((ROOT / 'environment.yaml').read_text())['dependencies']
    assert 'ete4=4.3.0' in dependencies
