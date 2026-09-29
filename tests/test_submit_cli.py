"""One-command submission freezes inputs once and preserves safe retry behavior."""
import json
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

import analysis
import dataset
from test_datasets import dataset_project, imported
from test_portable_build import completed_project


@pytest.fixture
def scheduler(monkeypatch):
    state = {'commands': [], 'running': False}
    def check_output(command, **kwargs):
        if command[0] == 'sbatch':
            state['commands'].append(command)
            return f"{1000 + len(state['commands'])};test\n"
        if command[0] == 'squeue':
            return ''.join(f'{1000 + i}\n' for i in range(1, len(state['commands']) + 1)) if state['running'] else ''
        pytest.fail(f'unexpected external command: {command}')
    monkeypatch.setattr(subprocess, 'check_output', check_output)
    return state


def invoke(monkeypatch, module, *args):
    monkeypatch.setattr(sys, 'argv', [module.__file__, *map(str, args)])
    return module.main()


def test_build_submit_prepares_and_retries_saved_inputs(dataset_project, scheduler, monkeypatch, capsys):
    root = dataset_project
    imported(root)
    config = root / 'config/build.yaml'
    cfg = yaml.safe_load(config.read_text()); cfg['name'] = 'combined'
    config.write_text(yaml.safe_dump(cfg))
    command = ('submit', '--root', root)
    invoke(monkeypatch, dataset, *command, '--dry-run')
    build = root / 'results/combined'
    before = (build / 'build.json').read_bytes()
    assert scheduler['commands'] == []
    assert 'sbatch' in capsys.readouterr().out
    assert (build / 'jobs/0001_mapping.sh').is_file()
    # Retry uses the snapshot, even when live inputs are gone or invalid.
    (root / 'input/metadata.tsv').unlink()
    cfg['translation']['table'] = 2
    config.write_text(yaml.safe_dump(cfg))
    resources = root / 'retry.yaml'
    resources.write_text('slurm:\n  stages:\n    controller:\n      mem_mb: 16000\n')
    invoke(monkeypatch, dataset, *command, '--resources', 'retry.yaml')
    assert len(scheduler['commands']) == 1
    assert '--mem=16000M' in scheduler['commands'][0]
    assert (build / 'build.json').read_bytes() == before
    assert dataset.load(build)['analysis']['translation']['table'] == 1
    assert 'Using saved build' in capsys.readouterr().out
    scheduler['running'] = True
    with pytest.raises(ValueError, match='queued/running'):
        invoke(monkeypatch, dataset, *command)
    assert len(scheduler['commands']) == 1
    # The original path-based submission still works without a source config.
    config.unlink(); scheduler['running'] = False
    invoke(monkeypatch, dataset, 'submit', '--build', build)
    assert len(scheduler['commands']) == 2


def test_build_submit_creates_and_submits_in_one_call(dataset_project, scheduler, monkeypatch):
    root = dataset_project
    imported(root)
    invoke(monkeypatch, dataset, 'submit', '--root', root, '--name', 'first', '--until', 'database')
    assert len(scheduler['commands']) == 1
    receipt = json.loads((root / 'results/first/jobs/submission_0001.json').read_text())
    assert receipt['jobs'][0]['state'] == 'submitted'
    assert (root / 'results/first/build.json').is_file()


def test_analysis_submit_prepares_and_retries_saved_conditions(completed_project, scheduler, monkeypatch, capsys):
    root, build, _, _ = completed_project
    cfg_path = root / 'config/analysis.yaml'
    cfg = yaml.safe_load(cfg_path.read_text())
    cfg['build'] = str(build)
    cfg['inputs']['species_trait'] = None
    cfg['phylogeny']['trees'] = []; cfg['phylogeny']['contrast_pairs']['enabled'] = False
    cfg_path.write_text(yaml.safe_dump(cfg))
    command = ('submit', '--root', root, '--name', 'combined')
    invoke(monkeypatch, analysis, *command)
    run = build / 'downstream/combined'
    before = (run / 'analysis.json').read_bytes()
    assert len(scheduler['commands']) == 1
    assert json.loads((run / 'jobs/submission_0001.json').read_text())['state'] == 'submitted'
    cfg['selection']['busco_threshold'] = 2  # Invalid for a new run, irrelevant to this saved run.
    cfg_path.write_text(yaml.safe_dump(cfg))
    resources = root / 'retry.yaml'; resources.write_text('slurm:\n  stages:\n    controller:\n      mem_mb: 12000\n')
    invoke(monkeypatch, analysis, *command, '--dry-run', '--resources', 'retry.yaml')
    assert len(scheduler['commands']) == 1
    assert '--mem=12000M' in capsys.readouterr().out
    assert (run / 'analysis.json').read_bytes() == before
    scheduler['running'] = True
    with pytest.raises(ValueError, match='queued/running'):
        invoke(monkeypatch, analysis, *command)
    assert len(scheduler['commands']) == 1
    cfg_path.unlink(); scheduler['running'] = False
    invoke(monkeypatch, analysis, 'submit', '--analysis', run)
    assert len(scheduler['commands']) == 2
    assert (run / 'analysis.json').read_bytes() == before


def test_analysis_submit_dry_run_prepares_without_sbatch(completed_project, scheduler, monkeypatch):
    root, build, _, _ = completed_project
    cfg = root / 'dry.yaml'
    cfg.write_text('inputs:\n  species_trait: null\nphylogeny:\n  trees: []\n  contrast_pairs:\n    enabled: false\n')
    invoke(monkeypatch, analysis, 'submit', '--root', root, '--config', cfg,
           '--build', build, '--name', 'dry', '--dry-run')
    run = build / 'downstream/dry'
    assert (run / 'analysis.json').is_file()
    assert (run / 'jobs/0001_analysis.sh').is_file()
    assert scheduler['commands'] == []
    assert not (run / 'jobs/submission_0001.json').exists()


@pytest.mark.parametrize('module,args', [
    (dataset, ['--build', 'results/existing', '--name', 'new']),
    (dataset, ['--build', 'results/existing', '--config', 'override.yaml']),
    (analysis, ['--analysis', 'results/existing/downstream/test', '--build', 'results/other']),
    (analysis, ['--analysis', 'results/existing/downstream/test', '--config', 'override.yaml']),
])
def test_submit_rejects_ambiguous_paths_and_new_settings(module, args, monkeypatch):
    with pytest.raises(SystemExit) as error:
        invoke(monkeypatch, module, 'submit', *args)
    assert error.value.code == 2


def test_named_submit_does_not_replace_unrelated_directories(dataset_project, scheduler):
    root = dataset_project
    target = root / 'results/existing'; target.mkdir(parents=True)
    marker = target / 'keep.txt'; marker.write_text('keep')
    with pytest.raises(ValueError, match='not a prepared build'):
        dataset.submit_named(root, 'existing')
    assert marker.read_text() == 'keep'
    assert scheduler['commands'] == []
