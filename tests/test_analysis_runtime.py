"""Completion follows native validation, without recreating an input bundle."""
import json
import subprocess
import sys

import pytest

import analysis
from analysis_runtime import read_state
from dataset_assets import verify
from test_datasets import dataset_project
from test_portable_build import completed_project


@pytest.mark.parametrize('invalid_output', [False, True])
def test_real_workflow_completion_and_validation_failure(completed_project, monkeypatch, invalid_output):
    root, build, environment, _ = completed_project
    config = root / 'minimal.yaml'
    config.write_text('inputs:\n  species_trait: null\nphylogeny:\n  trees: []\n  contrast_pairs:\n    enabled: false\n')
    path = analysis.prepare(root, 'runtime', config, build)
    phases = []
    execute = subprocess.run

    def without_container(command, **kwargs):
        # Exercise the production DAG locally with the fixture's
        # frozen taxonomy/ODB. Only the container launcher is replaced.
        assert command[0] == str(root / 'run_pipeline.sh')
        target = command[-1]
        phases.append(read_state(path / 'jobs/run_state.json')['phase'])
        assert not (path / 'completed.json').exists()
        result = execute([sys.executable, '-m', 'snakemake', '--snakefile', str(root / 'workflow/Snakefile'),
                         '--configfile', str(path / 'pipeline.yaml'), '--cores', '2',
                         '--resources', 'mem_mb=16000', '--shared-fs-usage', 'persistence', 'input-output',
                         'software-deployment', 'sources', '--', target],
                         cwd=root, env={**environment, **kwargs['env']}, capture_output=True, text=True)
        assert result.returncode == 0, result.stderr
        if invalid_output:
            table = path / 'orthogroups/expression/tpm.tsv'
            lines = table.read_text().splitlines()
            values = lines[1].split('\t')
            values[lines[0].split('\t').index('tpm')] = 'NaN'
            lines[1] = '\t'.join(values)
            table.write_text('\n'.join(lines) + '\n')
        return result

    monkeypatch.setattr(subprocess, 'run', without_container)
    if invalid_output:
        with pytest.raises(ValueError):
            analysis.run(path, local=True, cores=2, mem_mb=16000)
        assert read_state(path / 'jobs/run_state.json')['phase'] == 'failed'
        assert not (path / 'completed.json').exists()
    else:
        analysis.run(path, local=True, cores=2, mem_mb=16000)
        state = read_state(path / 'jobs/run_state.json')
        done = json.loads((path / 'completed.json').read_text())
        assert state['phase'] == 'complete'
        assert done['run_token'] == state['run_token']
        assert (path / 'orthogroups/expression/tpm.tsv').is_file()
        assert not (path / 'phenoradar_inputs').exists()
        assert done['output_layout'] == 'native'
        assert done['files']
        for entry in done['files']:
            verify(entry)
        assert analysis.status(path)['state'] == 'complete'
        checked = execute([sys.executable, str(root / 'workflow/scripts/analysis.py'),
                           'status', '--analysis', str(path), '--verify', '--json'],
                          cwd=root, env=environment, capture_output=True, text=True)
        assert checked.returncode == 0, checked.stderr
        assert json.loads(checked.stdout)['verification']['performed']
        assert not (path / 'phenoradar_inputs').exists()
    assert phases == ['computing']
