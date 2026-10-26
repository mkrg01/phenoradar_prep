import json
from pathlib import Path
import sys
import shutil
import subprocess
import os
import yaml

import pytest

from align_orthogroups import align, finish
from alignment_batches import plan, run_batch, batch_targets
from common import write_json


def fixture(tmp_path, count=7):
    inputs = tmp_path / 'inputs'; inputs.mkdir()
    for i in range(count):
        (inputs / f'OG{i}.faa').write_text(f'>Plant_g{i}\nMAA\n')
    write_json(inputs / 'provenance.json', {'orthogroups': [f'OG{i}' for i in range(count)]})
    return inputs, tmp_path / 'plan', tmp_path / 'out', tmp_path / 'reports'


def test_weighted_plan_bounds_and_singleton_recovery(tmp_path):
    inputs, plans, out, reports = fixture(tmp_path)
    plan(inputs, plans, batch_size=3)
    index = json.loads((plans / 'index.json').read_text())
    groups = [json.loads((plans / f'{b}.json').read_text())['orthogroups'] for b in index['batches']]
    assert len(groups) == 3 and all(len(g) <= 3 for g in groups)
    assert sorted(og for g in groups for og in g) == [f'OG{i}' for i in range(7)]
    for b in index['batches']:
        run_batch(plans / f'{b}.json', inputs, out, reports, out / 'batches' / f'{b}.json', cpus=2)
    finish(inputs, out, reports)
    times = {p: p.stat().st_mtime_ns for p in out.glob('*.faa')}
    (out / 'OG0.faa').unlink()
    targets = batch_targets(plans, out, reports)
    missing = [Path(p) for p in targets if not Path(p).exists()]
    assert len(missing) == 1
    name = missing[0].stem
    run_batch(plans / f'{name}.json', inputs, out, reports, missing[0], cpus=2)
    assert json.loads(missing[0].read_text())['computed'] == 1
    assert all(p.stat().st_mtime_ns == t for p, t in times.items() if p.name != 'OG0.faa')
    assert all(Path(p).exists() for p in batch_targets(plans, out, reports))


def test_interrupted_batch_keeps_completed_og_and_reuses_legacy_report(tmp_path):
    inputs, plans, out, reports = fixture(tmp_path, 2)
    (inputs / 'OG1.faa').write_text('>Plant_g1\nMAA\n>Plant_g2\nMA\n')
    align(inputs / 'OG0.faa', out / 'OG0.faa', reports / 'OG0.json')
    saved = (out / 'OG0.faa').stat().st_mtime_ns
    plan(inputs, plans)
    batch = plans / 'batch_000000.json'; receipt = out / 'batches/batch_000000.json'
    with pytest.raises(FileNotFoundError):
        run_batch(batch, inputs, out, reports, receipt, cpus=1, command='no-such-famsa')
    assert not receipt.exists() and (out / 'OG0.faa').stat().st_mtime_ns == saved
    command = tmp_path / 'famsa'
    command.write_text(f'#!{sys.executable}\nfrom pathlib import Path\nimport sys\n'
                       'Path(sys.argv[-1]).write_text(">Plant_g1\\nMAA\\n>Plant_g2\\nMA-\\n")\n')
    command.chmod(0o755)
    run_batch(batch, inputs, out, reports, receipt, cpus=1, command=str(command))
    assert json.loads(receipt.read_text())['reused'] == 1
    assert (out / 'OG0.faa').stat().st_mtime_ns == saved
    # Changed residues invalidate reuse even when the file size is unchanged.
    (out / 'OG1.faa').write_text('>Plant_g1\nMGA\n>Plant_g2\nMA-\n')
    assert not Path(batch_targets(plans, out, reports)[0]).exists()
    run_batch(batch, inputs, out, reports, receipt, cpus=1, command=str(command))
    assert 'MGA' not in (out / 'OG1.faa').read_text()


def test_restricted_batch_worker_does_not_open_database(workflow_project, tmp_path):
    snakemake = os.environ.get('SNAKEMAKE_BIN') or shutil.which('snakemake')
    if not snakemake:
        pytest.skip('Snakemake required')
    inputs = workflow_project / 'work/test/orthogroups/alignments/inputs'
    inputs.mkdir(parents=True)
    (inputs / 'OG1.faa').write_text('>Plant_g1\nMAA\n')
    write_json(inputs / 'provenance.json', {'orthogroups': ['OG1']})
    plan(inputs, inputs.parent / 'batches')
    database = tmp_path / 'unreadable-database.json'
    database.write_text('not database JSON; unrelated workers must never parse this')
    config = tmp_path / 'batch.yaml'
    config.write_text(yaml.safe_dump({'run_name': 'test', 'build_manifest': str(database)}))
    root = Path(__file__).resolve().parents[1]
    result = subprocess.run([snakemake, '--snakefile', str(root / 'workflow/Snakefile'),
                             '--configfile', str(config), '--cores', '12', '--dry-run',
                             '--allowed-rules', 'align_orthogroup_batch', '--target-jobs',
                             'align_orthogroup_batch:batch=batch_000000'],
                            cwd=workflow_project, text=True, capture_output=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'align_orthogroup_batch' in result.stdout
