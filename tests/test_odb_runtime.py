"""Real interpreter execution, native import checks, and mapping fail-fast."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from check_odb_runtime import check
from odb_environment import software_records

ROOT = Path(__file__).resolve().parents[1]
DEPLOY = ROOT / "workflow/envs/odb.post-deploy.sh"


def test_deployment_repairs_entry_points_preserves_code_and_records_originals(tmp_path):
    prefix = tmp_path / "env"
    (prefix / "bin").mkdir(parents=True)
    (prefix / "bin/python").symlink_to(sys.executable)
    helpers = prefix / "orthologer/bin"
    helpers.mkdir(parents=True)
    (helpers / "dependency.py").write_text('VALUE = "native helper works"\n')
    body = b'from dependency import VALUE\nprint(VALUE)\n'
    scripts = []
    for name, shebang in [("check_fasta.py", b"#!/missing/bin/python3"),
                          ("other.py", b"#!/bin/python"),
                          ("env.py", b"#!/usr/bin/env python3")]:
        path = helpers / name
        path.write_bytes(shebang + b'\n' + body)
        path.chmod(0o751)
        scripts.append(path)
    with pytest.raises(FileNotFoundError):
        subprocess.run([str(scripts[0])], check=True)
    env = {**os.environ, "CONDA_PREFIX": str(prefix)}
    subprocess.run(["bash", str(DEPLOY)], env=env, check=True)
    record = prefix / "share/phenoradar/odb-runtime.json"
    before = record.read_bytes()
    for path in scripts:
        assert path.read_bytes().partition(b'\n')[2] == body
        assert path.stat().st_mode & 0o777 == 0o751
        assert subprocess.check_output([str(path)], text=True).strip() == 'native helper works'
    data = json.loads(before)
    assert len(data['scripts']) == 3
    assert data['scripts'][0]['original_shebang'] == '#!/missing/bin/python3'
    subprocess.run(["bash", str(DEPLOY)], env=env, check=True)
    assert record.read_bytes() == before
    command = prefix / "bin/ODB-mapper"
    command.write_text("test command")
    records = software_records(str(command), str(prefix))
    assert any(entry['path'] == str(record) for entry in records)


def test_preflight_executes_helper_and_verifies_imported_content(tmp_path, fake_odb):
    result = check(command=str(fake_odb), work_dir=tmp_path)
    assert result['checks'] == {'fasta_type': 'PROT', 'fasta_import': 'passed'}
    assert len(result['commands']) == 3
    assert not list(tmp_path.glob('odb-runtime-*'))


def test_preflight_rejects_missing_interpreter_even_when_mapper_exists(tmp_path, fake_odb):
    helper = fake_odb.parent.parent / 'orthologer/bin/check_fasta.py'
    helper.write_text('#!/missing/bin/python3\nprint("PROT")\n')
    with pytest.raises(RuntimeError, match='Orthologer runtime check failed.*check_fasta.py'):
        check(command=str(fake_odb), work_dir=tmp_path)


def test_preflight_rejects_silent_failed_import(tmp_path, fake_odb):
    native = fake_odb.parent / 'orthologer'
    native.write_text(native.read_text().replace(
        "(project / 'Sequences' / (label + '.faa')).write_text(data)", "pass"))
    with pytest.raises(ValueError, match='did not import the expected protein'):
        check(command=str(fake_odb), work_dir=tmp_path)


def test_timeout_retains_partial_native_output(tmp_path, fake_odb, monkeypatch):
    def timeout(argv, **kwargs):
        raise subprocess.TimeoutExpired(argv, 60, output=b"native partial output")
    monkeypatch.setattr('check_odb_runtime.subprocess.run', timeout)
    with pytest.raises(RuntimeError, match='native partial output'):
        check(command=str(fake_odb), work_dir=tmp_path)


def test_failed_preflight_blocks_reference_and_mapping(
        workflow_project, tiny_inputs, seed_taxonomy, fake_odb, command_environment, tmp_path):
    snakemake = os.environ.get('SNAKEMAKE_BIN') or shutil.which('snakemake')
    if not snakemake:
        pytest.skip('Snakemake required')
    seed_taxonomy(tiny_inputs['taxonomy_db'])
    helper = fake_odb.parent.parent / 'orthologer/bin/check_fasta.py'
    helper.write_text('#!/missing/bin/python3\nprint("PROT")\n')
    env = command_environment({'python': sys.executable, 'ODB-mapper': fake_odb})
    config = tmp_path / 'runtime.yaml'
    config.write_text('run_name: runtime\n')
    result = subprocess.run([snakemake, '--snakefile', str(ROOT / 'workflow/Snakefile'),
        '--configfile', str(config), '--cores', '2', '--',
        'references'], cwd=workflow_project,
        env=env, text=True, capture_output=True, timeout=90)
    assert result.returncode != 0
    log = workflow_project / 'logs/runtime/orthogroups/mapping/reference.log'
    assert 'check_fasta.py' in log.read_text()
    assert not (workflow_project / 'resources/orthodb/v12_3193/odbmapper').exists()
    assert not list((workflow_project / 'work/runtime').rglob('status.json'))
