"""Exercise immutable imported references through the real Snakemake scheduler."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from common import sha256, write_json, write_tsv
from incremental_odb import plan
from odb_settings import pin_reference
from prepare_odb_reference import prepare

ROOT = Path(__file__).resolve().parents[1]


def test_forced_reference_preparation_preserves_imported_snapshot_and_pin(
        workflow_project, frozen_reference, fake_odb, command_environment):
    snakemake = os.environ.get('SNAKEMAKE_BIN') or shutil.which('snakemake')
    if not snakemake:
        pytest.skip('Snakemake required')
    reference = workflow_project / 'resources/orthodb/v12_3193'
    shutil.copytree(frozen_reference, reference)
    files = [reference / name for name in ('reference.json', 'files.json', 'dbinfo.txt', 'config.txt')]
    before = {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in files}
    pin = workflow_project / 'work/run001/odb_reference_pin.json'
    pin_reference(pin, 'v12', 3193, sha256(reference / 'reference.json'))
    env = command_environment({'python': sys.executable, 'ODB-mapper': fake_odb})
    command = [snakemake, '--snakefile', str(ROOT / 'workflow/Snakefile'), '--cores', '1',
               '--forcerun', 'prepare_odb_reference', '--', 'references']
    for _ in range(2):
        result = subprocess.run(command, cwd=workflow_project, env=env,
                                text=True, capture_output=True, timeout=90)
        assert result.returncode == 0, result.stdout + result.stderr
        assert {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in files} == before
        assert not (reference / 'setup_identity.json').exists()
    assert json.loads(pin.read_text())['reference_sha256'] == sha256(reference / 'reference.json')
    receipt = workflow_project / 'work/run001/orthogroups/mapping/reference-ready.json'
    assert json.loads(receipt.read_text())['reference']['sha256'] == sha256(reference / 'reference.json')


def test_plan_defers_reference_identity_until_preparation(tmp_path, frozen_reference):
    proteins = tmp_path / 'proteins'
    proteins.mkdir()
    (proteins / 'Alpha_protein.fa').write_text('>Alpha_g1\nMK\n')
    samples = tmp_path / 'samples.tsv'
    write_tsv(samples, ['species', 'odb_species'], [{'species': 'Alpha', 'odb_species': 'Alpha'}])
    pin = tmp_path / 'pin.json'
    result = plan(samples, proteins, tmp_path / 'plan', tmp_path / 'cache',
                  reference=frozen_reference / 'reference.json', reference_pin=pin)
    assert result['mapped_species'] == ['Alpha']
    assert json.loads(pin.read_text())['reference_sha256'] is None


def test_invalid_imported_inventory_is_rejected_without_rewriting(frozen_reference):
    marker = frozen_reference / 'reference.json'
    before = marker.read_bytes()
    (frozen_reference / 'files.json').write_text('[]\n')
    with pytest.raises(ValueError, match='inventory'):
        prepare(frozen_reference, command='/must/not/run')
    assert marker.read_bytes() == before


def test_changed_pinned_marker_fails_before_download_or_receipt(tmp_path, frozen_reference):
    pin = tmp_path / 'pin.json'
    pin_reference(pin, 'v12', 3193, sha256(frozen_reference / 'reference.json'))
    marker = frozen_reference / 'reference.json'
    data = json.loads(marker.read_text())
    write_json(marker, dict(data, changed=True))
    changed = marker.read_bytes()
    receipt = tmp_path / 'ready.json'
    with pytest.raises(ValueError, match='pinned build snapshot'):
        prepare(frozen_reference, command='/must/not/run', receipt=receipt, reference_pin=pin)
    assert not receipt.exists()
    assert marker.read_bytes() == changed


def test_new_download_is_validated_before_pin_and_receipt(tmp_path, fake_odb):
    root, pin, receipt = tmp_path / 'reference', tmp_path / 'pin.json', tmp_path / 'ready.json'
    result = prepare(root, command=str(fake_odb), receipt=receipt, reference_pin=pin)
    assert result['files_checked'] == 1
    assert json.loads(receipt.read_text()) == result
    assert json.loads(pin.read_text())['reference_sha256'] == sha256(root / 'reference.json')


def test_missing_pinned_reference_cannot_be_replaced_by_download(tmp_path, fake_odb):
    root, pin, receipt = tmp_path / 'reference', tmp_path / 'pin.json', tmp_path / 'ready.json'
    pin_reference(pin, 'v12', 3193, 'a' * 64)
    with pytest.raises(ValueError, match='pinned reference is missing'):
        prepare(root, command=str(fake_odb), receipt=receipt, reference_pin=pin)
    assert not root.exists()
    assert not receipt.exists()
    assert json.loads(pin.read_text())['reference_sha256'] == 'a' * 64


def test_missing_registered_data_cannot_be_pinned(tmp_path, frozen_reference):
    (frozen_reference / 'odbmapper/v12/data/tiny.db').unlink()
    pin, receipt = tmp_path / 'pin.json', tmp_path / 'ready.json'
    with pytest.raises(ValueError, match='reference file missing'):
        prepare(frozen_reference, command='/must/not/run', receipt=receipt, reference_pin=pin)
    assert json.loads(pin.read_text())['reference_sha256'] is None
    assert not receipt.exists()
