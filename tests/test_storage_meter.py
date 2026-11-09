"""Storage accounting includes inode quotas and never follows external symlinks."""
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from storage_meter import check_budget, measure, monitored, project_usage


def test_inventory_deduplicates_links_and_counts_directories(tmp_path):
    root = tmp_path / 'project'
    root.mkdir()
    (root / 'sub').mkdir()
    payload = root / 'payload'
    payload.write_bytes(b'x' * 10000)
    os.link(payload, root / 'sub/linked')
    external = tmp_path / 'external'
    external.mkdir()
    (external / 'large').write_bytes(b'y' * 1000000)
    (root / 'outside').symlink_to(external, target_is_directory=True)
    result = measure([root, root / 'sub'])
    assert not result['errors']
    assert result['inodes'] == 4
    assert result['files'] == 1 and result['directories'] == 2 and result['symlinks'] == 1
    assert result['logical_bytes'] == 10000
    assert result['allocated_bytes'] == sum(p.lstat().st_blocks * 512 for p in
                                            (root, root / 'sub', payload, root / 'outside'))


def test_missing_and_unreadable_paths_are_distinguished(tmp_path, monkeypatch):
    assert measure([tmp_path / 'absent'])['inodes'] == 0
    scan = os.scandir
    def denied(path):
        if Path(path) == tmp_path:
            raise PermissionError('test denial')
        return scan(path)
    monkeypatch.setattr(os, 'scandir', denied)
    usage = measure([tmp_path])
    with pytest.raises(ValueError, match='incomplete'):
        check_budget(usage, {'max_bytes': 10**10, 'max_inodes': 100})


def test_project_inventory_is_shared_between_short_tasks(tmp_path, monkeypatch):
    import storage_meter
    scan = storage_meter.measure
    calls = []
    def counted(roots):
        calls.append(roots)
        return scan(roots)
    monkeypatch.setattr(storage_meter, 'measure', counted)
    first = project_usage(tmp_path, 60)
    assert project_usage(tmp_path, 60) == first
    assert len(calls) == 1
    cached = tmp_path / 'work/storage_meter/latest.json'
    saved = json.loads(cached.read_text())
    saved['measured_at_seconds'] -= 61
    cached.write_text(json.dumps(saved))
    project_usage(tmp_path, 60)
    assert len(calls) == 2


@pytest.mark.parametrize('key,limit', [('allocated_bytes', 'max_bytes'), ('inodes', 'max_inodes')])
def test_start_reserve_blocks_before_limit(key, limit):
    settings = {'max_bytes': 10000, 'max_inodes': 100}
    usage = {'allocated_bytes': 0, 'inodes': 0, 'errors': []}
    usage[key] = int(settings[limit] * .85)
    with pytest.raises(ValueError, match='budget'):
        check_budget(usage, settings, .85)
    check_budget(usage, settings)


def settings():
    return {'max_bytes': 10**10, 'max_inodes': 100000, 'start_fraction': .85,
            'interval_seconds': .05, 'scratch_max_bytes': 10**9, 'scratch_min_free_bytes': 1}


def test_watchdog_records_success_and_failure(tmp_path):
    shared, local = tmp_path / 'shared', tmp_path / 'local'
    shared.mkdir(); local.mkdir()
    receipt = shared / 'metrics/job.json'
    command = [sys.executable, '-c', 'import pathlib,sys; pathlib.Path(sys.argv[1]).write_bytes(b"x"*10000)', str(local / 'result')]
    result = monitored(command, cwd=tmp_path, env=dict(os.environ), shared_root=shared,
                       scratch=local, sample_roots=[shared], receipt=receipt, settings=settings())
    assert result['state'] == 'complete' and result['observations'] >= 2
    assert result['scratch_peak']['logical_bytes'] == 10000
    assert receipt.with_suffix('.jsonl.gz').is_file()
    with pytest.raises(subprocess.CalledProcessError):
        monitored([sys.executable, '-c', 'raise SystemExit(7)'], cwd=tmp_path, env=dict(os.environ),
                  shared_root=shared, scratch=local, sample_roots=[shared], receipt=receipt, settings=settings())
    assert json.loads(receipt.read_text())['state'] == 'failed'


def test_watchdog_stops_child_when_shared_budget_exceeded(tmp_path):
    shared, local = tmp_path / 'shared', tmp_path / 'local'
    shared.mkdir(); local.mkdir()
    receipt = shared / 'metrics/job.json'
    cfg = dict(settings(), max_bytes=200000)
    command = [sys.executable, '-c',
               'import pathlib,sys,time; pathlib.Path(sys.argv[1]).write_bytes(b"x"*1000000); time.sleep(20)',
               str(shared / 'large')]
    with pytest.raises(ValueError, match='budget'):
        monitored(command, cwd=tmp_path, env=dict(os.environ), shared_root=shared,
                  scratch=local, sample_roots=[shared], receipt=receipt, settings=cfg)
    result = json.loads(receipt.read_text())
    assert result['state'] == 'failed' and result['shared_peak']['logical_bytes'] >= 1000000
