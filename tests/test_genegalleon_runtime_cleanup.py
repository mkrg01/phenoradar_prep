"""Auxiliary caches are disposable, including after an interrupted launcher."""
import fcntl
import os
import select
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from dataset import cleanup_genegalleon, workspace
from dataset_assets import locked
from genegalleon_run import runtime_directory


@pytest.fixture
def runtime_workspace(tmp_path):
    manifest = {'analysis': {'output_root': 'output'},
                'config': {'storage': {'keep_intermediates': False}}}
    item = {'species': 'Plant_R1', 'row': {'species_id': 'Plant'}}
    work = workspace(tmp_path, manifest, item)
    runtime = runtime_directory(work)
    for child in ('tmp', 'pycache'):
        directory = runtime / child
        directory.mkdir(parents=True)
        (directory / 'helper').write_bytes(b'cache')
    out = work / 'output/transcriptome_assembly'
    scratch = out / 'tmp/1_Plant'
    scratch.mkdir(parents=True)
    (scratch / 'large').write_bytes(b'scratch')
    reads = out / 'amalgkit_getfastq/Plant/R1/reads.fastq.gz'
    reads.parent.mkdir(parents=True)
    reads.write_bytes(b'reads for retry')
    checkpoint = out / 'longest_cds/Plant_longestCDS.fa.gz'
    checkpoint.parent.mkdir()
    checkpoint.write_bytes(b'published checkpoint')
    (out / 'assembly.log').write_text('diagnostics')
    return tmp_path, manifest, item, work, runtime, scratch, reads, checkpoint


def test_auxiliary_cleanup_preserves_checkpoints_reads_and_foreign_namespaces(runtime_workspace):
    build, manifest, item, work, runtime, scratch, reads, checkpoint = runtime_workspace
    foreign = work / f'.genegalleon-runtime-{os.getuid() + 1}'
    foreign.mkdir()
    (foreign / 'keep').write_text('other UID namespace')
    external = build / 'external'
    external.mkdir()
    (external / 'keep').write_text('external file')
    (runtime / 'tmp/external').symlink_to(external, target_is_directory=True)
    preview = cleanup_genegalleon(build, manifest, item, 'sample', {}, apply=False)
    assert str(runtime.relative_to(build)) in preview['targets']
    assert runtime.exists() and scratch.exists()
    result = cleanup_genegalleon(build, manifest, item, 'sample', {})
    assert result['state'] == 'complete'
    assert not runtime.exists() and not scratch.exists()
    assert (foreign / 'keep').read_text() == 'other UID namespace'
    assert (external / 'keep').read_text() == 'external file'
    assert reads.read_bytes() == b'reads for retry'
    assert checkpoint.read_bytes() == b'published checkpoint'
    assert (checkpoint.parent.parent / 'assembly.log').read_text() == 'diagnostics'


def test_linked_runtime_root_is_retained(runtime_workspace):
    build, manifest, item, _, runtime, scratch, _, _ = runtime_workspace
    shutil.rmtree(runtime)
    external = build / 'external'
    external.mkdir()
    (external / 'keep').write_text('external file')
    runtime.symlink_to(external, target_is_directory=True)
    with pytest.warns(RuntimeWarning, match='cleanup incomplete'):
        result = cleanup_genegalleon(build, manifest, item, 'sample', {})
    assert result['state'] == 'pending' and scratch.exists()
    assert runtime.is_symlink() and (external / 'keep').exists()


def test_runtime_owned_by_another_uid_is_retained(runtime_workspace, monkeypatch):
    build, manifest, item, _, runtime, scratch, _, _ = runtime_workspace
    stat = Path.stat
    def foreign_owner(path, **kwargs):
        value = stat(path, **kwargs)
        if path == runtime:
            fields = list(value)
            fields[4] = os.getuid() + 1
            return os.stat_result(fields)
        return value
    monkeypatch.setattr(Path, 'stat', foreign_owner)
    with pytest.warns(RuntimeWarning, match='cleanup incomplete'):
        result = cleanup_genegalleon(build, manifest, item, 'sample', {})
    assert result['state'] == 'pending' and runtime.exists() and scratch.exists()


def test_native_lock_survives_closed_parent_descriptor(runtime_workspace):
    build, manifest, item, work, runtime, scratch, _, _ = runtime_workspace
    shutil.rmtree(scratch.parent)  # Bootstrap/final reporting has no assembly lock.
    with locked(work / '.native.lock') as handle:
        child = subprocess.Popen([sys.executable, '-B', '-c',
                                  'import sys; print("ready", flush=True); sys.stdin.read()'],
                                 pass_fds=(handle.fileno(),), stdin=subprocess.PIPE,
                                 stdout=subprocess.PIPE, text=True)
    try:
        assert select.select([child.stdout], [], [], 10)[0]
        assert child.stdout.readline().strip() == 'ready'
        # The launching parent is gone, but its child still holds the native lock.
        with pytest.warns(RuntimeWarning, match='cleanup incomplete'):
            result = cleanup_genegalleon(build, manifest, item, 'sample', {})
        assert result['state'] == 'pending' and runtime.exists()
        child.communicate(timeout=10)
        assert child.returncode == 0
    finally:
        if child.poll() is None:
            child.kill()
        child.communicate(timeout=10)
    result = cleanup_genegalleon(build, manifest, item, 'sample', {})
    assert result['state'] == 'complete' and not runtime.exists()


def test_linked_native_lock_is_retained(runtime_workspace):
    build, manifest, item, work, runtime, scratch, _, _ = runtime_workspace
    external = build / 'external.lock'
    external.write_text('external lock')
    (work / '.native.lock').symlink_to(external)
    with pytest.warns(RuntimeWarning, match='cleanup incomplete'):
        result = cleanup_genegalleon(build, manifest, item, 'sample', {})
    assert result['state'] == 'pending' and runtime.exists() and scratch.exists()
    assert external.read_text() == 'external lock'
