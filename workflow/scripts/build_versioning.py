"""Timestamped build names and the latest successfully completed snapshot."""
import fcntl
import os
import shutil
import tempfile
from datetime import datetime, timezone
from pathlib import Path

from dataset_assets import SAFE


def naming_mode(config):
    mode = config.get('name_mode', 'fixed')
    if mode not in ('fixed', 'timestamp'):
        raise ValueError('name_mode must be fixed or timestamp')
    return mode


def timestamp_directory(root, prefix, prepared_at):
    """Choose a free name while the caller holds the dataset preparation lock."""
    base = datetime.fromisoformat(prepared_at).astimezone(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    directory = Path(root) / 'results' / f'{prefix}_{base}'
    number = 1
    while os.path.lexists(directory):
        number += 1
        directory = directory.with_name(f'{prefix}_{base}_{number:02d}')
    return directory


def metadata_history(provenance):
    """Expose known timestamps without inventing a fetch date for external tables."""
    stage = provenance.get('source_stage', {})
    while isinstance(stage, dict) and stage.get('kind') != 'fetch':
        origin = stage.get('source_origin')
        stage = origin.get('stage') if isinstance(origin, dict) else None
    return {'selected_at': provenance.get('created_at'),
            'accepted_at': provenance.get('accepted_at'),
            'fetched_at': stage.get('created_at') if isinstance(stage, dict) else None}


def publish_latest(build, completed):
    """Atomically advance a dataset alias only to a newer completed build."""
    if completed.get('name_mode') != 'timestamp':
        return None
    from build_products import load_complete
    build = Path(build).resolve()
    dataset = completed['dataset_name']
    if not isinstance(dataset, str) or not SAFE.fullmatch(dataset):
        raise ValueError('dataset name must be a simple directory name')
    latest = build.parent / f'{dataset}_latest'
    if latest == build:
        raise ValueError('build name is reserved for the latest dataset alias')

    def order(data):
        return datetime.fromisoformat(data['prepared_at']), data['build_id']

    with (build.parent / f'.{dataset}.latest.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if os.path.lexists(latest):
            if not latest.is_symlink():
                raise ValueError(f'latest path is not a managed build symlink: {latest}')
            previous = latest.resolve()
            if previous.parent != build.parent:
                raise ValueError(f'latest alias points outside the build directory: {latest}')
            old = load_complete(previous, verify_files=False)
            if old.get('dataset_name') != dataset or old.get('name_mode') != 'timestamp':
                raise ValueError(f'latest alias belongs to a different dataset: {latest}')
            if order(old) >= order(completed):
                return latest
        staging = Path(tempfile.mkdtemp(prefix=f'.{dataset}.latest-', dir=build.parent))
        try:
            link = staging / 'link'
            link.symlink_to(build.name, target_is_directory=True)
            os.replace(link, latest)
        finally:
            shutil.rmtree(staging)
    return latest
