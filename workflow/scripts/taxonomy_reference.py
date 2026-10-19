"""Freeze one ETE taxonomy database for all native sample jobs."""
import csv
import fcntl
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

import busco_reference as shared
from common import now, write_json
from dataset_assets import ProductIntegrityError, UncachedVerification, digest, record, verify


DESTINATION = '/workspace/downloads/ete_taxonomy'
REQUIRED = ('taxa.sqlite', 'taxa.sqlite.traverse.pkl', 'taxdump.tar.gz')


def settings(root, config):
    if set(config) - {'cache_dir'}:
        raise ValueError('unknown build.taxonomy settings')
    cache = config.get('cache_dir', 'resources/taxonomy')
    if cache is not None:
        if not isinstance(cache, str) or not cache.strip() or '\n' in cache or '\r' in cache:
            raise ValueError('taxonomy.cache_dir must be null or a directory path without newlines')
        cache = str((Path(root) / cache).resolve())
    config['cache_dir'] = cache


def read_reference(root, expected=None):
    root = Path(root)
    marker = root / 'reference.json'
    if expected is not None:
        if Path(expected['receipt']['path']) != marker:
            raise ProductIntegrityError('taxonomy reference path differs from the frozen snapshot')
        verify(expected['receipt'])
    if not marker.is_file():
        raise ProductIntegrityError(f'incomplete taxonomy reference: {root}; use a new cache directory')
    data = json.loads(marker.read_text())
    inventory = data.get('files', [])
    if (data.get('schema_version') != 1 or data.get('kind') != 'ete-taxonomy'
            or not inventory or digest(inventory) != data.get('sha256')):
        raise ProductIntegrityError('taxonomy reference inventory changed')
    if expected is not None and data['sha256'] != expected['sha256']:
        raise ProductIntegrityError('taxonomy reference differs from the frozen snapshot')
    directory = root / 'dataset'
    files = shared.tree_files(directory, 'taxonomy')
    if [str(p.relative_to(directory)) for p in files] != [entry['path'] for entry in inventory]:
        raise ProductIntegrityError('taxonomy reference files are missing or unexpected')
    for entry in inventory:
        verify(dict(entry, path=str(directory / entry['path'])))
    if not set(REQUIRED + ('identity.json',)) <= {entry['path'] for entry in inventory}:
        raise ProductIntegrityError('incomplete taxonomy database')
    if any(not (directory / name).is_dir() for name in ('ete', 'ete4')):
        raise ProductIntegrityError('taxonomy cache directories are missing')
    if data.get('taxdump_sha256') != next(entry['sha256'] for entry in inventory if entry['path'] == 'taxdump.tar.gz'):
        raise ProductIntegrityError('taxonomy archive identity differs from its inventory')
    return {'sha256': data['sha256'], 'dataset': str(directory),
            'taxdump_sha256': data['taxdump_sha256'], 'receipt': record(marker)}


def mount_spec(source, destination, readonly=False):
    text = io.StringIO()
    csv.writer(text, lineterminator='').writerow(
        ['type=bind', f'src={source}', f'dst={destination}'] + (['ro'] if readonly else []))
    return text.getvalue()


def prepare(config, software, *, expected=None, cached_only=False):
    if not config.get('cache_dir'):
        return None
    root = Path(config['cache_dir']) / 'ete'
    if expected is not None or root.exists():
        return read_reference(root, expected)
    if cached_only:
        raise UncachedVerification('uncached taxonomy reference preparation')
    command = shared.runtime()
    root.parent.mkdir(parents=True, exist_ok=True)
    with (root.parent / 'ete.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if root.exists():
            return read_reference(root, expected)
        staging = Path(tempfile.mkdtemp(prefix='.ete-', dir=root.parent))
        try:
            env = dict(os.environ)
            for prefix in ('APPTAINER', 'SINGULARITY'):
                for key in ('BIND', 'BINDPATH', 'MOUNT'):
                    env.pop(prefix + '_' + key, None)
                env[prefix + '_TMPDIR'] = str(staging / 'tmp')
                env[prefix + 'ENV_PYTHONPYCACHEPREFIX'] = '/taxonomy/tmp/pycache'
            (staging / 'tmp').mkdir()
            # Use the pinned native builder, including its ETE version and NCBI URL.
            script = ('if [ -f /opt/conda/etc/profile.d/conda.sh ]; then '
                      'source /opt/conda/etc/profile.d/conda.sh; conda activate base; fi; '
                      'source /genegalleon/workflow/support/gg_util.sh; '
                      'ensure_ete_taxonomy_db /taxonomy/work')
            print(f'Taxonomy reference: preparing shared ETE database under {root}', flush=True)
            subprocess.run([command, 'exec', '--contain', '--mount', mount_spec(staging, '/taxonomy'),
                            '--mount', mount_spec(software['repository'], '/genegalleon', True),
                            software['image'], 'bash', '-euc', script], env=env, check=True)
            dataset = staging / 'dataset'
            (staging / 'work/downloads/ete_taxonomy').rename(dataset)
            if any(not (dataset / name).is_file() or not (dataset / name).stat().st_size for name in REQUIRED):
                raise ValueError('incomplete ETE taxonomy database')
            taxdump = record(dataset / 'taxdump.tar.gz')['sha256']
            identity = {name: record(dataset / name)['sha256'] for name in REQUIRED}
            write_json(dataset / 'identity.json', identity, compact=True)
            inventory = [dict(path=str(p.relative_to(dataset)), bytes=p.stat().st_size, sha256=record(p)['sha256'])
                         for p in shared.tree_files(dataset, 'taxonomy')]
            write_json(staging / 'reference.json', {'schema_version': 1, 'kind': 'ete-taxonomy',
                       'created_at': now(), 'taxdump_sha256': taxdump,
                       'sha256': digest(inventory), 'files': inventory}, compact=True)
            shutil.rmtree(staging / 'work')
            shutil.rmtree(staging / 'tmp')
            staging.rename(root)
            return read_reference(root)
        finally:
            if staging.exists():
                shutil.rmtree(staging)


def mount_environment(env, reference):
    shared.mount_environment(env, reference, DESTINATION, 'taxonomy')


def verify_mount(reference, software, env, work):
    data = json.loads(verify(reference['receipt']).read_text())
    identity = next(entry for entry in data['files'] if entry['path'] == 'identity.json')
    script = ('import hashlib,pathlib,sys; from ete4 import NCBITaxa; '
              'from ete4.ncbi_taxonomy.ncbiquery import is_taxadb_up_to_date; '
              'p=pathlib.Path(sys.argv[1]); '
              'sys.exit("shared taxonomy mount differs") if '
              'hashlib.sha256((p/"identity.json").read_bytes()).hexdigest()!=sys.argv[2] else None; '
              'sys.exit("incompatible taxonomy database") if not is_taxadb_up_to_date(str(p/"taxa.sqlite")) else None; '
              'NCBITaxa(dbfile=str(p/"taxa.sqlite"), update=False).get_rank([1])')
    shell = ('if [ -f /opt/conda/etc/profile.d/conda.sh ]; then '
             'source /opt/conda/etc/profile.d/conda.sh; conda activate base; fi; exec python -c "$@"')
    env = dict(env)
    for prefix in ('APPTAINER', 'SINGULARITY'):
        for key in ('BIND', 'BINDPATH'):
            env.pop(prefix + '_' + key, None)
    subprocess.run([shared.runtime(), 'exec', '--mount', mount_spec(work, '/workspace'), software['image'],
                    'bash', '-euc', shell, 'taxonomy-probe', script, DESTINATION, identity['sha256']],
                   env=env, check=True)
