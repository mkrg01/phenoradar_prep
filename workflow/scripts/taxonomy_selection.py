"""Freeze taxonomy eligibility before importing products or scheduling samples."""
import json
from pathlib import Path
import sqlite3

from common import now, write_json, write_tsv
from dataset_assets import digest, record, verify, stat_identity

FIELDS = ['species', 'run', 'taxid', 'resolved_taxid', 'action', 'reason']


def enabled(config):
    # Old frozen builds have no policy and retain their original sample set.
    return config.get('taxonomy', {}).get('missing_taxid') is not None


def cached_database(config):
    settings = config.get('taxonomy', {})
    if settings.get('database'):
        path = Path(settings['database'])
        return path if path.is_file() else None
    cache = settings.get('cache_dir')
    if cache:
        for path in [Path(cache) / 'ete/dataset/taxa.sqlite', Path(cache) / 'taxa.sqlite']:
            if path.is_file():
                return path
    return None


def classify(items, database):
    """Resolve unique IDs offline; invalid databases are errors, not exclusions."""
    database = Path(database).resolve()
    before = stat_identity(database)
    resolved = {}
    with sqlite3.connect(database.as_uri() + '?mode=ro', uri=True) as db:
        db.execute('BEGIN')
        if db.execute('SELECT taxid FROM species WHERE taxid=1').fetchone() != (1,):
            raise ValueError('taxonomy snapshot has no root node')
        db.execute('SELECT taxid_old,taxid_new FROM merged LIMIT 0')
        for value in sorted({i['row']['taxid'] for i in items}, key=int):
            taxid = int(value)
            seen = set()
            while True:
                if taxid in seen:
                    raise ValueError(f'cycle in merged taxonomy IDs: {value}')
                seen.add(taxid)
                if db.execute('SELECT taxid FROM species WHERE taxid=?', (taxid,)).fetchone():
                    resolved[value] = str(taxid)
                    break
                target = db.execute('SELECT taxid_new FROM merged WHERE taxid_old=?', (taxid,)).fetchone()
                if target is None:
                    resolved[value] = ''
                    break
                taxid = int(target[0])
    if stat_identity(database) != before:
        raise ValueError('taxonomy snapshot changed during selection')
    rows = []
    for item in items:
        value = item['row']['taxid']
        target = resolved[value]
        rows.append(dict(species=item['species'], run=item['row']['run'], taxid=value,
                         resolved_taxid=target, action='include' if target else 'exclude',
                         reason=('taxid_merged' if target != value else '') if target else 'taxid_unresolved'))
    return rows


def partition(items, rows, policy):
    absent = [r for r in rows if r['action'] == 'exclude']
    if policy == 'error' and absent:
        raise ValueError('taxids absent from frozen taxonomy: ' + ', '.join(sorted({r['taxid'] for r in absent})))
    included = {r['species'] for r in rows if r['action'] == 'include'}
    return [i for i in items if i['species'] in included], [dict(r, source='taxonomy') for r in absent]


def preview(config, items):
    if not enabled(config):
        return items, [], False
    database = cached_database(config)
    if database is None:
        return items, [], True
    selected, excluded = partition(items, classify(items, database), config['taxonomy']['missing_taxid'])
    return selected, excluded, False


def apply(path, manifest):
    filename = Path(path) / 'taxonomy_selection.json'
    if not enabled(manifest['config']) or not filename.is_file():
        return
    data = json.loads(filename.read_text())
    if data.get('schema_version') != 1 or digest({k: v for k, v in data.items() if k != 'sha256'}) != data.get('sha256'):
        raise ValueError('taxonomy selection record changed')
    prepared = json.loads((Path(path) / 'checksums.json').read_text())['build.json']
    if data['prepared_sha256'] != prepared or data['policy'] != manifest['config']['taxonomy']['missing_taxid']:
        raise ValueError('taxonomy selection belongs to a different build or policy')
    verify(data['database'])
    link = Path(path) / 'taxonomy.sqlite'
    if not link.is_file() or link.resolve() != Path(data['database']['path']).resolve():
        raise ValueError('build taxonomy link differs from selected snapshot')
    for entry in data['reports']:
        verify(entry)
    partition([], data['rows'], data['policy'])
    names = [i['species'] for i in manifest['items']]
    if not manifest.get('_taxonomy_applied') and names != data['candidate_species']:
        raise ValueError('taxonomy candidates differ from frozen metadata')
    selected = set(data['included_species'])
    if not selected <= set(data['candidate_species']):
        raise ValueError('taxonomy selection contains unknown samples')
    manifest['items'] = [i for i in manifest['items'] if i['species'] in selected]
    exclusions = {r['species']: r for r in manifest.get('excluded', [])}
    exclusions.update({r['species']: dict(r, source='taxonomy') for r in data['rows'] if r['action'] == 'exclude'})
    manifest['excluded'] = list(exclusions.values())
    manifest['_taxonomy_applied'] = True


def ensure(path, manifest, *, software=None, cached_only=False):
    if not enabled(manifest['config']):
        return
    path = Path(path)
    filename = path / 'taxonomy_selection.json'
    if not filename.exists():
        if cached_only:
            from dataset_assets import UncachedVerification
            raise UncachedVerification('taxonomy eligibility has not been frozen')
        database = cached_database(manifest['config'])
        if database is None:
            cfg = manifest['config']
            if cfg['taxonomy'].get('database'):
                raise ValueError(f"taxonomy database missing: {cfg['taxonomy']['database']}")
            from dataset_software import resolve
            from taxonomy_reference import prepare
            if software is None:
                software = resolve(cfg['genegalleon'])
            cfg['genegalleon'], manifest['genegalleon'], manifest['software_lock'] = software
            settings = dict(cfg['taxonomy'])
            settings['cache_dir'] = settings.get('cache_dir') or str(Path(manifest['root']) / 'resources/taxonomy')
            reference = prepare(settings, cfg['genegalleon'], expected=manifest.get('taxonomy_reference'))
            manifest['taxonomy_reference'] = reference
            database = Path(reference['dataset']) / 'taxa.sqlite'
        frozen = record(database)
        rows = classify(manifest['items'], database)
        verify(frozen)
        link = path / 'taxonomy.sqlite'
        if link.exists() or link.is_symlink():
            if link.resolve() != Path(frozen['path']).resolve():
                raise ValueError('existing taxonomy link points to another snapshot')
        else:
            link.symlink_to(Path(frozen['path']).resolve())
        reports = []
        for name, selected in [('taxonomy_selection.tsv', rows),
                               ('taxonomy_excluded_runs.tsv', [r for r in rows if r['action'] == 'exclude'])]:
            write_tsv(path / name, FIELDS, selected)
            reports.append(record(path / name))
        data = {'schema_version': 1, 'created_at': now(),
                'prepared_sha256': json.loads((path / 'checksums.json').read_text())['build.json'],
                'policy': manifest['config']['taxonomy']['missing_taxid'], 'database': frozen,
                'candidate_species': [i['species'] for i in manifest['items']],
                'included_species': [r['species'] for r in rows if r['action'] == 'include'],
                'rows': rows, 'reports': reports}
        data['sha256'] = digest(data)
        write_json(filename, data)
    apply(path, manifest)
    data = json.loads(filename.read_text())
    # Keep failed strict-policy/all-excluded attempts inspectable and resumable.
    partition([], data['rows'], data['policy'])
    if not manifest['items']:
        raise ValueError('all candidate samples lack resolvable taxonomy; no jobs were submitted')
    from build_reports import read_context, save_context
    from dataset import excluded_report
    context = read_context(path)
    planned = {r['species']: r for r in context['items']}
    report = [planned[i['species']] for i in manifest['items']] + excluded_report(manifest['excluded'])
    save_context(path, manifest, report, deferred=context.get('deferred', False))


def require(path, manifest):
    if enabled(manifest.get('config', {})) and not manifest.get('_taxonomy_applied'):
        raise ValueError('taxonomy selection pending; run build verification before scheduling samples or mapping')
    if manifest.get('_taxonomy_applied') and not manifest['items']:
        raise ValueError('all candidate samples lack resolvable taxonomy; no jobs were submitted')
