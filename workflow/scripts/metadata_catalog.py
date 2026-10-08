#!/usr/bin/env python3
"""Curated metadata candidates, independent of database job submission."""
import argparse
import collections
import csv
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile

import yaml

from accession_exclusions import read_exclusions
from build_versioning import metadata_history
from common import file_record, now, read_tsv, sha256, write_json, write_tsv
from dataset_assets import SAFE, counts, locked, verify
from sample_identity import annotate

CONFIG_KEYS = {'schema_version', 'name', 'build_config', 'search_string', 'sample_group',
               'rule_set', 'rules', 'accepted_samples', 'previous_metadata',
               'excluded_accessions'}
ACCEPTED_FIELDS = ['taxid', 'run', 'bioproject', 'scientific_name',
                   'busco_complete', 'busco_total', 'source_build']
EXCLUSION_FIELDS = ['accession', 'taxid', 'bioproject', 'reason', 'source_build']
REPORT_FIELDS = ['taxid', 'scientific_name', 'previous_run', 'selected_run', 'bioproject',
                 'previous_bioproject', 'previous_total_bases', 'selected_total_bases',
                 'source_taxid', 'source_bioproject', 'source_scientific_name',
                 'decision', 'review', 'reason', 'eligible_candidates']
csv.field_size_limit(32 * 1024 * 1024)


def path_at(root, value):
    return (Path(root) / value).resolve()


def configuration(root, path):
    path = path_at(root, path)
    cfg = yaml.safe_load(path.read_text())
    if (not isinstance(cfg, dict) or not CONFIG_KEYS <= set(cfg)
            or set(cfg) - CONFIG_KEYS - {'overrides', 'busco_threshold'} or cfg['schema_version'] != 1):
        raise ValueError('selection config requires the documented schema_version: 1 keys')
    for key in CONFIG_KEYS - {'schema_version'}:
        if not isinstance(cfg[key], str) or not cfg[key].strip():
            raise ValueError(f'selection.{key} must be a nonempty string')
    if not SAFE.fullmatch(cfg['name']) or not SAFE.fullmatch(cfg['rule_set']):
        raise ValueError('dataset name and rule_set must be simple identifiers')
    cfg.setdefault('overrides', None)
    if cfg['overrides'] is not None and (not isinstance(cfg['overrides'], str) or not cfg['overrides'].strip()):
        raise ValueError('overrides must be null or a TSV path')
    # New samples below this threshold become run exclusions during QC recording.
    # Existing adoption records remain authoritative during later selection.
    cfg.setdefault('busco_threshold', 0.5)
    threshold = cfg['busco_threshold']
    if type(threshold) not in (int, float) or not 0 <= threshold <= 1:
        raise ValueError('busco_threshold must be between zero and one')
    cfg['_path'] = path
    return cfg


def table(path, required):
    with Path(path).open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fields = reader.fieldnames
        if not fields or len(fields) != len(set(fields)) or not set(required) <= set(fields):
            raise ValueError(f'{path}: required unique columns: {", ".join(required)}')
        rows = list(reader)
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise ValueError(f'malformed TSV: {path}')
    return fields, rows


def taxid(value):
    if not str(value).isdigit() or int(value) <= 0:
        raise ValueError(f'positive integer taxid required: {value!r}')
    return str(int(value))


def unique_by(rows, key):
    result = {}
    for row in rows:
        value = taxid(row[key]) if key == 'taxid' else row[key]
        if value in result:
            raise ValueError(f'duplicate {key}: {value}')
        result[value] = row
    return result


def resolved_exclusions(source, historical, source_build):
    decisions = []
    for run, decision in historical.items():
        row = {**source.get(run, {}), **{k: v for k, v in decision.items() if v}}
        decisions.append(dict(accession=run, taxid=taxid(row['taxid']) if row.get('taxid') else '',
                              bioproject=row.get('bioproject', ''),
                              reason=decision.get('reason', ''),
                              source_build=decision.get('source_build') or source_build))
    return sorted(decisions, key=lambda r: r['accession'])


def busco_exclusion(row, threshold, source_build):
    return dict(accession=row['run'], taxid=row['taxid'], bioproject=row['bioproject'],
                reason=f'busco_completeness_below_{threshold}', source_build=source_build)


def append_exclusions(path, decisions):
    """Keep existing decisions and custom columns while adding new run exclusions."""
    if not decisions:
        return
    fields, rows = table(path, ['accession'])
    existing = {row['accession'] for row in rows}
    additions = []
    for decision in decisions:
        if decision['accession'] not in existing:
            additions.append(decision)
            existing.add(decision['accession'])
    if additions:
        fields = list(dict.fromkeys([*fields, *EXCLUSION_FIELDS]))
        write_tsv(path, fields, [{key: row.get(key, '') for key in fields}
                                 for row in [*rows, *additions]])


def initialize_empty(root, cfg, exclusions=None, replace=False):
    """Prepare a fresh-query dataset, retaining only enriched exclusion decisions."""
    source = path_at(root, exclusions) if exclusions else None
    original = file_record(source) if source else None
    decisions = []
    if source:
        historical = read_exclusions(source)
        decisions = resolved_exclusions({}, historical, 'historical')
    outputs = [path_at(root, cfg[k]) for k in ('previous_metadata', 'accepted_samples', 'excluded_accessions')]
    with locked(outputs[0].parent / '.metadata.lock'):
        if not replace and any(p.exists() for p in outputs):
            raise ValueError('dataset already exists; --replace deliberately clears representatives and success evidence')
        write_tsv(outputs[0], ['scientific_name', 'taxid', 'run', 'bioproject'], [])
        write_tsv(outputs[1], ACCEPTED_FIELDS, [])
        write_tsv(outputs[2], EXCLUSION_FIELDS, decisions)
        write_json(outputs[0].parent / 'provenance.json', {
            'schema_version': 1, 'kind': 'awaiting_fresh_query', 'created_at': now(),
            'selection_config': file_record(cfg['_path']), 'original_exclusions': original,
            'metadata': file_record(outputs[0]), 'accepted_samples': file_record(outputs[1]),
            'excluded_accessions': file_record(outputs[2]), 'samples': 0,
            'note': 'Start with amalgkit metadata; no historical representatives or products are reused.'})
    return {'samples': 0, 'accepted_samples': 0, 'exclusions': len(decisions)}


def seed_history(root, cfg, metadata, exclusions, source_build='historical', replace=False):
    """Import comparison metadata and exclusions without claiming container/QC success."""
    source = path_at(root, metadata)
    fields, rows = table(source, ['scientific_name', 'run', 'taxid', 'bioproject'])
    unique_by(rows, 'taxid')
    by_run = unique_by(rows, 'run')
    decisions = resolved_exclusions(by_run, read_exclusions(path_at(root, exclusions)), source_build)
    outputs = [path_at(root, cfg[k]) for k in ('previous_metadata', 'accepted_samples', 'excluded_accessions')]
    with locked(outputs[0].parent / '.metadata.lock'):
        if not replace and any(p.exists() for p in outputs):
            raise ValueError('metadata history already exists; --replace deliberately resets success evidence')
        write_tsv(outputs[0], fields, sorted(rows, key=lambda r: int(r['taxid'])))
        write_tsv(outputs[1], ACCEPTED_FIELDS, [])
        write_tsv(outputs[2], EXCLUSION_FIELDS, decisions)
        write_json(outputs[0].parent / 'provenance.json', {
            'schema_version': 1, 'kind': 'historical_metadata_seed', 'created_at': now(),
            'source_metadata': file_record(source), 'original_exclusions': file_record(path_at(root, exclusions)),
            'selection_config': file_record(cfg['_path']), 'metadata': file_record(outputs[0]),
            'accepted_samples': file_record(outputs[1]), 'excluded_accessions': file_record(outputs[2]),
            'samples': len(rows), 'note': 'Historical results are not retained as container-verified successes.'})
    return {'historical_samples': len(rows), 'accepted_samples': 0, 'exclusions': len(decisions)}


def normalized_row(row):
    """Port the two historical notebooks' name cleanup; preserve submitted names."""
    row = dict(row)
    original = row['scientific_name']
    name = re.sub(r'[\[\]().,]', '', original)
    name = re.sub(r'[:/#]', '-', name).strip()
    name = re.sub(r'^([^\s]+)$', r'\1_XXX', name)
    name = re.sub(r'^(\S+)\s+', r'\1_', name)
    name = re.sub(r'\s+', '-', name).replace('--', '-').replace('_', ' ')
    row['scientific_name_original'] = row.get('scientific_name_original') or original
    row['scientific_name'] = name
    row['taxid'] = taxid(row['taxid'])
    # Input identities are recomputed only for newly selected rows.
    row.pop('species_id', None)
    row.pop('analysis_sample_id', None)
    return annotate(row)


def previous_state(root, cfg):
    """Read optional metadata and validate the authoritative adoption registry."""
    fields, rows = history_table(path_at(root, cfg['previous_metadata']),
                                ['scientific_name', 'taxid', 'run', 'bioproject'])
    previous = unique_by(rows, 'taxid')
    unique_by(rows, 'run')
    _, accepted_rows = history_table(path_at(root, cfg['accepted_samples']), ACCEPTED_FIELDS)
    accepted = unique_by(accepted_rows, 'taxid')
    unique_by(accepted_rows, 'run')
    for tid, row in accepted.items():
        if not SAFE.fullmatch(row['run']) or not row['bioproject'].strip():
            raise ValueError(f'invalid accepted sample identity: taxid {tid}')
        annotate(row)
        complete, total = int(row['busco_complete']), int(row['busco_total'])
        if total <= 0 or not 0 <= complete <= total:
            raise ValueError(f'invalid accepted BUSCO counts: taxid {tid}')
    return fields, previous, accepted


def retained_row(source, accepted):
    """Refresh run attributes while preserving the reviewed sample identity."""
    row = dict(source)
    row['scientific_name_original'] = source.get('scientific_name_original') or source['scientific_name']
    row.update({key: accepted[key] for key in ACCEPTED_FIELDS[:4]})
    row['taxid'] = taxid(row['taxid'])
    row.pop('species_id', None)
    row.pop('analysis_sample_id', None)
    return annotate(row)


def retained_source_changes(source, accepted, cfg):
    changes = []
    try:
        if taxid(source['taxid']) != taxid(accepted['taxid']):
            changes.append('taxid_changed')
    except ValueError:
        changes.append('invalid_taxid')
    if source['bioproject'] != accepted['bioproject']:
        changes.append('bioproject_changed')
    try:
        source_name = normalized_row(source)['scientific_name']
    except ValueError:
        source_name = None
    if source_name != accepted['scientific_name']:
        changes.append('scientific_name_changed')
    if source['sample_group'] != cfg['sample_group']:
        changes.append('sample_group')
    if source['exclusion'] != 'no':
        changes.append('amalgkit_exclusion')
    if not source['bioproject'].strip():
        changes.append('missing_bioproject')
    try:
        bases = Decimal(source['total_bases'])
        if not bases.is_finite() or bases <= 0:
            raise InvalidOperation
    except InvalidOperation:
        changes.append('invalid_total_bases')
    return changes


def history_table(path, fields):
    """A dataset has no representative or success history before its first run."""
    if not path.exists():
        return fields, []
    return table(path, fields)


def selection_inputs(root, cfg):
    """Record absent first-run history as well as existing input checksums."""
    inputs = {}
    for key in ('build_config', 'rules', 'previous_metadata', 'accepted_samples',
                'excluded_accessions', 'overrides'):
        if cfg[key] is None:
            inputs[key] = None
            continue
        path = path_at(root, cfg[key])
        if key in ('previous_metadata', 'accepted_samples', 'overrides') and not path.exists():
            inputs[key] = {'path': str(path), 'missing': True}
        else:
            inputs[key] = file_record(path)
    return inputs


def matching_container(root, cfg, build):
    requested = yaml.safe_load(path_at(root, cfg['build_config']).read_text())['genegalleon']
    used = build.get('config', {}).get('genegalleon', {})
    lock = build.get('software_lock') or {}
    images = lock.get('container', {}).get('files', [])
    if (any(used.get(k) != requested.get(k) for k in ('image_uri', 'version', 'revision'))
            or not images or any(not e.get('sha256') or not e.get('bytes') for e in images)):
        raise ValueError('build has no matching pinned-container receipt; use seed for historical metadata')
    return lock


def record_successes(root, cfg, build_path, runs=None, dry_run=False, *, exclude_failed=False, plan=None):
    """Record adoption, reviewed exclusions, and cleanup through one operation."""
    from sample_review import record
    return record(root, cfg, build_path, runs, dry_run, exclude_failed=exclude_failed, plan=plan)


def initialize(root, cfg, database, source_metadata, exclusions, replace=False):
    """Import completed-product evidence and enrich historical exclusion decisions."""
    from build_products import load_complete
    from portable_build import completion_path
    database = path_at(root, database)
    manifest = completion_path(database)
    data = load_complete(database, verify_files=False)
    # Check the portable metadata itself, without hashing all large sample products.
    inventory = {entry['path']: entry for entry in data['files']}
    metadata_path = Path(data['bundle_root']) / 'metadata.tsv'
    verify(inventory[str(metadata_path)])
    provenance_path = Path(data['bundle_root']) / 'provenance/build.json'
    if str(provenance_path) not in inventory:
        raise ValueError('completed database lacks pinned-container build provenance; use seed for historical metadata')
    verify(inventory[str(provenance_path)])
    build = json.loads(provenance_path.read_text())
    lock = matching_container(root, cfg, build)
    fields, metadata = table(metadata_path, ['scientific_name', 'taxid', 'run', 'bioproject'])
    unique_by(metadata, 'taxid')
    products = {p['row']['run']: p for p in data['products'].values()}
    _, source = table(path_at(root, source_metadata), ['run', 'taxid', 'bioproject'])
    source = unique_by(source, 'run')
    original_exclusions = file_record(path_at(root, exclusions))
    historical = read_exclusions(path_at(root, exclusions))
    decisions = resolved_exclusions(source, historical, data['build_id'])
    accepted, selected, excluded_runs = [], [], []
    for row in metadata:
        product = products.get(row['run'])
        if product is None or any(row[k] != product['row'][k] for k in ('run', 'taxid', 'scientific_name')):
            raise ValueError(f'completed database metadata/product mismatch: {row["run"]}')
        qc = counts(product['counts'])
        complete = qc['busco_cds_single'] + qc['busco_cds_duplicated']
        total = qc['busco_cds_total']
        if row['run'] in historical:
            continue
        if complete / total < cfg['busco_threshold']:
            decisions.append(busco_exclusion(row, cfg['busco_threshold'], data['build_id']))
            excluded_runs.append(row['run'])
            continue
        accepted.append({**{k: row[k] for k in ACCEPTED_FIELDS[:4]},
                         'busco_complete': complete, 'busco_total': total,
                         'source_build': data['build_id']})
        selected.append(row)
    outputs = [path_at(root, cfg[k]) for k in ('previous_metadata', 'accepted_samples', 'excluded_accessions')]
    provenance = outputs[0].parent / 'provenance.json'
    with locked(outputs[0].parent / '.metadata.lock'):
        if not replace and any(p.exists() for p in [*outputs, provenance]):
            raise ValueError('baseline already exists; use --replace only to deliberately record a new completed baseline')
        write_tsv(outputs[0], fields, sorted(selected, key=lambda r: int(r['taxid'])))
        write_tsv(outputs[1], ACCEPTED_FIELDS, sorted(accepted, key=lambda r: int(r['taxid'])))
        write_tsv(outputs[2], EXCLUSION_FIELDS, sorted(decisions, key=lambda r: r['accession']))
        write_json(provenance, {'schema_version': 1, 'kind': 'imported_completed_baseline',
                   'created_at': now(), 'database_manifest': file_record(manifest),
                   'source_metadata': file_record(path_at(root, source_metadata)),
                   'original_exclusions': original_exclusions,
                   'selection_config': file_record(cfg['_path']),
                   'busco_threshold': cfg['busco_threshold'],
                   'excluded_runs': excluded_runs,
                   'metadata': file_record(outputs[0]), 'accepted_samples': file_record(outputs[1]),
                   'excluded_accessions': file_record(outputs[2]), 'samples': len(selected),
                   'software_lock': lock})
    return {'accepted_samples': len(accepted), 'exclusions': len(decisions),
            'excluded_runs': excluded_runs}


def select_candidates(root, cfg, source, output):
    """Keep QC-approved accessions, then fill unrepresented taxids deterministically."""
    source, output = path_at(root, source), path_at(root, output)
    if output.exists():
        raise ValueError(f'candidate directory already exists; choose a new attempt: {output}')
    if not output.is_relative_to(Path(root).resolve() / 'work'):
        raise ValueError('candidate outputs must be under the project work/ directory')
    inputs = selection_inputs(root, cfg)
    config_record = file_record(cfg['_path'])
    fields, previous, accepted = previous_state(root, cfg)
    exclusions = read_exclusions(path_at(root, cfg['excluded_accessions']))
    override_rows = []
    if cfg['overrides'] is not None:
        _, override_rows = history_table(path_at(root, cfg['overrides']), ['taxid', 'run', 'reason'])
    overrides = unique_by(override_rows, 'taxid')
    if any(not r['reason'].strip() or not SAFE.fullmatch(r['run']) for r in overrides.values()):
        raise ValueError('overrides require a valid run and an explicit reason')
    active_accepted = {tid: row for tid, row in accepted.items() if row['run'] not in exclusions}
    accepted_by_run = {row['run']: row for row in active_accepted.values()}
    previous_by_run = {row['run']: row for row in previous.values()}
    kept = {tid: retained_row(previous_by_run[row['run']], row) for tid, row in active_accepted.items()
            if row['run'] in previous_by_run}
    best, forced, observed, candidate_counts = {}, {}, {}, collections.Counter()
    source_changes = {}
    reasons, seen_runs, observed_taxids = collections.Counter(), set(), set()
    before = file_record(source)
    with source.open(newline='') as handle:
        reader = csv.reader(handle, delimiter='\t')
        source_fields = next(reader, [])
        required = ['scientific_name', 'run', 'taxid', 'bioproject', 'sample_group', 'exclusion', 'total_bases']
        if len(source_fields) != len(set(source_fields)) or not set(required) <= set(source_fields):
            raise ValueError('curated source requires unique columns: ' + ', '.join(required))
        index = {key: source_fields.index(key) for key in required}
        for number, values in enumerate(reader, 2):
            if len(values) != len(source_fields):
                raise ValueError(f'malformed source row: {number}')
            run = values[index['run']]
            if not run or not SAFE.fullmatch(run):
                reasons['invalid_run'] += 1
                continue
            if run in seen_runs:
                raise ValueError(f'duplicate source run: {run}')
            seen_runs.add(run)
            if run in accepted_by_run:
                evidence = accepted_by_run[run]
                tid = taxid(evidence['taxid'])
                incoming = dict(zip(source_fields, values))
                observed[tid] = incoming
                observed_taxids.add(tid)
                kept[tid] = retained_row(incoming, evidence)
                source_changes[tid] = retained_source_changes(incoming, evidence, cfg)
                rejection = next((reason for reason in source_changes[tid]
                                  if reason in {'invalid_taxid', 'sample_group', 'amalgkit_exclusion',
                                                'missing_bioproject', 'invalid_total_bases'}), '')
                if rejection:
                    reasons[rejection] += 1
                else:
                    candidate_counts[tid] += 1
                if tid in overrides and overrides[tid]['run'] == run:
                    forced[tid] = kept[tid]
                continue
            try:
                tid = taxid(values[index['taxid']])
            except ValueError:
                reasons['invalid_taxid'] += 1
                continue
            project = values[index['bioproject']]
            reason = ''
            if values[index['sample_group']] != cfg['sample_group']:
                reason = 'sample_group'
            elif values[index['exclusion']] != 'no':
                reason = 'amalgkit_exclusion'
            elif run in exclusions:
                reason = 'excluded_run'
            elif not project.strip():
                reason = 'missing_bioproject'
            if values[index['sample_group']] == cfg['sample_group']:
                observed_taxids.add(tid)
            if reason:
                reasons[reason] += 1
                continue
            try:
                bases = Decimal(values[index['total_bases']])
                if not bases.is_finite() or bases <= 0:
                    raise InvalidOperation
            except InvalidOperation:
                reasons['invalid_total_bases'] += 1
                continue
            candidate_counts[tid] += 1
            # Reserve reviewed taxids even when their accession is absent from this source.
            if tid in overrides and run == overrides[tid]['run']:
                forced[tid] = normalized_row(dict(zip(source_fields, values)))
            if tid in active_accepted or tid in overrides:
                continue
            rank = (-bases, run)
            if tid not in best or rank < best[tid][0]:
                best[tid] = (rank, normalized_row(dict(zip(source_fields, values))))
    if file_record(source) != before:
        raise ValueError('source metadata changed during selection')
    if file_record(cfg['_path']) != config_record or selection_inputs(root, cfg) != inputs:
        raise ValueError('selection inputs changed during selection')
    for tid, override in overrides.items():
        if tid not in forced:
            if tid in kept and override['run'] == kept[tid]['run']:
                forced[tid] = kept[tid]
            else:
                raise ValueError(f'override is absent or excluded from eligible source candidates: {tid}/{override["run"]}')
    selected = {tid: row for tid, (_, row) in best.items()}
    selected.update(kept)
    selected.update(forced)
    unresolved = sorted(row['run'] for tid, row in active_accepted.items()
                        if tid not in kept and tid not in forced)
    reports = []
    for tid in sorted(set(previous) | set(active_accepted) | observed_taxids | set(selected), key=int):
        old, new = active_accepted.get(tid) or previous.get(tid), selected.get(tid)
        cached = previous_by_run.get(old['run'], {}) if old else {}
        incoming = observed.get(tid, {})
        review, reason = False, ''
        if tid in overrides:
            decision, reason = 'override', overrides[tid]['reason']
        elif tid in active_accepted:
            decision = 'retained' if tid in kept else 'unresolved'
            reason = 'previous_completed_sample'
            if tid not in observed:
                review, reason = True, 'retained_run_absent_from_source'
            elif source_changes[tid]:
                review, reason = True, 'retained_run_source_status:' + ';'.join(source_changes[tid])
        elif new:
            decision = ('reselected' if old['run'] == new['run'] else 'replacement') if old else 'new'
            reason = 'largest_total_bases'
        else:
            decision, reason = 'no_candidate', 'no_eligible_run'
        reports.append({'taxid': tid, 'scientific_name': (new or old or {}).get('scientific_name', ''),
                        'previous_run': old['run'] if old else '', 'selected_run': new['run'] if new else '',
                        'bioproject': new['bioproject'] if new else '', 'decision': decision,
                        'previous_bioproject': old.get('bioproject', '') if old else '',
                        'previous_total_bases': cached.get('total_bases', ''),
                        'selected_total_bases': new.get('total_bases', '') if new else '',
                        'source_taxid': incoming.get('taxid', ''),
                        'source_bioproject': incoming.get('bioproject', ''),
                        'source_scientific_name': incoming.get('scientific_name', ''),
                        'review': 'yes' if review else 'no', 'reason': reason,
                        'eligible_candidates': candidate_counts[tid]})
    merged_fields = list(dict.fromkeys([*fields, *source_fields, 'scientific_name_original',
                                      'species_id', 'analysis_sample_id']))
    rows = [annotate(selected[tid]) for tid in sorted(selected, key=int)]
    unique_by(rows, 'run')
    # Match the workflow's normalized identity/collision checks before publishing a candidate.
    names = [r['analysis_sample_id'].replace('-', '_') for r in rows]
    if len(set(names)) != len(names):
        raise ValueError('selected metadata has colliding normalized sample identities')
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix='.metadata-candidate-', dir=output.parent))
    try:
        write_tsv(staging / 'metadata.tsv', merged_fields, [{k: r.get(k, '') for k in merged_fields} for r in rows])
        write_tsv(staging / 'selection.tsv', REPORT_FIELDS, reports)
        summary = {'samples': len(rows), 'decisions': dict(collections.Counter(r['decision'] for r in reports)),
                   'review_required': sum(r['review'] == 'yes' for r in reports), 'source_rejections': dict(reasons),
                   'unresolved_accepted_runs': unresolved}
        provenance = {'schema_version': 1, 'kind': 'metadata_candidate', 'created_at': now(),
                      'dataset': cfg['name'], 'source_metadata': before, 'selection_config': config_record,
                      'inputs': inputs, 'generator': file_record(Path(__file__)), 'summary': summary,
                      'metadata': file_record(staging / 'metadata.tsv'), 'selection': file_record(staging / 'selection.tsv')}
        stage_info = source.parent / 'provenance.json'
        if stage_info.exists():
            provenance['source_provenance'] = file_record(stage_info)
            provenance['source_stage'] = json.loads(stage_info.read_text())
        else:
            provenance['source_stage'] = {'kind': 'external_curated_metadata'}
            query_info = source.parent / 'metadata.query_info.json'
            if query_info.exists():
                provenance['source_query_info'] = file_record(query_info)
        for key in ('metadata', 'selection'):
            provenance[key]['path'] = str(output / Path(provenance[key]['path']).name)
        write_json(staging / 'provenance.json', provenance)
        staging.rename(output)
    finally:
        if staging.exists(): shutil.rmtree(staging)
    return summary


def accept_candidate(root, cfg, candidate, allow_review=False):
    candidate = path_at(root, candidate)
    data = json.loads((candidate / 'provenance.json').read_text())
    if data.get('kind') != 'metadata_candidate' or data.get('dataset') != cfg['name']:
        raise ValueError('candidate does not belong to this dataset')
    if sha256(cfg['_path']) != data['selection_config']['sha256']:
        raise ValueError('selection config changed; regenerate the candidate')
    current_inputs = selection_inputs(root, cfg)
    for key, entry in data['inputs'].items():
        if current_inputs[key] != entry:
            raise ValueError(f'{key} changed; regenerate the candidate')
    for key in ('metadata', 'selection'):
        if sha256(candidate / f'{key}.tsv') != data[key]['sha256']:
            raise ValueError(f'candidate {key} changed; use overrides and regenerate')
    unresolved = data.get('summary', {}).get('unresolved_accepted_runs', [])
    if unresolved:
        raise ValueError('accepted runs lack source metadata and a local fallback: ' + ', '.join(unresolved)
                         + '; restore their metadata or explicitly exclude/replace them, then regenerate')
    # Source changes remain visible in selection.tsv. An accepted representative
    # is retained until an explicit exclusion or override changes the decision.
    target = path_at(root, cfg['previous_metadata'])
    # Keep success evidence only for unchanged representatives; additions need actual build/QC results.
    _, rows = table(candidate / 'metadata.tsv', ['taxid', 'run'])
    selected = unique_by(rows, 'taxid')
    _, evidence = history_table(path_at(root, cfg['accepted_samples']), ACCEPTED_FIELDS)
    evidence = [r for r in evidence if r['taxid'] in selected and
                all(r[k] == selected[r['taxid']][k] for k in ('run', 'scientific_name', 'bioproject'))]
    with locked(target.parent / '.metadata.lock'):
        # Recheck all decisions inside the lock, before replacing current metadata.
        if sha256(cfg['_path']) != data['selection_config']['sha256'] or selection_inputs(root, cfg) != data['inputs']:
            raise ValueError('selection inputs changed while accepting')
        run = candidate_run(root, cfg, candidate, data)
        with tempfile.TemporaryDirectory(prefix='.accept-', dir=target.parent) as temporary:
            staged = Path(temporary) / target.name
            shutil.copyfile(candidate / 'metadata.tsv', staged)
            os.replace(staged, target)
        write_tsv(path_at(root, cfg['accepted_samples']), ACCEPTED_FIELDS, evidence)
        write_tsv(target.parent / 'selection.tsv', REPORT_FIELDS, read_tsv(candidate / 'selection.tsv'))
        data.update(kind='accepted_metadata', accepted_at=now(), candidate=str(candidate))
        data['metadata'] = file_record(target)
        data['accepted_samples'] = file_record(path_at(root, cfg['accepted_samples']))
        write_json(target.parent / 'provenance.json', data)
        if run is not None:
            run.update(status='complete', accepted_at=data['accepted_at'], metadata_history=metadata_history(data))
            if run.get('finished_at') is None or run.get('error'):
                run['finished_at'] = data['accepted_at']
            run.pop('error', None)
            write_json(candidate.parent / 'run.json', run)
    return {'samples': len(rows), 'retained_success_evidence': len(evidence)}


def container(root, cfg):
    from dataset_software import fetch_image, inspect_image, validate
    build = yaml.safe_load(path_at(root, cfg['build_config']).read_text())
    settings = dict(build['genegalleon'])
    validate(settings)
    settings['cache_dir'] = str(path_at(root, settings.get('cache_dir', 'resources/software/genegalleon')))
    runtime = shutil.which('apptainer') or shutil.which('singularity')
    if not runtime:
        raise ValueError('metadata preparation requires Apptainer or Singularity')
    if settings.get('image'):
        image = path_at(root, settings['image'])
        labels = inspect_image(runtime, image, settings)
    else:
        image, receipt = fetch_image(settings)
        labels = receipt['labels']
    identity = {'image': file_record(image), 'image_uri': settings['image_uri'], 'labels': labels}
    return runtime, Path(image), identity


def execute(runtime, image, root, arguments, binds=(), capture=False, log_path=None):
    paths = sorted({str(Path(root).resolve()), *(str(Path(p).resolve()) for p in binds)})
    command = [runtime, 'exec', '--cleanenv', '--env', 'PYTHONUNBUFFERED=1']
    for path in paths:
        if ':' in path or ',' in path:
            raise ValueError(f'container bind path contains unsupported separators: {path}')
        command += ['--bind', path + ':' + path]
    command += ['--pwd', str(Path(root).resolve()), str(image), *arguments]
    if capture:
        return subprocess.check_output(command, text=True)
    if log_path:
        with Path(log_path).open('w') as log:
            subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
    else:
        subprocess.run(command, check=True)
    return command


INSPECT_AMALGKIT = '''import importlib.metadata, json
from pathlib import Path
revision = None
p = Path('/opt/pg/logs/source_revisions.tsv')
if p.exists():
    for line in p.read_text().splitlines():
        fields = line.split('\\t')
        if fields[0] == 'amalgkit' and len(fields) == 2: revision = fields[1]
print(json.dumps({'version': importlib.metadata.version('amalgkit'), 'revision': revision}))
'''


def software_identity(root, cfg):
    runtime, image, identity = container(root, cfg)
    identity['amalgkit'] = json.loads(execute(runtime, image, root,
                                             ['python', '-c', INSPECT_AMALGKIT], capture=True))
    return runtime, image, identity


def export_rules(root, cfg, replace=False):
    rules = path_at(root, cfg['rules'])
    if rules.exists() and not replace:
        raise ValueError('rules already exist; --replace explicitly replaces the active rules')
    runtime, image, identity = software_identity(root, cfg)
    code = '''import importlib.resources, sys
print(importlib.resources.files('amalgkit.select_rule_sets').joinpath(sys.argv[1], 'select_rules.tsv').read_text(), end='')
'''
    text = execute(runtime, image, root, ['python', '-c', code, cfg['rule_set']], capture=True)
    if 'rule_id\t' not in text.splitlines()[0]:
        raise ValueError('container did not provide a select_rules.tsv rule set')
    rules.parent.mkdir(parents=True, exist_ok=True)
    rules.write_text(text)
    return {'rules': str(rules), 'software': identity}


def metadata_stage(root, cfg, action, work, source=None):
    threads = os.environ.get('SLURM_CPUS_PER_TASK')
    if threads is not None and (not threads.isdigit() or int(threads) < 1):
        raise ValueError('SLURM_CPUS_PER_TASK must be a positive integer')
    work = path_at(root, work)
    if not work.is_relative_to(Path(root).resolve() / 'work'):
        raise ValueError('metadata work directories must be under project work/')
    target = work / action
    if target.exists():
        raise ValueError(f'stage already exists; use a new work directory: {target}')
    inputs = {'selection_config': file_record(cfg['_path']),
              'build_config': file_record(path_at(root, cfg['build_config']))}
    runtime, image, identity = software_identity(root, cfg)
    if any(file_record(entry['path']) != entry for entry in inputs.values()):
        raise ValueError('metadata configuration changed while resolving the container')
    target.mkdir(parents=True)
    if action == 'fetch':
        arguments = ['amalgkit', 'metadata', '--out_dir', str(target), '--search_string', cfg['search_string']]
        binds = []
    else:
        source = path_at(root, source) if source else work / 'fetch/metadata/metadata.tsv'
        inputs['source_metadata'] = file_record(source)
        source_origin = {}
        source_provenance = source.parent / 'provenance.json'
        if source_provenance.exists():
            inputs['source_provenance'] = file_record(source_provenance)
            source_origin['stage'] = json.loads(source_provenance.read_text())
        query_info = source.parent / 'metadata.query_info.json'
        if query_info.exists():
            inputs['source_query_info'] = file_record(query_info)
            source_origin['query_info'] = json.loads(query_info.read_text())
        snapshot = target / 'source_metadata.tsv'
        shutil.copyfile(source, snapshot)
        if sha256(snapshot) != inputs['source_metadata']['sha256']:
            raise ValueError('source metadata changed while copying')
        rules = path_at(root, cfg['rules'])
        inputs['rules'] = file_record(rules)
        shutil.copyfile(rules, target / 'select_rules.tsv')
        if sha256(target / 'select_rules.tsv') != inputs['rules']['sha256']:
            raise ValueError('effective rules changed while copying')
        arguments = ['amalgkit', 'select', '--out_dir', str(target), '--metadata', str(snapshot),
                     '--select_rules_tsv', str(target / 'select_rules.tsv')]
        binds = []
    if threads is not None:
        arguments += ['--threads', threads]
    # select parameters, including sample_group, are supplied by select_rules.tsv.
    # The pinned AMALGKIT version does not accept a --sample_group CLI argument.
    provenance = {'schema_version': 1, 'kind': action, 'created_at': now(), 'software': identity,
                  'inputs': inputs, 'arguments': arguments, 'status': 'running',
                  'generator': file_record(Path(__file__))}
    if action == 'curate':
        provenance['source_origin'] = source_origin
    write_json(target / 'provenance.json', provenance)
    log_path = target / 'amalgkit.log'
    print(f'metadata preparation: {action}; log: {log_path}', file=sys.stderr, flush=True)
    try:
        command = execute(runtime, image, root, arguments, binds, log_path=log_path)
    except subprocess.CalledProcessError as error:
        provenance.update(status='failed', finished_at=now(), exit_code=error.returncode)
        write_json(target / 'provenance.json', provenance)
        raise ValueError(f'AMALGKIT {action} failed; see {log_path}') from error
    output = target / 'metadata/metadata.tsv'
    if not output.is_file():
        provenance.update(status='failed', finished_at=now(), error='expected metadata output missing')
        write_json(target / 'provenance.json', provenance)
        raise ValueError(f'AMALGKIT did not create {output}')
    provenance.update(status='complete', finished_at=now(), command=command, metadata=file_record(output))
    write_json(target / 'provenance.json', provenance)
    # Keep the provenance alongside the table for the candidate-selection step.
    write_json(output.parent / 'provenance.json', provenance)
    return {'metadata': str(output), 'log': str(log_path), 'software': identity}


def run_reference(work):
    return {'run_id': work.name, 'work': str(work)}


def candidate_run(root, cfg, candidate, provenance):
    """Only update the receipt belonging to this candidate's actual directory."""
    reference = provenance.get('metadata_run')
    if reference is None:
        return None
    work = candidate.parent
    if not work.is_relative_to(Path(root).resolve() / 'work') or reference != run_reference(work):
        raise ValueError('candidate metadata run does not match its work directory')
    run = json.loads((work / 'run.json').read_text())
    if (run.get('kind') != 'metadata_run' or run.get('dataset') != cfg['name']
            or run.get('work') != str(work) or run.get('candidate') != str(candidate)):
        raise ValueError('candidate does not belong to its metadata run')
    return run


def source_run(root, cfg, source):
    for directory in Path(source).parents:
        if not directory.is_relative_to(Path(root).resolve() / 'work'):
            break
        receipt = directory / 'run.json'
        if receipt.is_file():
            run = json.loads(receipt.read_text())
            if run.get('kind') == 'metadata_run' and run.get('dataset') == cfg['name']:
                return run_reference(directory)
    return None


def start_metadata_run(root, cfg, work, source, dry_run):
    root = Path(root).resolve()
    started_at = now()
    automatic = work is None
    if automatic:
        stamp = datetime.fromisoformat(started_at).astimezone(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
        work = Path(root) / 'work/datasets' / cfg['name'] / stamp
    work = work.parent.resolve() / work.name if automatic else path_at(root, work)
    if not work.is_relative_to(Path(root).resolve() / 'work'):
        raise ValueError('metadata work directories must be under project work/')
    run = {'schema_version': 1, 'kind': 'metadata_run', 'dataset': cfg['name'],
           **run_reference(work), 'mode': 'refresh' if source is None else 'reselect',
           'status': 'running', 'started_at': started_at, 'finished_at': None,
           'slurm_job_id': os.environ.get('SLURM_JOB_ID'), 'dry_run': dry_run,
           'selection_config': file_record(cfg['_path']), 'source_metadata': None,
           'source_run': None, 'candidate': str(work / 'candidate'), 'accepted_at': None}
    try:
        # mkdir reserves automatic names atomically, including across Slurm jobs.
        work.mkdir(parents=True, exist_ok=not automatic)
    except FileExistsError as error:
        raise ValueError(f'metadata work directory already exists: {work}; specify a new --work or run again later') from error
    if (work / 'candidate').exists():
        raise ValueError(f'candidate directory already exists; choose a new attempt: {work / "candidate"}')
    try:
        # Explicit work paths may contain separately prepared fetch/curate stages.
        # Claim their run receipt without overwriting a previous update.
        with (work / 'run.json').open('x') as handle:
            json.dump(run, handle, indent=2, ensure_ascii=False, allow_nan=False)
            handle.write('\n')
    except FileExistsError as error:
        raise ValueError(f'metadata run already exists; choose a new --work: {work}') from error
    return work, run


def update_metadata(root, cfg, work=None, source=None, dry_run=False):
    """Fetch/curate/select/publish in one command, or reuse a curated download."""
    work, run = start_metadata_run(root, cfg, work, source, dry_run)
    candidate = work / 'candidate'
    print(f'metadata work: {work}\nmetadata candidate: {candidate}', file=sys.stderr, flush=True)
    try:
        if source is None:
            metadata_stage(root, cfg, 'fetch', work)
            curated = metadata_stage(root, cfg, 'curate', work)
            source = curated['metadata']
        source = path_at(root, source)
        run['source_metadata'] = file_record(source)
        if run['mode'] == 'reselect':
            run['source_run'] = source_run(root, cfg, source)
            if run['source_run'] == run_reference(work):
                run['source_run'] = None
        write_json(work / 'run.json', run)
        summary = select_candidates(root, cfg, source, candidate)
        provenance = json.loads((candidate / 'provenance.json').read_text())
        provenance['metadata_run'] = run_reference(work)
        write_json(candidate / 'provenance.json', provenance)
        run.update(summary=summary, source_metadata=provenance['source_metadata'],
                   metadata_history=metadata_history(provenance))
        write_json(work / 'run.json', run)
        result = {**run_reference(work), 'candidate': str(candidate), 'summary': summary, 'dry_run': dry_run}
        unresolved = summary['unresolved_accepted_runs']
        if not dry_run and not unresolved:
            result['updated'] = accept_candidate(root, cfg, candidate)
            run = json.loads((work / 'run.json').read_text())
        run.update(status='review_required' if unresolved else 'complete', finished_at=now())
        write_json(work / 'run.json', run)
        if unresolved:
            result['status'] = 'review_required'
            print('metadata review required: ' + ', '.join(unresolved)
                  + '; restore their metadata or explicitly exclude/replace them, then rerun update',
                  file=sys.stderr, flush=True)
        elif dry_run:
            result['accept_command'] = shlex.join([str(Path(root).resolve() / 'run_metadata.sh'), 'accept',
                '--config', str(cfg['_path']), '--candidate', str(candidate)])
            print(f'accept candidate: {result["accept_command"]}', file=sys.stderr, flush=True)
        return result
    except (Exception, KeyboardInterrupt) as error:
        run.update(status='failed', finished_at=now(), error=str(error))
        write_json(work / 'run.json', run)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    for name in ('update', 'record', 'init', 'seed', 'baseline', 'rules', 'fetch', 'curate', 'select', 'accept'):
        p = commands.add_parser(name)
        p.add_argument('--root', default='.')
        p.add_argument('--config', default='datasets/angiosperm_leaf/selection.yaml')
        if name == 'update':
            p.add_argument('--work', help='Work directory; default: work/datasets/<name>/<UTC timestamp>')
            p.add_argument('--metadata', help='Optional already curated AMALGKIT table; otherwise fetch and curate')
            p.add_argument('--dry-run', action='store_true', help='Write a candidate without publishing dataset metadata')
        elif name == 'record':
            source = p.add_mutually_exclusive_group(required=True)
            source.add_argument('--build', help='Build to review')
            source.add_argument('--apply', metavar='PLAN.tsv', help='Apply a saved, optionally edited review plan')
            p.add_argument('--runs', nargs='+', help='Reviewed run accessions; default: all current representatives in the build')
            p.add_argument('--dry-run', action='store_true', help='Preview decisions and cleanup without changing tables or deleting files')
            p.add_argument('--exclude-failed', action='store_true', help='Exclude all inactive incomplete samples, including unstarted runs, and remove their reads/scratch; --dry-run previews without applying')
            p.add_argument('--plan', metavar='PLAN.tsv', help='Save an editable review plan and evidence; does not apply decisions')
        elif name == 'init':
            p.add_argument('--exclusions', help='Optional historical run exclusion table')
            p.add_argument('--replace', action='store_true')
        elif name == 'seed':
            p.add_argument('--metadata', required=True)
            p.add_argument('--exclusions', required=True)
            p.add_argument('--source-build', default='historical')
            p.add_argument('--replace', action='store_true')
        elif name == 'baseline':
            p.add_argument('--database', required=True)
            p.add_argument('--source-metadata', required=True)
            p.add_argument('--exclusions', required=True)
            p.add_argument('--replace', action='store_true')
        elif name == 'rules':
            p.add_argument('--replace', action='store_true', help='Replace the active rules with the pinned container rule set')
        elif name in ('fetch', 'curate'):
            p.add_argument('--work', required=True)
            if name == 'curate': p.add_argument('--metadata')
        elif name == 'select':
            p.add_argument('--metadata', required=True, help='Already curated AMALGKIT metadata')
            p.add_argument('--output', required=True, help='New candidate directory under work/')
        else:
            p.add_argument('--candidate', required=True)
            p.add_argument('--allow-review', action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args()
    root = Path(args.root).resolve()
    try:
        cfg = configuration(root, args.config)
        if args.command == 'update':
            result = update_metadata(root, cfg, args.work, args.metadata, args.dry_run)
        elif args.command == 'record':
            if args.apply:
                if args.runs or args.exclude_failed or args.plan:
                    raise ValueError('--apply cannot be combined with --runs, --exclude-failed or --plan')
                from sample_review import apply_plan
                result = apply_plan(root, cfg, args.apply, dry_run=args.dry_run)
            else:
                result = record_successes(root, cfg, args.build, args.runs, args.dry_run,
                                          exclude_failed=args.exclude_failed, plan=args.plan)
        elif args.command == 'init':
            result = initialize_empty(root, cfg, args.exclusions, args.replace)
        elif args.command == 'seed':
            result = seed_history(root, cfg, args.metadata, args.exclusions, args.source_build, args.replace)
        elif args.command == 'baseline':
            result = initialize(root, cfg, args.database, args.source_metadata, args.exclusions, args.replace)
        elif args.command == 'rules': result = export_rules(root, cfg, args.replace)
        elif args.command in ('fetch', 'curate'):
            result = metadata_stage(root, cfg, args.command, args.work, getattr(args, 'metadata', None))
        elif args.command == 'select': result = select_candidates(root, cfg, args.metadata, args.output)
        else: result = accept_candidate(root, cfg, args.candidate, args.allow_review)
        print(json.dumps(result, indent=2))
    except (ValueError, OSError, KeyError, subprocess.CalledProcessError) as error:
        parser.exit(1, f'metadata preparation: {error}\n')


if __name__ == '__main__':
    main()
