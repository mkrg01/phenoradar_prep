#!/usr/bin/env python3
"""Curated metadata candidates, independent of database job submission."""
import argparse
import collections
import csv
from decimal import Decimal, InvalidOperation
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

import yaml

from accession_exclusions import read_exclusions
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
    fields, rows = history_table(path_at(root, cfg['previous_metadata']),
                                ['scientific_name', 'taxid', 'run', 'bioproject'])
    previous = unique_by(rows, 'taxid')
    unique_by(rows, 'run')
    _, accepted_rows = history_table(path_at(root, cfg['accepted_samples']), ACCEPTED_FIELDS)
    accepted = unique_by(accepted_rows, 'taxid')
    for tid, row in accepted.items():
        prior = previous.get(tid)
        if not prior or any(prior[key] != row[key] for key in ('run', 'bioproject', 'scientific_name')):
            raise ValueError(f'accepted sample differs from previous metadata: taxid {tid}')
        complete, total = int(row['busco_complete']), int(row['busco_total'])
        if total <= 0 or not 0 <= complete <= total:
            raise ValueError(f'invalid accepted BUSCO counts: taxid {tid}')
    return fields, previous, accepted


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


def record_successes(root, cfg, build_path, runs=None, dry_run=False):
    """Adopt completed samples and exclude new BUSCO failures from partial builds."""
    from dataset import load, item_products
    build_path = path_at(root, build_path)
    if build_path.name == 'build.json':
        build_path = build_path.parent
    source = file_record(build_path / 'build.json')
    build = load(build_path)
    matching_container(root, cfg, build)
    before, config_record = selection_inputs(root, cfg), file_record(cfg['_path'])
    _, previous, accepted = previous_state(root, cfg)
    exclusion_path = path_at(root, cfg['excluded_accessions'])
    exclusions = read_exclusions(exclusion_path)
    items = {item['row']['run']: item for item in build['items']}
    requested = set(runs) if runs is not None else set(items)
    if requested - items.keys():
        raise ValueError('runs absent from build: ' + ', '.join(sorted(requested - items.keys())))
    accepted = {tid: row for tid, row in accepted.items() if row['run'] not in exclusions}
    report, recorded, decisions = [], [], []
    for run in sorted(requested):
        item = items[run]
        row = item['row']
        tid = taxid(row['taxid'])
        current = previous.get(tid)
        status = 'excluded' if run in exclusions else ''
        if not status and (not current or any(current[k] != row[k] for k in ACCEPTED_FIELDS[:4])):
            status = 'not_current_representative'
        if status:
            report.append({'run': run, 'status': status})
            continue
        products = item_products(build, item)
        if any(not products[key] for key in ('reference', 'busco')):
            report.append({'run': run, 'status': 'incomplete'})
            continue
        qc = counts(products['busco']['counts'])
        evidence = {**{k: row[k] for k in ACCEPTED_FIELDS[:4]},
                    'busco_complete': qc['busco_cds_single'] + qc['busco_cds_duplicated'],
                    'busco_total': qc['busco_cds_total'], 'source_build': build['name']}
        if evidence['busco_complete'] / evidence['busco_total'] < cfg['busco_threshold']:
            if tid not in accepted:
                decisions.append(busco_exclusion(row, cfg['busco_threshold'], build['name']))
            report.append({'run': run, 'status': 'busco_below_threshold',
                           'busco_complete': evidence['busco_complete'], 'busco_total': evidence['busco_total'],
                           'retained_accepted': tid in accepted})
            continue
        if not products['quant']:
            report.append({'run': run, 'status': 'incomplete'})
            continue
        accepted[tid] = evidence
        recorded.append(run)
        report.append({'run': run, 'status': 'reviewed_complete',
                       'busco_complete': evidence['busco_complete'], 'busco_total': evidence['busco_total']})
    excluded_runs = [decision['accession'] for decision in decisions]
    result = {'recorded_runs': recorded, 'excluded_runs': excluded_runs,
              'accepted_samples': len(accepted), 'samples': report,
              'dry_run': dry_run, 'busco_threshold': cfg['busco_threshold']}
    if dry_run:
        return result
    target = path_at(root, cfg['accepted_samples'])
    with locked(path_at(root, cfg['previous_metadata']).parent / '.metadata.lock'):
        if (file_record(cfg['_path']) != config_record or selection_inputs(root, cfg) != before
                or file_record(build_path / 'build.json') != source):
            raise ValueError('sample adoption inputs changed; repeat the review')
        append_exclusions(exclusion_path, decisions)
        write_tsv(target, ACCEPTED_FIELDS, [accepted[tid] for tid in sorted(accepted, key=int)])
        write_json(target.with_suffix('.provenance.json'), {'kind': 'reviewed_samples', 'created_at': now(),
                   'build_manifest': source, 'selection_config': config_record, 'inputs': before,
                   'recorded_runs': recorded, 'busco_threshold': cfg['busco_threshold'],
                   'excluded_runs': excluded_runs, 'samples': report,
                   'excluded_accessions': file_record(exclusion_path),
                   'accepted_samples': file_record(target)})
    return result


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
    kept = {tid: previous[tid] for tid, row in accepted.items()
            if row['run'] not in exclusions}
    best, forced, observed, candidate_counts = {}, {}, {}, collections.Counter()
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
            if tid in kept and run == kept[tid]['run']:
                observed[tid] = reason or ('identity_changed' if project != kept[tid]['bioproject'] else '')
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
            # Successful prior rows are kept verbatim, even when a larger run appears.
            if tid in overrides and run == overrides[tid]['run']:
                forced[tid] = normalized_row(dict(zip(source_fields, values)))
            if tid in kept or tid in overrides:
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
    reports = []
    for tid in sorted(set(previous) | observed_taxids | set(selected), key=int):
        old, new = previous.get(tid), selected.get(tid)
        review, reason = False, ''
        if tid in overrides:
            decision, reason = 'override', overrides[tid]['reason']
        elif tid in kept:
            decision, reason = 'retained', 'previous_completed_sample'
            if tid not in observed:
                review, reason = True, 'retained_run_absent_from_source'
            elif observed[tid]:
                review, reason = True, 'retained_run_source_status:' + observed[tid]
        elif new:
            decision = ('reselected' if old['run'] == new['run'] else 'replacement') if old else 'new'
            reason = 'largest_total_bases'
        else:
            decision, reason = 'no_candidate', 'no_eligible_run'
        reports.append({'taxid': tid, 'scientific_name': (new or old or {}).get('scientific_name', ''),
                        'previous_run': old['run'] if old else '', 'selected_run': new['run'] if new else '',
                        'bioproject': new['bioproject'] if new else '', 'decision': decision,
                        'previous_bioproject': old.get('bioproject', '') if old else '',
                        'previous_total_bases': old.get('total_bases', '') if old else '',
                        'selected_total_bases': new.get('total_bases', '') if new else '',
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
                   'review_required': sum(r['review'] == 'yes' for r in reports), 'source_rejections': dict(reasons)}
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


def snapshot_rules(root, cfg, refresh=False):
    runtime, image, identity = software_identity(root, cfg)
    code = '''import importlib.resources, sys
print(importlib.resources.files('amalgkit.select_rule_sets').joinpath(sys.argv[1], 'select_rules.tsv').read_text(), end='')
'''
    text = execute(runtime, image, root, ['python', '-c', code, cfg['rule_set']], capture=True)
    if 'rule_id\t' not in text.splitlines()[0]:
        raise ValueError('container did not provide a select_rules.tsv rule set')
    effective = path_at(root, cfg['rules'])
    upstream = effective.parent / 'upstream/select_rules.tsv'
    if upstream.exists() and not refresh:
        raise ValueError('upstream rules already exist; --refresh-upstream explicitly refreshes the snapshot')
    upstream.parent.mkdir(parents=True, exist_ok=True)
    upstream.write_text(text)
    if not effective.exists():
        shutil.copyfile(upstream, effective)
    write_json(upstream.with_suffix('.json'), {'created_at': now(), 'rule_set': cfg['rule_set'],
               'software': identity, 'rules': file_record(upstream)})
    return {'upstream': str(upstream), 'effective_rules_preserved': str(effective), 'software': identity}


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


def update_metadata(root, cfg, work, source=None, dry_run=False):
    """Fetch/curate/select/publish in one command, or reuse a curated download."""
    work = path_at(root, work)
    if source is None:
        metadata_stage(root, cfg, 'fetch', work)
        curated = metadata_stage(root, cfg, 'curate', work)
        source = curated['metadata']
    candidate = work / 'candidate'
    summary = select_candidates(root, cfg, source, candidate)
    result = {'candidate': str(candidate), 'summary': summary, 'dry_run': dry_run}
    if not dry_run:
        result['updated'] = accept_candidate(root, cfg, candidate)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    for name in ('update', 'record', 'init', 'seed', 'baseline', 'rules', 'fetch', 'curate', 'select', 'accept'):
        p = commands.add_parser(name)
        p.add_argument('--root', default='.')
        p.add_argument('--config', default='datasets/angiosperm_leaf/selection.yaml')
        if name == 'update':
            p.add_argument('--work', required=True)
            p.add_argument('--metadata', help='Optional already curated AMALGKIT table; otherwise fetch and curate')
            p.add_argument('--dry-run', action='store_true', help='Write a candidate without publishing dataset metadata')
        elif name == 'record':
            p.add_argument('--build', required=True)
            p.add_argument('--runs', nargs='+', help='Reviewed run accessions; default: all current representatives in the build')
            p.add_argument('--dry-run', action='store_true', help='Preview adoption and BUSCO exclusions without changing tables')
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
            p.add_argument('--refresh-upstream', action='store_true')
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
            result = record_successes(root, cfg, args.build, args.runs, args.dry_run)
        elif args.command == 'init':
            result = initialize_empty(root, cfg, args.exclusions, args.replace)
        elif args.command == 'seed':
            result = seed_history(root, cfg, args.metadata, args.exclusions, args.source_build, args.replace)
        elif args.command == 'baseline':
            result = initialize(root, cfg, args.database, args.source_metadata, args.exclusions, args.replace)
        elif args.command == 'rules': result = snapshot_rules(root, cfg, args.refresh_upstream)
        elif args.command in ('fetch', 'curate'):
            result = metadata_stage(root, cfg, args.command, args.work, getattr(args, 'metadata', None))
        elif args.command == 'select': result = select_candidates(root, cfg, args.metadata, args.output)
        else: result = accept_candidate(root, cfg, args.candidate, args.allow_review)
        print(json.dumps(result, indent=2))
    except (ValueError, OSError, KeyError, subprocess.CalledProcessError) as error:
        parser.exit(1, f'metadata preparation: {error}\n')


if __name__ == '__main__':
    main()
