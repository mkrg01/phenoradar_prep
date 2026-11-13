"""Publish a filtered copy of a portable database from a completed taxonomy screen."""
import csv
import json
import os
from pathlib import Path
import shutil
import tempfile

from accession_exclusions import read_exclusions
from build_products import load_complete
from common import now, read_tsv, write_json, write_tsv
from dataset_assets import digest, locked, record, stage_file, verify
from mapping_tables import load_tables, subset
from portable_build import PRODUCT_FILES, completion_path


def filter_database(database, screen, output):
    source = record(completion_path(database))
    data = load_complete(source['path'], verify_files=False)
    evidence = record(screen)
    qc = json.loads(verify(evidence).read_text())
    if (qc.get('kind') != 'angiosperms353_qc' or qc.get('schema_version') != 3
            or not qc.get('database_manifest') or qc['database_manifest']['sha256'] != source['sha256']
            or qc['settings']['restricted_loci'] or qc['settings']['omit_same_species']):
        raise ValueError('filtering requires a completed production screen of the same portable database')
    by_run = {product['row']['run']: product for product in data['products'].values()}
    if len(by_run) != len(data['products']) or set(qc['samples']) != by_run.keys():
        raise ValueError('the screen must cover every sample in the portable database')
    for run, observed in qc['samples'].items():
        product = by_run[run]
        if observed['taxid'] != str(product['row']['taxid']) or observed['cds']['sha256'] != product['cds']['sha256']:
            raise ValueError(f'screen sample identity differs from the source database: {run}')
    decisions = read_exclusions(verify(qc['exclusion_candidates']))
    summaries = read_tsv(verify(qc['summary']))
    if (len(summaries) != len(by_run) or {row['run'] for row in summaries} != by_run.keys()
            or {row['run'] for row in summaries if row['status'] == 'discordant'} != decisions.keys()):
        raise ValueError('exclusion candidates do not agree with the completed screen')
    selected = {name: product for name, product in data['products'].items() if product['row']['run'] not in decisions}
    if not selected:
        raise ValueError('the taxonomy screen excludes every sample')
    output = Path(output).resolve()
    source_root = Path(data['bundle_root'])
    if output.is_relative_to(source_root):
        raise ValueError('write the filtered database outside the input bundle')
    provenance = dict(source_manifest_sha256=source['sha256'], screen_sha256=evidence['sha256'],
                      samples=len(by_run), retained=len(selected), excluded=len(decisions),
                      max_registered_support=qc['settings']['max_registered_support'], min_loci=qc['settings']['min_loci'])
    with locked(output.with_name('.' + output.name + '.lock')):
        if output.exists():
            existing = load_complete(output)
            if existing.get('taxonomy_filter') != provenance:
                raise ValueError(f'existing database belongs to a different filter run: {output}')
        else:
            publish(output, source, evidence, qc, data, selected, decisions, provenance)
    return dict(output=str(output), **provenance)


def publish(output, source, evidence, qc, data, selected, decisions, provenance):
    original = Path(data['bundle_root'])
    staging = Path(tempfile.mkdtemp(prefix='.' + output.name + '-', dir=output.parent))
    files = {}

    def register(relative, saved):
        files[str(relative)] = dict(saved, path=str(relative))
        return files[str(relative)]

    def copy(entry, relative):
        return register(relative, stage_file(entry, staging / relative))

    def table(relative, fields, rows):
        write_tsv(staging / relative, fields, rows)
        register(relative, record(staging / relative))

    try:
        products = {}
        for name, product in selected.items():
            products[name] = dict(product, **{
                key: copy(product[key], Path(product[key]['path']).relative_to(original))
                for key in PRODUCT_FILES})
        for entry in data['files']:
            relative = Path(entry['path']).relative_to(original)
            if relative.parts[0] == 'provenance' or relative == Path('source_metadata.tsv'):
                copy(entry, relative)
        for name, entry in [('taxonomy_screen.json', evidence), ('taxonomy_summary.tsv', qc['summary']),
                            ('taxonomy_exclusion_candidates.tsv', qc['exclusion_candidates'])]:
            copy(entry, Path('provenance') / name)
        excluded = [*data['excluded_runs'],
                    *(dict(species=name, run=product['row']['run'],
                           reason=decisions[product['row']['run']]['reason'])
                      for name, product in data['products'].items() if product['row']['run'] in decisions)]
        indexed = {entry['path']: entry for entry in data['files']}
        policy = read_exclusions(verify(indexed[str(original / 'excluded_accessions.tsv')]))
        policy.update(decisions)
        fields = list(dict.fromkeys(['accession', 'reason', *(field for row in policy.values() for field in row)]))
        table('excluded_accessions.tsv', fields, policy.values())
        table('excluded_runs.tsv', ['species', 'run', 'reason'],
              ({key: row.get(key, '') for key in ('species', 'run', 'reason')} for row in excluded))
        with verify(indexed[str(original / 'metadata.tsv')]).open(newline='') as handle:
            rows = csv.DictReader(handle, delimiter='\t')
            retained = [row for row in rows if row['run'] not in decisions]
            if len(retained) != len(selected) or {row['run'] for row in retained} != {
                    product['row']['run'] for product in selected.values()}:
                raise ValueError('source metadata does not match the retained products')
            table('metadata.tsv', rows.fieldnames, retained)
        from dataset_assets import COUNTS
        table('busco/summary.tsv', ['Species', *COUNTS],
              (dict(Species=name, **product['counts']) for name, product in products.items()))
        mapping_path = subset(verify(data['mapping']), products, staging / 'odb')
        mapping = load_tables(mapping_path, verify_files=False)
        for entry in mapping['tables'].values():
            register('odb/' + entry['table']['path'], entry['table'])
        from sample_processing import write_summary
        write_summary(staging / 'sample_processing.tsv', dict(data, products=products))
        register('sample_processing.tsv', record(staging / 'sample_processing.tsv'))
        verify(source)
        verify(evidence)
        stamp = now()
        final = {key: value for key, value in data.items()
                 if key not in {'sha256', 'input', 'bundle_root', 'metadata_history', 'prepared_at'}}
        mapping_record = register('odb/snapshot.json', record(mapping_path))
        final.update(build_id=output.name, created_at=stamp, completed_at=stamp, dataset_name=output.name,
                     name_mode='fixed', products=products, excluded_runs=excluded, files=list(files.values()),
                     mapping=mapping_record, taxonomy_filter=provenance)
        final['sha256'] = digest(final)
        write_json(staging / 'manifest.json', final)
        load_complete(staging)
        os.rename(staging, output)
    finally:
        if staging.exists():
            shutil.rmtree(staging)
