"""Representative updates preserve successful runs and remain separate from builds."""
import json
from pathlib import Path
import subprocess

import pytest
import yaml

import metadata_catalog as catalog
from common import file_record, read_tsv, write_json, write_tsv
from sample_identity import annotate


FIELDS = ['scientific_name', 'taxid', 'run', 'bioproject', 'sample_group',
          'exclusion', 'total_bases', 'lib_layout']


def row(tid, run, project='P1', bases=100, **extra):
    return dict(scientific_name=f'Plant species{tid}', taxid=str(tid), run=run,
                bioproject=project, sample_group='leaf', exclusion='no',
                total_bases=str(bases), lib_layout='paired', **extra)


@pytest.fixture
def curated_project(tmp_path):
    dataset = tmp_path / 'datasets/leaf'
    (dataset / 'rules').mkdir(parents=True)
    (dataset / 'rules/select_rules.tsv').write_text('rule_id\tenabled\nleaf\tyes\n')
    cfg = {'schema_version': 1, 'name': 'leaf', 'build_config': 'datasets/leaf/build.yaml',
           'search_string': 'plants', 'sample_group': 'leaf', 'rule_set': 'plantae',
           'rules': 'datasets/leaf/rules/select_rules.tsv',
           'accepted_samples': 'datasets/leaf/accepted_samples.tsv',
           'previous_metadata': 'datasets/leaf/metadata.tsv',
           'excluded_accessions': 'datasets/leaf/excluded_accessions.tsv',
           'overrides': 'datasets/leaf/overrides.tsv', 'busco_threshold': 0.5}
    config = dataset / 'selection.yaml'
    config.write_text(yaml.safe_dump(cfg))
    (dataset / 'build.yaml').write_text(yaml.safe_dump({'genegalleon': {
        'version': '0.7.77', 'revision': 'a' * 40,
        'image_uri': 'docker://example/image@sha256:' + 'a' * 64}}))
    previous = [annotate(row(1, 'OLD1')), annotate(row(2, 'BAD2'))]
    write_tsv(dataset / 'metadata.tsv', list(previous[0]), previous)
    write_tsv(dataset / 'accepted_samples.tsv', catalog.ACCEPTED_FIELDS,
              [dict(taxid='1', run='OLD1', bioproject='P1', scientific_name='Plant species1',
                    busco_complete=80, busco_total=100, source_build='completed')])
    write_tsv(dataset / 'excluded_accessions.tsv', catalog.EXCLUSION_FIELDS,
              [dict(accession='BAD2', taxid='2', bioproject='P1', reason='poor_quality', source_build='old')])
    write_tsv(dataset / 'overrides.tsv', ['taxid', 'run', 'reason'], [])
    cfg = catalog.configuration(tmp_path, config)
    return tmp_path, dataset, cfg


@pytest.fixture
def fresh_project(curated_project):
    root, dataset, cfg = curated_project
    (dataset / 'metadata.tsv').unlink()
    (dataset / 'accepted_samples.tsv').unlink()
    return root, dataset, cfg


def select(project, rows, name='attempt'):
    root, dataset, cfg = project
    source = root / f'{name}-source.tsv'
    fields = list(dict.fromkeys(k for r in rows for k in r))
    write_tsv(source, fields, [{k: r.get(k, '') for k in fields} for r in rows])
    output = root / 'work' / name
    summary = catalog.select_candidates(root, cfg, source, output)
    return summary, output


def test_first_selection_and_acceptance_work_without_placeholder_tables(fresh_project):
    root, dataset, cfg = fresh_project
    exclusions = (dataset / 'excluded_accessions.tsv').read_bytes()
    summary, output = select(fresh_project, [row(1, 'OLD1'), row(1, 'BIG1', bases=5000),
                                            row(2, 'RELATED2'), row(2, 'GOOD2', project='P2')])
    assert summary['decisions'] == {'new': 2}
    assert {r['run'] for r in read_tsv(output / 'metadata.tsv')} == {'BIG1', 'GOOD2'}
    for name in ('metadata.tsv', 'accepted_samples.tsv', 'selection.tsv', 'provenance.json'):
        assert not (dataset / name).exists()
    assert catalog.accept_candidate(root, cfg, output) == {
        'samples': 2, 'retained_success_evidence': 0}
    assert {r['run'] for r in read_tsv(dataset / 'metadata.tsv')} == {'BIG1', 'GOOD2'}
    assert not read_tsv(dataset / 'accepted_samples.tsv')
    assert (dataset / 'excluded_accessions.tsv').read_bytes() == exclusions
    assert not (root / 'results').exists()


@pytest.mark.parametrize('name', ['metadata.tsv', 'accepted_samples.tsv'])
def test_history_created_after_first_selection_requires_regeneration(fresh_project, name):
    root, dataset, cfg = fresh_project
    _, output = select(fresh_project, [row(1, 'NEW1')])
    fields = catalog.ACCEPTED_FIELDS if name == 'accepted_samples.tsv' else FIELDS
    write_tsv(dataset / name, fields, [])
    with pytest.raises(ValueError, match='changed'):
        catalog.accept_candidate(root, cfg, output)
    assert not (dataset / 'provenance.json').exists()


def test_success_evidence_cannot_be_used_without_previous_metadata(curated_project):
    (curated_project[1] / 'metadata.tsv').unlink()
    with pytest.raises(ValueError, match='accepted sample differs from previous metadata'):
        select(curated_project, [row(1, 'OLD1')])
    assert not (curated_project[0] / 'work').exists()


def test_retains_accepted_run_and_excludes_only_listed_accession(curated_project):
    rows = [row(1, 'OLD1'), row(1, 'BIGGER1', bases=1000000),
            row(2, 'BAD2', bases=1000000), row(2, 'SAME_PROJECT2', bases=900000),
            row(2, 'OTHER_PROJECT2', project='P2', bases=300),
            row(3, 'OTHER_TAXID_SAME_PROJECT', bases=900000)]
    old = (curated_project[1] / 'metadata.tsv').read_bytes()
    summary, output = select(curated_project, rows)
    selected = {r['taxid']: r for r in read_tsv(output / 'metadata.tsv')}
    assert selected['1']['run'] == 'OLD1'
    assert selected['2']['run'] == 'SAME_PROJECT2'
    assert selected['3']['run'] == 'OTHER_TAXID_SAME_PROJECT'
    assert summary['decisions'] == {'retained': 1, 'replacement': 1, 'new': 1}
    assert summary['source_rejections'] == {'excluded_run': 1}
    assert (curated_project[1] / 'metadata.tsv').read_bytes() == old
    assert not (curated_project[0] / 'results').exists()
    assert selected['1']['scientific_name'] == 'Plant species1'


def test_ties_are_order_independent_and_missing_bases_are_skipped(curated_project):
    rows = [row(1, 'OLD1'), row(3, 'Z3'), row(3, 'A3'), row(4, 'MISSING4', bases='NA')]
    _, first = select(curated_project, rows, 'first')
    _, second = select(curated_project, list(reversed(rows)), 'second')
    assert (first / 'metadata.tsv').read_bytes() == (second / 'metadata.tsv').read_bytes()
    assert {r['taxid']: r['run'] for r in read_tsv(first / 'metadata.tsv')} == {'1': 'OLD1', '3': 'A3'}
    assert any(r['taxid'] == '4' and r['decision'] == 'no_candidate' for r in read_tsv(first / 'selection.tsv'))


def test_new_exclusion_overrides_previous_success(curated_project):
    root, dataset, cfg = curated_project
    rows = read_tsv(dataset / 'excluded_accessions.tsv')
    rows.append(dict(accession='OLD1', taxid='1', bioproject='P1', reason='misidentified', source_build='old'))
    write_tsv(dataset / 'excluded_accessions.tsv', catalog.EXCLUSION_FIELDS, rows)
    _, output = select(curated_project, [row(1, 'OLD1'), row(1, 'RELATED1', bases=500), row(1, 'GOOD1', project='P2')])
    selected = read_tsv(output / 'metadata.tsv')
    assert [r['run'] for r in selected] == ['RELATED1']
    catalog.accept_candidate(root, cfg, output)
    assert not read_tsv(dataset / 'accepted_samples.tsv')


def test_absent_success_is_preserved_and_reported_without_blocking_update(curated_project):
    root, dataset, cfg = curated_project
    summary, output = select(curated_project, [row(1, 'NEW1', bases=5000)])
    assert summary['review_required'] == 1
    assert read_tsv(output / 'metadata.tsv')[0]['run'] == 'OLD1'
    catalog.accept_candidate(root, cfg, output)
    assert read_tsv(dataset / 'metadata.tsv')[0]['run'] == 'OLD1'


def test_acceptance_keeps_success_evidence_only_for_old_representatives(curated_project):
    root, dataset, cfg = curated_project
    _, output = select(curated_project, [row(1, 'OLD1'), row(2, 'NEW2', project='P2')])
    result = catalog.accept_candidate(root, cfg, output)
    assert result == {'samples': 2, 'retained_success_evidence': 1}
    assert [r['run'] for r in read_tsv(dataset / 'accepted_samples.tsv')] == ['OLD1']
    assert not (root / 'results').exists()
    assert json.loads((dataset / 'provenance.json').read_text())['kind'] == 'accepted_metadata'


@pytest.mark.parametrize('changed', ['rules/select_rules.tsv', 'overrides.tsv', 'metadata.tsv', 'build.yaml'])
def test_inputs_changed_after_selection_cannot_be_accepted(curated_project, changed):
    root, dataset, cfg = curated_project
    _, output = select(curated_project, [row(1, 'OLD1')])
    with (dataset / changed).open('a') as f: f.write('\n')
    with pytest.raises(ValueError, match='changed'):
        catalog.accept_candidate(root, cfg, output)


def test_edited_candidate_requires_regeneration(curated_project):
    root, dataset, cfg = curated_project
    _, output = select(curated_project, [row(1, 'OLD1')])
    with (output / 'metadata.tsv').open('a') as f: f.write('\n')
    with pytest.raises(ValueError, match='candidate metadata changed'):
        catalog.accept_candidate(root, cfg, output)


def test_override_requires_eligible_run_and_keeps_reason(curated_project):
    root, dataset, cfg = curated_project
    write_tsv(dataset / 'overrides.tsv', ['taxid', 'run', 'reason'],
              [dict(taxid='1', run='MANUAL1', reason='verified tissue')])
    _, output = select(curated_project, [row(1, 'OLD1'), row(1, 'MANUAL1', project='P2')])
    assert read_tsv(output / 'metadata.tsv')[0]['run'] == 'MANUAL1'
    assert read_tsv(output / 'selection.tsv')[0]['reason'] == 'verified tissue'
    catalog.accept_candidate(root, cfg, output)
    assert not read_tsv(dataset / 'accepted_samples.tsv')
    write_tsv(dataset / 'overrides.tsv', ['taxid', 'run', 'reason'],
              [dict(taxid='2', run='BAD2', reason='manual')])
    with pytest.raises(ValueError, match='override is absent or excluded'):
        select(curated_project, [row(2, 'BAD2')], 'blocked')


def test_accepted_run_is_retained_independently_of_busco_completeness(curated_project):
    root, dataset, cfg = curated_project
    evidence = read_tsv(dataset / 'accepted_samples.tsv')
    evidence[0]['busco_complete'] = '40'
    write_tsv(dataset / 'accepted_samples.tsv', catalog.ACCEPTED_FIELDS, evidence)
    _, output = select(curated_project, [row(1, 'OLD1'), row(1, 'LARGER1', bases=500)])
    assert read_tsv(output / 'metadata.tsv')[0]['run'] == 'OLD1'


def test_source_duplicates_fail_and_exclusion_context_is_optional(curated_project):
    with pytest.raises(ValueError, match='duplicate source run'):
        select(curated_project, [row(1, 'OLD1'), row(1, 'OLD1')])
    root, dataset, cfg = curated_project
    decisions = read_tsv(dataset / 'excluded_accessions.tsv')
    decisions[0]['bioproject'] = ''
    write_tsv(dataset / 'excluded_accessions.tsv', catalog.EXCLUSION_FIELDS, decisions)
    _, output = select(curated_project, [row(1, 'OLD1'), row(2, 'BAD2'), row(2, 'RELATED2')], 'missing-project')
    assert {r['run'] for r in read_tsv(output / 'metadata.tsv')} == {'OLD1', 'RELATED2'}


def test_container_rule_refresh_preserves_user_edits(curated_project, monkeypatch):
    root, dataset, cfg = curated_project
    effective = dataset / 'rules/select_rules.tsv'
    before = effective.read_bytes()
    monkeypatch.setattr(catalog, 'software_identity', lambda *a: ('runtime', Path('image'), {'amalgkit': {'version': 'fixed'}}))
    monkeypatch.setattr(catalog, 'execute', lambda *a, **k: 'rule_id\tenabled\nupstream\tyes\n')
    catalog.snapshot_rules(root, cfg)
    assert effective.read_bytes() == before
    upstream = dataset / 'rules/upstream/select_rules.tsv'
    assert 'upstream' in upstream.read_text()
    with pytest.raises(ValueError, match='already exist'):
        catalog.snapshot_rules(root, cfg)
    catalog.snapshot_rules(root, cfg, refresh=True)
    assert effective.read_bytes() == before


def test_curation_copies_input_and_rules_without_running_builds(curated_project, monkeypatch):
    root, dataset, cfg = curated_project
    source = root / 'source.tsv'
    write_tsv(source, FIELDS, [row(1, 'OLD1')])
    original = source.read_bytes()
    monkeypatch.setattr(catalog, 'software_identity', lambda *a: ('runtime', Path('image'), {'amalgkit': {'version': 'fixed'}}))
    commands = []
    def run(runtime, image, root, arguments, binds=(), **kwargs):
        commands.append(arguments)
        out = Path(arguments[arguments.index('--out_dir') + 1]) / 'metadata'
        out.mkdir()
        (out / 'metadata.tsv').write_bytes(original)
        return ['runtime', *arguments]
    monkeypatch.setattr(catalog, 'execute', run)
    result = catalog.metadata_stage(root, cfg, 'curate', 'work/curation', source)
    assert commands[0][:2] == ['amalgkit', 'select']
    assert '--sample_group' not in commands[0]
    assert source.read_bytes() == original
    assert Path(result['metadata']).is_file()
    assert json.loads((root / 'work/curation/curate/provenance.json').read_text())['status'] == 'complete'
    assert not (root / 'results').exists()


def test_failed_acquisition_records_attempt_and_requires_a_new_directory(curated_project, monkeypatch):
    root, dataset, cfg = curated_project
    monkeypatch.setattr(catalog, 'software_identity', lambda *a: ('runtime', Path('image'), {'amalgkit': {'version': 'fixed'}}))
    def fail(*args, **kwargs):
        raise subprocess.CalledProcessError(1, ['amalgkit', 'metadata'])
    monkeypatch.setattr(catalog, 'execute', fail)
    with pytest.raises(ValueError, match='AMALGKIT fetch failed'):
        catalog.metadata_stage(root, cfg, 'fetch', 'work/failed')
    provenance = json.loads((root / 'work/failed/fetch/provenance.json').read_text())
    assert provenance['status'] == 'failed'
    assert provenance['exit_code'] == 1
    assert not (root / 'results').exists()
    with pytest.raises(ValueError, match='stage already exists'):
        catalog.metadata_stage(root, cfg, 'fetch', 'work/failed')


def test_historical_name_cleanup_retains_original_name():
    r = row(1, 'A1'); r['scientific_name'] = '[Elymus] breviaristatus'
    normalized = catalog.normalized_row(r)
    assert normalized['scientific_name'] == 'Elymus breviaristatus'
    assert normalized['scientific_name_original'] == '[Elymus] breviaristatus'
    assert normalized['analysis_sample_id'] == 'Elymus_breviaristatus_A1'


def test_historical_seed_carries_exclusions_without_retaining_old_success(curated_project):
    root, dataset, cfg = curated_project
    source = root / 'historical.tsv'
    write_tsv(source, FIELDS, [row(1, 'OLD1'), row(2, 'BAD2')])
    exclusions = root / 'historical_exclusions.tsv'
    write_tsv(exclusions, ['accession', 'reason'], [dict(accession='BAD2', reason='poor_quality')])
    with pytest.raises(ValueError, match='already exists'):
        catalog.seed_history(root, cfg, source, exclusions)
    result = catalog.seed_history(root, cfg, source, exclusions, 'native_build', replace=True)
    assert result == {'historical_samples': 2, 'accepted_samples': 0, 'exclusions': 1}
    assert not read_tsv(dataset / 'accepted_samples.tsv')
    assert read_tsv(dataset / 'excluded_accessions.tsv') == [
        dict(accession='BAD2', taxid='2', bioproject='P1', reason='poor_quality', source_build='native_build')]
    _, output = select(curated_project, [row(1, 'OLD1'), row(1, 'NEW1', bases=5000),
                                         row(2, 'RELATED2'), row(2, 'GOOD2', project='P2')])
    assert {r['run'] for r in read_tsv(output / 'metadata.tsv')} == {'NEW1', 'GOOD2'}


def test_fresh_initialization_preserves_only_enriched_exclusions(curated_project):
    root, dataset, cfg = curated_project
    exclusions = dataset / 'excluded_accessions.tsv'
    before = exclusions.read_bytes()
    with pytest.raises(ValueError, match='dataset already exists'):
        catalog.initialize_empty(root, cfg, exclusions)
    assert catalog.initialize_empty(root, cfg, exclusions, replace=True) == {
        'samples': 0, 'accepted_samples': 0, 'exclusions': 1}
    assert not read_tsv(dataset / 'metadata.tsv')
    assert not read_tsv(dataset / 'accepted_samples.tsv')
    assert exclusions.read_bytes() == before
    assert json.loads((dataset / 'provenance.json').read_text())['kind'] == 'awaiting_fresh_query'
    summary, output = select(curated_project, [row(1, 'OLD1'), row(1, 'BIG1', bases=5000),
                                              row(2, 'BAD2'), row(2, 'GOOD2', project='P2')])
    assert summary['decisions'] == {'new': 2}
    assert {r['run'] for r in read_tsv(output / 'metadata.tsv')} == {'BIG1', 'GOOD2'}


def test_fresh_fetch_uses_ncbi_query_without_historical_input(fresh_project, monkeypatch):
    root, dataset, cfg = fresh_project
    monkeypatch.setattr(catalog, 'software_identity', lambda *a: ('runtime', Path('image'), {'amalgkit': {'version': 'fixed'}}))
    commands = []
    def run(runtime, image, root, arguments, binds=(), **kwargs):
        commands.append(arguments)
        target = Path(arguments[arguments.index('--out_dir') + 1]) / 'metadata/metadata.tsv'
        write_tsv(target, FIELDS, [row(3, 'FRESH3')])
        return ['runtime', *arguments]
    monkeypatch.setattr(catalog, 'execute', run)
    result = catalog.metadata_stage(root, cfg, 'fetch', 'work/new_query')
    assert commands[0][:2] == ['amalgkit', 'metadata']
    assert commands[0][commands[0].index('--search_string') + 1] == cfg['search_string']
    assert '--metadata' not in commands[0]
    assert read_tsv(result['metadata'])[0]['run'] == 'FRESH3'
    assert not (dataset / 'metadata.tsv').exists()
    assert not (dataset / 'accepted_samples.tsv').exists()
    assert not (root / 'results').exists()


@pytest.fixture
def completed_baseline(curated_project, monkeypatch):
    """Supply completed product evidence; portable manifest validation has its own tests."""
    import build_products
    import portable_build
    root, dataset, cfg = curated_project
    database = root / 'database'
    (database / 'provenance').mkdir(parents=True)
    rows = [row(1, 'OLD1'), row(2, 'BAD2'), row(3, 'LOW3')]
    write_tsv(database / 'metadata.tsv', FIELDS, rows)
    write_tsv(root / 'original.tsv', FIELDS, rows)
    image = root / 'fixed.sif'
    image.write_bytes(b'fixed container test receipt')
    build = {'config': yaml.safe_load((dataset / 'build.yaml').read_text()),
             'software_lock': {'container': {'files': [file_record(image)]}}}
    write_json(database / 'provenance/build.json', build)
    manifest = database / 'manifest.json'
    write_json(manifest, {})
    products = {}
    for r in rows:
        complete = 40 if r['taxid'] == '3' else 80
        products[r['run']] = {'row': r, 'counts': {'busco_cds_single': complete,
            'busco_cds_duplicated': 0, 'busco_cds_fragmented': 0,
            'busco_cds_missing': 100 - complete, 'busco_cds_total': 100}}
    def load(*args, **kwargs):
        return {'bundle_root': str(database), 'build_id': 'fixed_build', 'products': products,
                'files': [file_record(database / p) for p in ('metadata.tsv', 'provenance/build.json')]}
    monkeypatch.setattr(build_products, 'load_complete', load)
    monkeypatch.setattr(portable_build, 'completion_path', lambda _: manifest)
    return curated_project, database, build


@pytest.mark.parametrize('problem', ['no_receipt', 'different_image', 'different_revision'])
def test_baseline_rejects_native_or_different_container(completed_baseline, problem):
    (root, dataset, cfg), database, build = completed_baseline
    if problem == 'no_receipt':
        build['software_lock'] = None
    elif problem == 'different_image':
        build['config']['genegalleon']['image_uri'] = 'docker://old/image@sha256:' + 'b' * 64
    else:
        build['config']['genegalleon']['revision'] = 'b' * 40
    write_json(database / 'provenance/build.json', build)
    before = (dataset / 'accepted_samples.tsv').read_bytes()
    with pytest.raises(ValueError, match='matching pinned-container receipt'):
        catalog.initialize(root, cfg, database, root / 'original.tsv',
                           dataset / 'excluded_accessions.tsv', replace=True)
    assert (dataset / 'accepted_samples.tsv').read_bytes() == before


@pytest.mark.parametrize('threshold, expected', [(0.5, ['OLD1']), (0.4, ['OLD1', 'LOW3']), (0.4001, ['OLD1'])])
def test_baseline_applies_adoption_busco_threshold(completed_baseline, threshold, expected):
    (root, dataset, cfg), database, build = completed_baseline
    definition = yaml.safe_load(cfg['_path'].read_text())
    definition['busco_threshold'] = threshold
    cfg['_path'].write_text(yaml.safe_dump(definition))
    cfg = catalog.configuration(root, cfg['_path'])
    result = catalog.initialize(root, cfg, database, root / 'original.tsv',
                                dataset / 'excluded_accessions.tsv', replace=True)
    assert result == {'accepted_samples': len(expected), 'exclusions': 1}
    accepted = read_tsv(dataset / 'accepted_samples.tsv')
    assert [(r['run'], r['source_build']) for r in accepted] == [(run, 'fixed_build') for run in expected]
    assert [r['run'] for r in read_tsv(dataset / 'metadata.tsv')] == expected


@pytest.mark.parametrize('threshold', [-0.1, 1.1, None, True, '0.5', float('nan'), float('inf')])
def test_invalid_adoption_busco_threshold_is_rejected(curated_project, threshold):
    root, dataset, cfg = curated_project
    definition = yaml.safe_load(cfg['_path'].read_text())
    definition['busco_threshold'] = threshold
    cfg['_path'].write_text(yaml.safe_dump(definition))
    with pytest.raises(ValueError, match='busco_threshold must be'):
        catalog.configuration(root, cfg['_path'])


def test_adoption_busco_threshold_defaults_to_half(curated_project):
    root, dataset, cfg = curated_project
    definition = yaml.safe_load(cfg['_path'].read_text())
    definition.pop('busco_threshold')
    cfg['_path'].write_text(yaml.safe_dump(definition))
    assert catalog.configuration(root, cfg['_path'])['busco_threshold'] == 0.5


@pytest.mark.parametrize('mode', ['omitted', 'null', 'missing_file'])
def test_overrides_are_optional(curated_project, mode):
    root, dataset, cfg = curated_project
    definition = yaml.safe_load(cfg['_path'].read_text())
    if mode == 'omitted': definition.pop('overrides')
    elif mode == 'null': definition['overrides'] = None
    else: (dataset / 'overrides.tsv').unlink()
    cfg['_path'].write_text(yaml.safe_dump(definition))
    cfg = catalog.configuration(root, cfg['_path'])
    _, output = select((root, dataset, cfg), [row(1, 'OLD1')])
    catalog.accept_candidate(root, cfg, output)
    assert read_tsv(dataset / 'metadata.tsv')[0]['run'] == 'OLD1'


def test_single_update_previews_then_publishes_without_success_registration(curated_project):
    root, dataset, cfg = curated_project
    source = root / 'curated.tsv'
    write_tsv(source, FIELDS, [row(1, 'OLD1'), row(2, 'REPLACEMENT2', bases=500)])
    before = (dataset / 'metadata.tsv').read_bytes()
    preview = catalog.update_metadata(root, cfg, 'work/preview', source, dry_run=True)
    assert (dataset / 'metadata.tsv').read_bytes() == before
    assert preview['summary']['samples'] == 2
    catalog.update_metadata(root, cfg, 'work/update', source)
    assert {r['run'] for r in read_tsv(dataset / 'metadata.tsv')} == {'OLD1', 'REPLACEMENT2'}
    assert [r['run'] for r in read_tsv(dataset / 'accepted_samples.tsv')] == ['OLD1']


def test_excluding_last_representative_publishes_empty_metadata(curated_project):
    root, dataset, cfg = curated_project
    write_tsv(dataset / 'excluded_accessions.tsv', ['accession', 'reason'],
              [{'accession': 'OLD1', 'reason': 'misidentified'}, {'accession': 'BAD2', 'reason': 'unusable'}])
    source = root / 'curated.tsv'
    write_tsv(source, FIELDS, [row(1, 'OLD1'), row(2, 'BAD2')])
    result = catalog.update_metadata(root, cfg, 'work/empty', source)
    assert result['summary']['samples'] == 0
    assert read_tsv(dataset / 'metadata.tsv') == []
    assert read_tsv(dataset / 'accepted_samples.tsv') == []


def test_update_runs_fetch_and_curation_before_publishing(fresh_project, monkeypatch):
    root, dataset, cfg = fresh_project
    actions = []
    def stage(root, cfg, action, work):
        actions.append(action)
        target = Path(work) / action / 'metadata/metadata.tsv'
        write_tsv(target, FIELDS, [row(3, 'NEW3')])
        return {'metadata': str(target)}
    monkeypatch.setattr(catalog, 'metadata_stage', stage)
    catalog.update_metadata(root, cfg, 'work/full_update')
    assert actions == ['fetch', 'curate']
    assert [r['run'] for r in read_tsv(dataset / 'metadata.tsv')] == ['NEW3']
    assert not read_tsv(dataset / 'accepted_samples.tsv')
