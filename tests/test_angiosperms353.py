"""Validate rank voting, reference masking and the DB review boundary."""
from contextlib import contextmanager
from decimal import Decimal
import io
import sqlite3
import sys
from types import SimpleNamespace
import xml.etree.ElementTree as ET

import pytest
import yaml

import angiosperms353 as qc
import angiosperms353_reference as reference
from common import read_tsv, write_json, write_tsv
from dataset_assets import record
import metadata_catalog as catalog
import taxonomy_screen_jobs as jobs
from test_sample_curation import finish, sample_project
from test_datasets import dataset_project


def test_kew_download_and_names(tmp_path, monkeypatch):
    class Response(io.BytesIO):
        headers = {'Content-Length': '3'}

    monkeypatch.setattr(reference, 'urlopen', lambda *_, **__: Response(b'abc'))
    downloaded = reference.cached_download('https://example.invalid/x', tmp_path / 'marker.gz', compress=True)
    with reference.text_open(downloaded) as handle:
        assert handle.read() == 'abc'
    monkeypatch.setattr(reference, 'urlopen', lambda *_, **__: Response(b'x'))
    with pytest.raises(ValueError, match='incomplete'):
        reference.cached_download('https://example.invalid/x', tmp_path / 'partial.gz', compress=True)
    assert not (tmp_path / 'partial.gz').exists()
    tree = tmp_path / 'tree.nwk'
    tree.write_text('(Order_Family_Genus_species_ERR1_ERR2,Order_Family_Genus_other_ERR3)1;')
    assert reference.kew_families(tree) == ({'Genus_species_ERR1_ERR2': 'Family', 'Genus_other_ERR3': 'Family'},
                                            {'Genus': 'Family'})
    assert reference.accessions('ERR1_ERR2 SRR3 unrelated') == {'ERR1', 'ERR2', 'SRR3'}


def test_taxonomy_reconciles_synonyms_and_rejects_family_conflict(tmp_path):
    db = tmp_path / 'tax.sqlite'
    with sqlite3.connect(db) as connection:
        connection.executescript('''
            CREATE TABLE species(taxid INT, spname TEXT COLLATE NOCASE, rank TEXT, track TEXT);
            CREATE TABLE synonym(taxid INT, spname TEXT COLLATE NOCASE);
            CREATE TABLE merged(taxid_old INT, taxid_new INT);
            INSERT INTO species VALUES (1,'root','no rank','1'), (2,'Family','family','2,7,1'),
              (3,'Otheraceae','family','3,1'), (4,'Genus','genus','4,6,5,2,7,1'),
              (5,'Subfamily','subfamily','5,2,7,1'), (6,'Tribe','tribe','6,5,2,7,1'),
              (7,'Order','order','7,1'), (8,'Genus species','species','8,4,6,5,2,7,1'),
              (9,'Genus second','species','9,4,6,5,2,7,1'), (10,'Othergenus','genus','10,2,7,1');
            INSERT INTO synonym VALUES (4,'Synonymus'), (2,'Oldaceae'), (11,'Othergenus species'),
              (8,'Ambiguous species'), (9,'Ambiguous species');
            INSERT INTO merged VALUES (11,8);
        ''')
    tax = reference.ReferenceTaxonomy(db)
    try:
        assert tax.reconcile('synonymus', 'Oldaceae') == (2, 4)
        assert tax.reconcile('Genus', 'Otheraceae') == (None, None)
        assert tax.reconcile('Unknown', 'Family') == (2, None)
        assert dict(zip(qc.RANKS, tax.reference_ranks('Synonymus', 'Oldaceae')))['tribe'] == 6
        unknown = dict(zip(qc.RANKS, tax.reference_ranks('Unknown', 'Family')))
        assert unknown['family'] == 2 and unknown['order'] == 7 and unknown['genus'] is None
        assignment = tax.resolve_reference('Othergenus_species', 'Family')
        assert assignment['method'] == 'species' and assignment['taxid'] == 8
        assert assignment['ranks'][qc.RANKS.index('genus')] == 4
        assert tax.resolve_reference('Othergenus_species', 'Otheraceae')['family_conflict']
        ambiguous = tax.resolve_reference('Ambiguous_species', 'Family')
        assert ambiguous['candidates'] == (8, 9)
        assert ambiguous['ranks'][0] == 2 and ambiguous['ranks'][qc.RANKS.index('genus')] is None
        assert tax.resolve_reference('Othergenus_unlisted', 'Family')['method'] == 'genus'
    finally:
        tax.db.close()


def test_common_identity_breaks_mmseqs_ties_and_splits_equal_alternatives():
    refs = [('One', 'INSDC', 'ERR1', (1, 10)), ('Two', 'INSDC', 'ERR2', (2, 20)),
            ('Three', 'INSDC', 'ERR3', (2, 30))]
    query = 'A' * 200
    original = qc.Hit(0, 390, 1, 200, query, 'CC' + query[2:])
    alternative = qc.Hit(1, 400, 1, 200, query, query)
    lower_identity = qc.Hit(2, 400, 1, 200, query, 'C' + query[1:])
    equal_identity = qc.Hit(2, 400, 200, 1, query, query)

    family = qc.vote('A', 'family', '1', 200, 1, {1: original, 2: alternative}, refs)
    genus = qc.vote('A', 'genus', '1', 200, 10, {10: original, 20: lower_identity, 30: alternative}, refs)
    split = qc.vote('A', 'genus', '1', 200, 10, {10: original, 20: alternative, 30: equal_identity}, refs)
    tie = qc.vote('A', 'genus', '1', 200, 10, {10: alternative, 20: equal_identity}, refs)
    assert (family['vote'], family['registered_matches'], family['alternative_matches']) == ('alternative', 198, 200)
    assert genus['winner'] == 30 and genus['alternative_winners'] == [30]
    assert split['winner'] is None and split['alternative_winners'] == [20, 30]
    assert tie['vote'] == 'tie' and tie['alternative_winners'] == [20]
    assert qc.vote('A', 'family', '1', 200, 9, {2: alternative}, refs)['vote'] == 'no_registered_hit'

    # Both alternatives individually overlap 170 bases, but share only 140
    # with each other and the registered hit.
    left = qc.Hit(1, 400, 1, 170, query[:170], query[:170])
    right = qc.Hit(2, 400, 31, 200, query[:170], query[:170])
    assert qc.vote('A', 'genus', '1', 200, 10, {10: original, 20: left, 30: right}, refs)['vote'] == 'insufficient_overlap'


@pytest.mark.parametrize('rank', ['genus', 'tribe', 'family', 'order'])
def test_same_strict_registered_threshold_at_all_ranks(rank):
    sample = dict(run='A', taxid='10')
    split = dict(vote='alternative', alternative_winners=[20, 30])
    single = dict(vote='alternative', alternative_winners=[20])
    registered = dict(vote='registered')
    tie = dict(vote='tie', alternative_winners=[20, 30])
    cutoff = qc.registered_threshold(50)
    assert cutoff == Decimal('.10')
    assert qc.registered_threshold(50, support='.75') == Decimal('.25')

    rows = [split] * 31 + [single] * 9 + [tie] * 9 + [registered]
    result = qc.summarize(sample, rank, 10, rows, 50, 50, cutoff)
    assert result['status'] == 'discordant' and result['tied_loci'] == 9
    assert result['comparable_loci'] == 50 and result['voted_loci'] == 50
    assert result['alternative'] == 20 and result['alternative_votes'] == 27.5
    assert result['registered_votes'] == 4 and result['alternative_votes_total'] == 46
    assert Decimal(result['registered_support']) == Decimal('.08')
    assert Decimal(result['discordance_support']) == Decimal('.92')
    boundary = qc.summarize(sample, rank, 10, [split] * 45 + [registered] * 5 + [dict(vote='no_registered_hit')],
                            51, 50, cutoff)
    assert boundary['status'] == 'no_discordance'
    assert boundary['comparable_loci'] == 50
    assert Decimal(boundary['discordance_support']) == Decimal('.90')
    assert qc.summarize(sample, rank, 10, [tie] * 49, 49, 50,
                        cutoff)['status'] == 'insufficient_comparison'


def test_closest_reference_respects_masking_and_can_reach_order():
    def ref(source, **lineage):
        return ('Species_x', 'INSDC', source, tuple(lineage.get(r) for r in qc.RANKS), {source}, 'A')

    sample = dict(run='ERR1', taxid='10', source_runs={'ERR1'}, exclude_accessions=set())
    lineage = {'genus': 10, 'tribe': 20, 'family': 30, 'order': 40}
    masked = ref('ERR1', genus=10, tribe=20, family=30, order=40)
    order = ref('ERR2', genus=11, tribe=21, family=31, order=40)
    tribe = ref('ERR3', genus=12, tribe=20, family=30, order=40)
    assert qc.closest_reference_rank(lineage, [masked, order], sample, False) == 'order'
    assert qc.closest_reference_rank(lineage, [masked, order, tribe], sample, False) == 'tribe'
    assert qc.closest_reference_rank(lineage, [masked], sample, False) is None


def test_later_genus_reference_prevents_fallback_on_insufficient_overlap(tmp_path, monkeypatch):
    sample = dict(run='ERR1', taxid='10', source_runs=set(), exclude_accessions=set())
    lineage = {'genus': 10, 'family': 20}
    marker_path = tmp_path / 'markers.fna'
    marker_path.write_text('>sample__1\n' + 'A' * 200 + '\n>sample__2\n' + 'A' * 200 + '\n')
    @contextmanager
    def search(locus, *args):
        genus = 10 if locus == '2' else 11
        ranks = tuple({'family': 20, 'genus': genus}.get(r) for r in qc.RANKS)
        yield [('Species_x', 'INSDC', 'ERR2', ranks, set(), 'A')], marker_path

    def hits(*args, **kwargs):
        return {0: {'family': {20: qc.Hit(0, 400, 1, 200, 'A' * 200, 'A' * 200),
                               21: qc.Hit(0, 300, 1, 200, 'A' * 200, 'C' * 200)}}}

    monkeypatch.setattr(qc, 'search_locus', search)
    monkeypatch.setattr(qc, 'read_hits', hits)
    evidence = SimpleNamespace(rows=[])
    evidence.writerows = evidence.rows.extend
    markers = qc.markers(marker_path)
    summaries = qc.score_batch([sample], [lineage], {locus: {0: sequence} for locus, sequence in markers.items()},
                               dict.fromkeys(markers, marker_path), None, tmp_path, tmp_path,
                               1, 2, False, 1, Decimal('.20'), evidence)
    assert len(summaries) == 1
    assert summaries[0]['rank'] == 'genus'
    assert summaries[0]['status'] == 'insufficient_comparison'
    assert {r['rank'] for r in evidence.rows} == {'genus'}


def test_cigar_reconstructs_gaps_and_reverse_strands():
    assert qc.alignment_from_cigar('AACCG','AACGG',1,5,1,5,'2M1I1M1D1M') == ('AACC-G','AA-CGG')
    assert qc.alignment_from_cigar('AACG','AACG',4,1,4,1,'0M4M') == ('CGTT','CGTT')
    with pytest.raises(ValueError, match='disagrees'):
        qc.alignment_from_cigar('AACG','AACG',1,4,1,4,'3M')


def test_mask_same_source_before_picking_best_hit(tmp_path):
    seq = 'A' * 160
    refs = [('One_species', 'INSDC', 'ERR1_ERR2', tuple({'family': 1, 'genus': 10}.get(r) for r in qc.RANKS), {'ERR1', 'ERR2'}, seq),
            ('One_other', 'INSDC', 'ERR3', tuple({'family': 1, 'genus': 10}.get(r) for r in qc.RANKS), {'ERR3'}, seq)]
    hits = tmp_path / 'hits.tsv'
    hits.write_text('0\t0\t320\t1\t160\t1\t160\t160M\n0\t1\t315\t1\t160\t1\t160\t160M\n')
    sample = dict(source_runs={'ERR2'}, exclude_accessions=set(), species_name='one species')
    assert [hit.target for hit in qc.read_hits(hits, {0: seq}, refs, [sample], omit_same_species=False)[0]['family'].values()] == [1]


def test_screen_config_skips_work_and_replots_exact_cutoff(sample_project):
    root, cfg, _ = sample_project
    assert cfg['taxonomy_screen'] == {'enabled': False, 'method': 'angiosperms353', 'max_registered_support': 0.1, 'batch_size': 24}
    src = root / 'source_screen'
    src.mkdir()
    args = SimpleNamespace(command='report-taxonomy', screen=str(src / 'screen.json'),
                           output=str(root / 'report'))
    assert catalog.taxonomy_screen_command(root, cfg, args) == {'status': 'disabled'}
    assert not (root / 'report').exists()

    build_path = cfg['_path']
    definition = yaml.safe_load(build_path.read_text())
    definition['taxonomy_screen'] = {'enabled': True, 'method': 'angiosperms353', 'max_registered_support': 0.2}
    build_path.write_text(yaml.safe_dump(definition))
    cfg = catalog.configuration(root, build_path)
    assert cfg['previous_metadata'] == definition['metadata']
    assert cfg['excluded_accessions'] == definition['excluded_accessions']
    summaries = []
    evidence = []
    for run, rank, registered_loci in [('A1', 'genus', 2), ('B1', 'family', 1), ('C1', 'family', 0)]:
        rows = [dict(run=run, rank=rank, locus=locus, vote='registered') if locus < registered_loci else
                dict(run=run, rank=rank, locus=locus, vote='alternative', alternative_winners=[20])
                for locus in range(10 if run != 'C1' else 1)]
        evidence.extend(rows)
        summaries.append(qc.summarize(dict(run=run, taxid='10'), rank, 10, rows, len(rows), 10, Decimal('.3')))
    write_tsv(src / 'summary.tsv', qc.SUMMARY_FIELDS, summaries)
    write_tsv(src / 'loci.tsv', qc.LOCUS_FIELDS, evidence)
    write_json(src / 'screen.json', dict(kind='angiosperms353_qc', schema_version=3,
               settings=dict(min_loci=10, max_registered_support='0.30'),
               summary=record(src / 'summary.tsv'), loci=record(src / 'loci.tsv')))
    result = catalog.taxonomy_screen_command(root, cfg, args)
    assert result['statuses'] == {'no_discordance': 1, 'discordant': 1, 'insufficient_comparison': 1}
    assert [row['accession'] for row in read_tsv(root / 'report/exclusion_candidates.tsv')] == ['B1']
    receipt = qc.json.loads((root / 'report/screen.json').read_text())
    assert receipt['settings']['max_registered_support'] == '0.2'
    assert receipt['source_screen']['path'] == str(src / 'screen.json')
    for name, titles in zip(qc.PLOT_NAMES, [('All ranks',), ('family', 'genus')], strict=True):
        plot = (root / 'report' / name).read_text()
        ET.fromstring(plot)
        assert 'Threshold: 20%' in plot
        assert all(title + ' · ' in plot for title in titles)
    assert read_tsv(src / 'summary.tsv')[0]['status'] == 'discordant'


def test_screen_portable_database_and_reuse_retained_products(sample_project, monkeypatch):
    import dataset
    from database_fixtures import database_from_stages

    root, _, build = sample_project
    for index in range(1, 4):
        finish(build, index)
    prepared = dataset.load(build)
    source = database_from_stages(root, prepared['config']['store'], prepared['items'])
    rows, manifest = qc.database_panel(source, ['B1'])
    assert [row['run'] for row in rows] == ['B1']
    assert qc.verify(rows[0]['cds']).is_relative_to(source)
    before = record(manifest['path'])
    config = root / 'config/build.yaml'
    selection = yaml.safe_load(config.read_text())
    selection['taxonomy_screen'] = {'enabled': True, 'method': 'angiosperms353', 'max_registered_support': 0.15}
    config.write_text(yaml.safe_dump(selection))
    taxonomy = root / 'taxonomy.sqlite'
    taxonomy.touch()
    monkeypatch.setattr(qc, 'screen', lambda samples, **kwargs: dict(samples=samples, **kwargs))
    args = SimpleNamespace(command='screen-taxonomy', database=str(source), runs=['B1'], taxonomy=str(taxonomy),
                           output='screen', cache='cache', loci=None, omit_same_species=False,
                           threads=1, ram=4, min_loci=50, blat=None)
    call = catalog.taxonomy_screen_command(root, catalog.configuration(root, config, require_selection=False), args)
    assert call['samples'] == rows and call['database'] == manifest and call['build'] is None
    assert call['registered_support'] == 0.15

    exclusions = root / 'reviewed_exclusions.tsv'
    write_tsv(exclusions, ['accession', 'reason'], [dict(accession='A1', reason='reviewed taxonomy')])
    settings = yaml.safe_load(config.read_text())
    settings.update(metadata=str(source / 'metadata.tsv'), metadata_provenance=None,
                    reuse_from=str(source), excluded_accessions=str(exclusions))
    config.write_text(yaml.safe_dump(settings))
    report = dataset.plan(root, config)[-1]
    retained = [row for row in report if row['assembly'] != 'excluded']
    assert len(retained) == 2
    assert all(row['mapping'] == 'reuse' for row in retained)
    filtered = dataset.prepare(root, 'filtered', config)
    assert dataset.submit(filtered, until='quant', dry_run=True) == []
    assert {row['run'] for row in read_tsv(dataset.materialize(filtered) / 'metadata.tsv')} == {'B1', 'G1'}
    assert record(manifest['path']) == before

    from build_products import load_complete
    from mapping_tables import load_tables
    from taxonomy_screen_filter import filter_database
    rows, manifest = qc.database_panel(source)
    evidence = root / 'evidence'
    summaries = [dict(run=row['run'], status='discordant' if row['run'] == 'A1' else 'no_discordance')
                 for row in rows]
    write_tsv(evidence / 'summary.tsv', ['run', 'status'], summaries)
    write_tsv(evidence / 'exclusion_candidates.tsv', ['accession', 'reason'],
              [dict(accession='A1', reason='taxonomy evidence')])
    screen = evidence / 'screen.json'
    data = dict(kind='angiosperms353_qc', schema_version=3, database_manifest=manifest,
                settings=dict(restricted_loci=False, omit_same_species=False, min_loci=50,
                              max_registered_support='0.1'),
                summary=record(evidence / 'summary.tsv'), exclusion_candidates=record(evidence / 'exclusion_candidates.tsv'),
                samples={row['run']: dict(taxid=str(row['taxid']), cds=row['cds']) for row in rows})
    write_json(screen, data)
    output = root / 'filtered_database'
    result = filter_database(source, screen, output)
    assert (result['retained'], result['excluded']) == (2, 1)
    derived = load_complete(output)
    assert {p['row']['run'] for p in derived['products'].values()} == {'B1', 'G1'}
    assert set(load_tables(derived['mapping']['path'])['tables']) == set(derived['products'])
    assert read_tsv(output / 'excluded_accessions.tsv') == [dict(accession='A1', reason='taxonomy evidence')]
    assert filter_database(source, screen, output) == result
    assert record(manifest['path']) == before
    data['samples'].pop('A1')
    write_json(screen, data)
    with pytest.raises(ValueError, match='every sample'):
        filter_database(source, screen, root / 'incomplete_filter')


@pytest.mark.parametrize('version,excluded_rank', [(1, 'family'), (2, 'genus'), (3, 'auto')])
def test_record_uses_bound_screen_and_only_excludes_on_review_apply(sample_project, version, excluded_rank):
    root, cfg, build = sample_project
    finish(build, 1)
    rows, build_record = qc.build_panel(build, ['A1'])
    row = rows[0]
    summary = root / 'summary.tsv'
    ranks = {1: qc.LEGACY_RANKS, 2: qc.INTERMEDIATE_RANKS, 3: ('order',)}[version]
    write_tsv(summary, qc.SUMMARY_FIELDS,
              [dict(run='A1', taxid=row['taxid'], rank=rank, status='discordant', registered='10', alternative='20',
                    alternative_loci=40 if version >= 2 else '', voted_loci=50 if version == 3 else '')
               for rank in ranks])
    screen = root / 'screen.json'
    write_json(screen, dict(schema_version=version, kind='angiosperms353_qc', build_manifest=build_record,
                            settings=dict(restricted_loci=False, omit_same_species=False),
                            samples={'A1': {'taxid': row['taxid'], 'cds': row['cds']}}, summary=record(summary)))
    plan = root / 'review.tsv'
    observed = catalog.record_successes(root, cfg, build, ['A1'], dry_run=True, taxonomy_screen=screen)
    assert observed['recorded_runs'] == ['A1'] and observed['excluded_runs'] == []
    with pytest.raises(ValueError, match='require --taxonomy-screen and --plan or --dry-run'):
        catalog.record_successes(root, cfg, build, ['A1'], taxonomy_screen=screen, exclude_taxonomy=excluded_rank)
    proposed = catalog.record_successes(root, cfg, build, ['A1'], plan=plan,
                                        taxonomy_screen=screen, exclude_taxonomy=excluded_rank)
    assert proposed['excluded_runs'] == ['A1']
    assert read_tsv(root / cfg['excluded_accessions']) == []
    from sample_review import apply_plan
    assert apply_plan(root, cfg, plan)['excluded_runs'] == ['A1']
    reason = {1: 'angiosperms353_family_discordant_10_to_20',
              2: 'angiosperms353_genus_discordant_10_outside_registered_genus',
              3: 'angiosperms353_order_discordant_10_outside_registered_order'}[version]
    assert read_tsv(root / cfg['excluded_accessions'])[0]['reason'] == reason


def test_screen_batches_of_50_reuse_markers_after_failed_search(tmp_path, monkeypatch):
    cds = tmp_path / 'cds.fna'
    cds.write_text('>sample__1\n' + 'A' * 200 + '\n')
    tree = tmp_path / 'kew.nwk'
    tree.write_text('(Order_Family_Genus_species_ERR999);')
    reference_fasta = tmp_path / 'kew.fna'
    reference_fasta.write_text('>1 Species:Genus_species Repository:INSDC Sequence_ID:ERR999\nAAAA\n')
    taxonomy = tmp_path / 'taxonomy.sqlite'
    with sqlite3.connect(taxonomy) as connection:
        connection.executescript("""
            CREATE TABLE species(taxid INT, spname TEXT, rank TEXT, track TEXT);
            CREATE TABLE merged(taxid_old INT, taxid_new INT);
            CREATE TABLE synonym(taxid INT, spname TEXT COLLATE NOCASE);
            INSERT INTO species VALUES (10,'Genus','genus','10');
        """)
    monkeypatch.setattr(jobs, 'tool_versions', lambda: {'captus': 'test', 'mmseqs': 'test'})
    monkeypatch.setattr(jobs.KewReference, 'manifests', lambda self: [tree])
    monkeypatch.setattr(jobs.KewReference, 'loci', lambda self, selected: {'1': reference_fasta})
    extractions, scores = [], []

    def extract(samples, *_args, **_kwargs):
        extractions.append(len(samples))
        return [cds] * len(samples)

    def score(batch, ancestry, queries, paths, assignments, scratch, work, threads, ram, omit, minimum, cutoff,
              evidence):
        scores.append(len(batch))
        assert len(queries['1']) == len(batch)
        if len(scores) == 1:
            raise ValueError('simulated search failure')
        evidence.writerows(dict(run=row['run'], locus='1') for row in batch)
        return [qc.summarize(row, None, None, [], 1, minimum, cutoff) for row in batch]

    monkeypatch.setattr(qc, 'extract', extract)
    monkeypatch.setattr(qc, 'score_batch', score)
    samples = [dict(run=f'ERR{n}', taxid='10', cds=record(cds)) for n in range(101)]
    output = tmp_path / 'screen'
    kwargs = dict(cache=tmp_path, taxonomy=taxonomy, output=output, batch_size=50, blat=sys.executable,
                  config={'storage': {'scratch_root': 'workspace'}})
    with pytest.raises(ValueError, match='simulated search failure'):
        qc.screen(samples, **kwargs)
    assert (output / 'batches/0001/markers.json.gz').is_file()
    assert not (output / 'screen.json').exists()
    with pytest.raises(ValueError, match='incomplete taxonomy screen batches'):
        jobs.finalize(output, qc.json.loads((output / 'plan.json').read_text()))
    result = qc.screen(samples, **kwargs)
    assert result['batches'] == 3
    assert extractions == [50, 50, 1]
    assert scores == [50, 50, 50, 1]
    receipt = qc.json.loads((output / 'batches/0001/completed.json').read_text())
    assert receipt['reused_markers']
    assert receipt['recovery_seconds'] >= 0 and receipt['search_and_vote_seconds'] >= 0
    assert [r['run'] for r in read_tsv(output / 'loci.tsv')] == [r['run'] for r in samples]
    assert qc.screen(samples, **kwargs) == result
    assert extractions == [50, 50, 1]
    with pytest.raises(ValueError, match='different inputs'):
        qc.screen(samples, **dict(kwargs, batch_size=24))


@pytest.mark.parametrize('backend', ['slurm', 'age'])
def test_screen_scheduler_serializes_array_chunks_with_global_limits(tmp_path, monkeypatch, backend):
    monkeypatch.delenv('PHENORADAR_SCHEDULER', raising=False)
    execution = {'scheduler': backend, 'total_limits': {'cpus': 9, 'mem_gb': 40}, 'concurrency': 3}
    if backend == 'age':
        execution['age'] = {'parallel_environment': 'smp', 'memory_resource': 's_vmem'}
    monkeypatch.setattr(jobs.scheduler, 'array_size', lambda _: 2)
    monkeypatch.setattr(jobs.scheduler, 'active_ids', lambda _: set())
    commands = []

    def enqueue(command, **_kwargs):
        commands.append(command)
        return str(100 + len(commands))

    monkeypatch.setattr(jobs.subprocess, 'check_output', enqueue)
    result = jobs.submit(tmp_path, {'root': str(tmp_path)}, [2, 4, 7], execution,
                         {'cpus': 4, 'mem_gb': 16, 'time': '08:00:00'}, 8, {})
    assert result['arrays'] == ['101', '102']
    assert result['finalize_job'] == '103'
    assert result['concurrency'] == 2
    first, second, finalize = commands
    if backend == 'slurm':
        assert '--array=1,2%2' in first
        assert '--dependency=afterany:101' in second
        assert '--dependency=afterany:102' in finalize
    else:
        assert first[first.index('-t') + 1] == '1-2'
        assert first[first.index('-tc') + 1] == '2'
        assert second[second.index('-hold_jid') + 1] == '101'
        assert finalize[finalize.index('-hold_jid') + 1] == '102'
    mapping = qc.json.loads((tmp_path / 'jobs/0001_batch_001.json').read_text())
    assert mapping['indices'] == [2, 4]
    assert mapping['threads'] == 4
