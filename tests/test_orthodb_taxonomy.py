import gzip
import json
import sqlite3
import subprocess

from common import read_tsv
from dataset_assets import record
from angiosperms353 import screen
from angiosperms353_reference import fasta


def test_orthodb_screen_reuses_busco_and_selects_rank_with_enough_votes(tmp_path, monkeypatch):
    taxonomy = tmp_path / 'taxonomy.sqlite'
    with sqlite3.connect(taxonomy) as db:
        db.executescript("""
            CREATE TABLE species(taxid INT, spname TEXT, rank TEXT, track TEXT);
            CREATE TABLE merged(taxid_old INT, taxid_new INT);
            INSERT INTO species VALUES (30,'Order','order','30');
            INSERT INTO species VALUES (20,'Family A','family','20,30'), (21,'Family B','family','21,30');
            INSERT INTO species VALUES (10,'Genus A','genus','10,20,30'),
                (11,'Genus B','genus','11,20,30'), (12,'Genus C','genus','12,21,30');
            INSERT INTO species VALUES (101,'Genus A queried','species','101,10,20,30'),
                (102,'Genus B related','species','102,11,20,30'),
                (103,'Genus C alternative','species','103,12,21,30');
        """)

    lineage = 'vertebrata_odb12.2'
    dataset = tmp_path / lineage
    dataset.mkdir()
    (dataset / 'dataset.cfg').write_text(f'name={lineage}\ndomain=eukaryota\nOrthoDB_version=12.2\nnumber_of_BUSCOs=2\n')
    loci = ['1at7742', '2at7742']
    (dataset / 'scores_cutoff').write_text(''.join(f'{locus}\t20\n' for locus in loci))
    with gzip.open(dataset / 'refseq_db.faa.gz', 'wt') as handle:
        for locus in loci:
            for taxid, aa in [(101, 'F'), (102, 'F'), (103, 'M')]:
                handle.write(f'>{locus}_{taxid}_0:01\n{aa * 80}\n')
        handle.write(f'>{loci[1]}_101_0:02\n' + 'L' * 80 + '\n')
    cds = tmp_path / 'sample.fna'
    cds.write_text('>cds1\n' + 'ATG' * 80 + '\n>cds2\n' + 'ATG' * 80 + '\n')
    table = tmp_path / 'full_table.tsv'
    table.write_text(
        f'# The lineage dataset is: {lineage} (Creation date: 2026-05-13)\n'
        '# Busco id\tStatus\tSequence\tScore\tLength\n'
        + ''.join(f'{locus}\tComplete\tcds{index}:1-240\t80\t80\n'
                  for index, locus in enumerate(loci, start=1)))
    monkeypatch.setattr('taxonomy_screen_jobs.tool_versions', lambda method: {'busco': 'test', 'mmseqs': 'test'})
    run = subprocess.run

    def search(argv, **kwargs):
        assert argv[0] != 'busco', 'the compatible BUSCO table should be reused'
        if argv[0] != 'mmseqs':
            return run(argv, **kwargs)
        assert argv[:2] == ['mmseqs', 'easy-search'] and argv[argv.index('--search-type') + 1] == '1'
        [(query, sequence)] = list(fasta(argv[2]))
        with open(argv[4], 'w') as handle:
            for target, protein in fasta(argv[3]):
                score = sum(a == b for a, b in zip(sequence, protein, strict=True))
                handle.write(f'{query}\t{target}\t{score}\t1\t80\t1\t80\t80M\n')

    monkeypatch.setattr('angiosperms353.subprocess.run', search)
    samples = [dict(run='ERR1', taxid='101', cds=record(cds), busco=record(table), busco_lineage=lineage)]
    output = tmp_path / 'screen'
    options = dict(cache=tmp_path / 'cache', taxonomy=taxonomy, output=output, busco_lineage=lineage,
                   busco_dataset=dataset, min_loci=2, config={'storage': {'temporary_directory': 'workspace'}})
    result = screen(samples, **options)
    assert screen(samples, **options) == result
    assert result['statuses'] == {'discordant': 1}
    [summary] = read_tsv(output / 'summary.tsv')
    assert (summary['rank'], summary['alternative'], summary['voted_loci']) == ('family', '21', '2')
    [row] = read_tsv(output / 'taxonomy_summary.tsv')
    assert row['Assessment rank'] == 'family'
    assert not (output / 'taxonomy_by_rank.tsv').exists()
    assert read_tsv(output / 'loci.tsv')[0]['shared_residues'] == '80'
    manifest = json.loads((output / 'screen.json').read_text())
    assert manifest['kind'] == 'orthodb_qc' and manifest['reference']['release'] == '12.2'
    assert read_tsv(output / 'exclusion_candidates.tsv')[0]['reason'].startswith('orthodb_family_discordant')
