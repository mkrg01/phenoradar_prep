"""Build a taxonomy reference from the OrthoDB proteins shipped by BUSCO."""
from collections import defaultdict
from contextlib import closing
from itertools import groupby
import json
from pathlib import Path
import re
import sqlite3
import subprocess
import tempfile

from angiosperms353_reference import RANKS, fasta
from busco_phylogeny import COORDINATES, busco_table
from common import atomic_writer, write_json
from dataset_assets import SAFE, digest, locked, record, verify
from taxonomy_check import Taxonomy

BUSCO_URL = 'https://busco-data.ezlab.org/v5/data/lineages/'
MARKER = re.compile(r'[0-9]+at[0-9]+')
GENE = re.compile(r'([1-9][0-9]*_[0-9]+):[0-9a-fA-F]+')
PROTEIN_ALPHABET = 'ACDEFGHIKLMNPQRSTVWY'
MIN_RESIDUES = 50


def lineage_directory(cache, lineage, supplied=None):
    if not isinstance(lineage, str) or not SAFE.fullmatch(lineage):
        raise ValueError('the OrthoDB screen requires a BUSCO lineage dataset name')
    if supplied is not None:
        return Path(supplied).resolve()
    cache = Path(cache).resolve() / 'busco'
    target = cache / 'lineages' / lineage
    with locked(cache / (lineage + '.lock')):
        if not target.is_dir():
            with tempfile.TemporaryDirectory(prefix='.busco-', dir=cache) as temporary:
                subprocess.run(['busco', '--download_path', temporary, '--download', lineage], check=True)
                downloaded = Path(temporary) / 'lineages' / lineage
                if not (downloaded / 'dataset.cfg').is_file():
                    raise ValueError(f'BUSCO did not download the requested dataset: {lineage}')
                target.parent.mkdir(parents=True, exist_ok=True)
                downloaded.rename(target)
    return target


def dataset_metadata(directory, lineage):
    metadata = dict(line.split('=', 1) for line in (directory / 'dataset.cfg').read_text().splitlines() if '=' in line)
    if metadata.get('name') != lineage or metadata.get('domain') != 'eukaryota':
        raise ValueError(f'expected a eukaryotic BUSCO dataset named {lineage}: {directory}')
    release = metadata.get('OrthoDB_version')
    if not release or not re.fullmatch(r'[0-9]+(?:\.[0-9]+)?', release):
        raise ValueError(f'BUSCO dataset does not specify an OrthoDB release: {directory}')
    markers = sorted({line.split()[0] for line in (directory / 'scores_cutoff').read_text().splitlines() if line.strip()})
    if (not markers or any(not MARKER.fullmatch(marker) for marker in markers)
            or str(len(markers)) != metadata.get('number_of_BUSCOs')):
        raise ValueError(f'invalid BUSCO marker catalogue: {directory}')
    return release, markers


def build_panel(source, markers, taxonomy_path, output):
    assignments = {}
    markers = set(markers)
    taxonomy = Taxonomy(taxonomy_path)

    def proteins():
        for header, sequence in fasta(source):
            marker, _, gene = header.split()[0].partition('_')
            match = GENE.fullmatch(gene)
            if marker not in markers or match is None or not sequence:
                raise ValueError(f'invalid reference protein header or sequence: {header}')
            organism = match[1]
            if organism not in assignments:
                lineage = taxonomy.lineage(int(organism.split('_')[0]))
                species = lineage.get('species')
                assignments[organism] = (dict(species=species[0], name=species[1],
                                             ranks=[lineage.get(rank, (None,))[0] for rank in RANKS])
                                         if species else None)
            if assignments[organism]:
                yield marker, organism, gene, sequence

    with closing(taxonomy.db), tempfile.TemporaryDirectory(prefix='orthodb-panel-') as temp:
        with closing(sqlite3.connect(Path(temp) / 'panel.sqlite')) as db:
            db.execute('PRAGMA journal_mode=OFF')
            db.execute('CREATE TABLE proteins (marker TEXT, organism TEXT, gene TEXT, sequence TEXT)')
            db.executemany('INSERT INTO proteins VALUES (?, ?, ?, ?)', proteins())
            db.execute('CREATE INDEX protein_groups ON proteins (marker, organism)')
            db.execute('''CREATE TABLE single_copy AS
                SELECT min(rowid) AS protein, organism FROM proteins GROUP BY marker, organism HAVING count(*) = 1''')
            representatives = defaultdict(list)
            for organism, count in db.execute('SELECT organism, count(*) FROM single_copy GROUP BY organism'):
                species = assignments[organism]['species']
                representatives[species].append((-count, int(organism.split('_')[0]) != species, organism))
            selected = sorted(min(candidates)[2] for candidates in representatives.values())
            db.execute('CREATE TABLE selected (organism TEXT PRIMARY KEY)')
            db.executemany('INSERT INTO selected VALUES (?)', ((name,) for name in selected))
            rows = db.execute('''SELECT p.marker, p.gene, p.sequence FROM single_copy s
                JOIN selected a ON a.organism = s.organism JOIN proteins p ON p.rowid = s.protein
                ORDER BY p.marker, p.gene''')
            loci = []
            for marker, hits in groupby(rows, key=lambda row: row[0]):
                with atomic_writer(output / f'{marker}.faa.gz') as handle:
                    for _, gene, sequence in hits:
                        handle.write(f'>{gene}\n{sequence}\n')
                loci.append(marker)
    if not loci:
        raise ValueError('no OrthoDB marker proteins resolved to species in the taxonomy snapshot')
    write_json(output / 'taxonomy.json', {name: assignments[name] for name in selected}, compact=True)
    return loci, len(selected)


def prepare(cache, lineage, taxonomy_path, supplied=None):
    directory = lineage_directory(cache, lineage, supplied)
    release, markers = dataset_metadata(directory, lineage)
    sources = [record(directory / name) for name in ('dataset.cfg', 'scores_cutoff', 'refseq_db.faa.gz')]
    key = digest([1, record(taxonomy_path)['sha256'], *(entry['sha256'] for entry in sources)])
    root = Path(cache).resolve() / 'orthodb' / lineage
    output = root / key
    with locked(root / (key + '.lock')):
        if not output.exists():
            with tempfile.TemporaryDirectory(prefix='.reference-', dir=root) as temp:
                staging = Path(temp) / 'reference'
                staging.mkdir()
                loci, species = build_panel(verify(sources[-1]), markers, taxonomy_path, staging)
                write_json(staging / 'panel.json', dict(loci=loci, species=species))
                staging.rename(output)
    panel = json.loads((output / 'panel.json').read_text())
    return dict(release=release, url=BUSCO_URL, busco_dataset=str(directory), markers=markers,
                species=panel['species'], manifests=[*sources, record(output / 'panel.json')],
                taxonomy_assignments=record(output / 'taxonomy.json'),
                loci={locus: record(output / f'{locus}.faa.gz') for locus in panel['loci']})


def reference_rows(path, locus, assignments):
    for header, sequence in fasta(path):
        organism = header.split(':', 1)[0]
        entry = assignments[organism]
        yield entry['name'], 'OrthoDB', header, entry['ranks'], sequence, entry['species'], organism


def protein(sequence, table):
    from Bio.Seq import Seq

    sequence = re.sub(r'[^ACGTURYSWKMBDHVN]', 'N', sequence[:len(sequence) // 3 * 3])
    return str(Seq(sequence).translate(table=table)).rstrip('*').replace('*', 'X')


def recover(sample, reference, lineage, table, work, log_dir, threads, only_loci):
    cds = verify(sample['cds'])
    markers = set(reference['markers'])
    complete = None
    if sample.get('busco') and sample.get('busco_lineage') == lineage:
        hits, observed = busco_table(verify(sample['busco']), lineage)
        if markers == observed:
            complete = hits
    if complete is None:
        proteins = work / f"{sample['run']}.faa"
        with proteins.open('w') as handle:
            for header, sequence in fasta(cds):
                translated = protein(sequence, table)
                if translated:
                    handle.write(f'>{header.split()[0]}\n{translated}\n')
        if not proteins.stat().st_size:
            return {}
        dataset = Path(reference['busco_dataset'])
        with (log_dir / 'busco.log').open('ab') as log:
            subprocess.run(['busco', '--offline', '--opt-out-run-stats', '-m', 'proteins',
                            '-l', str(dataset), '-i', str(proteins),
                            '-o', sample['run'], '--out_path', str(work), '-c', str(threads)],
                           check=True, cwd=work, stdout=log, stderr=subprocess.STDOUT)
        full = work / sample['run'] / f'run_{lineage}' / 'full_table.tsv'
        complete, observed = busco_table(full, lineage)
        if observed != markers:
            raise ValueError(f'{sample["run"]}: BUSCO results differ from the reference marker catalogue')
    wanted = reference['loci'].keys() if only_loci is None else only_loci
    targets = {COORDINATES.sub('', hit[1]): locus for locus, hit in complete.items() if locus in wanted}
    recovered, seen = {}, set()
    for header, sequence in fasta(cds):
        gene = header.split()[0]
        if gene in targets:
            if gene in seen:
                raise ValueError(f'{sample["run"]}: duplicate CDS identifier {gene}')
            seen.add(gene)
            translated = protein(sequence, table)
            if len(translated) >= MIN_RESIDUES:
                recovered[targets[gene]] = translated
    if targets.keys() - seen:
        raise ValueError(f'{sample["run"]}: BUSCO table refers to CDS identifiers absent from the assembly')
    return recovered
