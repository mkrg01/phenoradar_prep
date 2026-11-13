"""Cache Kew Angiosperms353 markers and resolve reference taxonomy offline."""
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from contextlib import nullcontext
from functools import lru_cache
import fcntl
import gzip
import re
from pathlib import Path
from urllib.request import urlopen

from common import atomic_writer, read_tsv, text_open
from taxonomy_check import Taxonomy

RELEASE = '4.0'
BASE_URL = f'https://sftp.kew.org/pub/paftol/releases/{RELEASE}'
TREE = f'tree/species/treeoflife.{RELEASE}.tree'
LOCUS = re.compile(r'[0-9]+')
RUN = re.compile(r'(?:[DES]RR)[0-9]+', re.I)
LEGACY_RANKS = ('family', 'genus')
INTERMEDIATE_RANKS = ('family', 'subfamily', 'supertribe', 'tribe', 'subtribe', 'genus')
RANKS = (*INTERMEDIATE_RANKS, 'superfamily', 'parvorder', 'infraorder', 'suborder',
         'order', 'superorder', 'infraclass', 'subclass', 'class', 'superclass',
         'subphylum', 'phylum', 'superphylum', 'subkingdom', 'kingdom', 'superkingdom', 'domain')


def fasta(path):
    """Stream full headers, including Kew's source accession and species."""
    with text_open(path) as handle:
        header, sequence = None, []
        for raw in handle:
            line = raw.strip()
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(sequence).upper()
                header, sequence = line[1:], []
                if not header:
                    raise ValueError(f'empty FASTA header in {path}')
            elif line:
                if header is None:
                    raise ValueError(f'sequence before FASTA header in {path}')
                sequence.append(line)
        if header is not None:
            yield header, ''.join(sequence).upper()


def accessions(value):
    return {run.upper() for run in re.split(r'[,;_\s]+', value) if RUN.fullmatch(run)}


def cached_download(url, path, *, compress=False):
    path = Path(path)
    if path.is_file():
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.with_name('.' + path.name + '.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if path.is_file():
            return path
        with urlopen(url, timeout=120) as response, atomic_writer(path, 'wb') as output:
            expected = response.headers.get('Content-Length')
            size = 0
            stream = gzip.GzipFile(fileobj=output, mode='wb', filename='', mtime=0, compresslevel=3) if compress else nullcontext(output)
            with stream as destination:
                for block in iter(lambda: response.read(1024 * 1024), b''):
                    destination.write(block)
                    size += len(block)
            if not size or (expected is not None and size != int(expected)):
                raise ValueError(f'incomplete Kew download: {url}')
    return path


class KewReference:
    def __init__(self, cache):
        self.cache = Path(cache).resolve() / f'kew-{RELEASE}'

    def download(self, relative, *, compress=False):
        filename = Path(relative).name + ('.gz' if compress else '')
        return cached_download(f'{BASE_URL}/{relative}', self.cache / filename, compress=compress)

    def manifests(self):
        return [self.download(name) for name in ('gene_manifest.txt', TREE)]

    def gene_ids(self):
        rows = read_tsv(self.download('gene_manifest.txt'))
        genes = {row['Gene_ID'] for row in rows}
        if len(rows) != 353 or len(genes) != 353 or any(not LOCUS.fullmatch(gene) for gene in genes):
            raise ValueError('Kew gene manifest must contain the 353 unique numeric Angiosperms353 loci')
        return genes

    def loci(self, selected=None, workers=4):
        available = self.gene_ids()
        wanted = available if selected is None else set(selected)
        if wanted - available:
            raise ValueError(f'unknown Kew loci: {sorted(wanted - available)}')
        ordered = sorted(wanted)
        with ThreadPoolExecutor(max_workers=workers) as pool:
            return dict(zip(ordered, pool.map(self.marker, ordered), strict=True))

    def marker(self, locus):
        return self.download(f'fasta/by_gene/{locus}.dna.fasta', compress=True)


def kew_families(tree_path):
    """Read the family annotation in the released leaf labels, not tree edges."""
    # Kew uses bare labels with order_family_species_accession; ignore internal labels.
    by_recovery, by_genus = {}, defaultdict(set)
    for label in re.findall(r'(?<=[(,])([^(),:;\s]+)', Path(tree_path).read_text()):
        fields = label.split('_', 3)
        if len(fields) != 4:
            raise ValueError(f'unexpected Kew leaf label: {label}')
        _, family, genus, suffix = fields
        key = f'{genus}_{suffix}'
        if key in by_recovery and by_recovery[key] != family:
            raise ValueError(f'conflicting Kew family annotations: {key}')
        by_recovery[key] = family
        by_genus[genus].add(family)
    if not by_recovery:
        raise ValueError('no Kew family annotations found')
    return by_recovery, {genus: next(iter(families)) for genus, families in by_genus.items()
                         if len(families) == 1}


class ReferenceTaxonomy(Taxonomy):
    @lru_cache(maxsize=None)
    def name_ids(self, name, rank):
        return {row[0] for row in self.db.execute('''
            SELECT taxid FROM species WHERE spname=? AND rank=?
            UNION
            SELECT s.taxid FROM synonym n JOIN species s ON s.taxid=n.taxid
            WHERE n.spname=? COLLATE NOCASE AND s.rank=?
        ''', (name, rank, name, rank))}

    @lru_cache(maxsize=None)
    def reconcile(self, genus, family):
        genera = self.name_ids(genus, 'genus')
        families = self.name_ids(family, 'family') if family else set()
        if family and not families:
            return (None, None)
        if families:
            matching = {g for g in genera if self.lineage(g).get('family', (None,))[0] in families}
            if genera and not matching:
                return (None, None)
            genera = matching
        if len(genera) == 1:
            gid = next(iter(genera))
            return (self.lineage(gid).get('family', (None,))[0], gid)
        return (next(iter(families)) if len(families) == 1 else None, None)

    @lru_cache(maxsize=None)
    def reference_ranks(self, genus, family):
        family_id, genus_id = self.reconcile(genus, family)
        anchor = genus_id if genus_id is not None else family_id
        lineage = self.lineage(anchor) if anchor is not None else {}
        return tuple(family_id if rank == 'family' else lineage.get(rank, (None,))[0]
                     for rank in RANKS)


def reference_rows(path, locus, taxonomy, families):
    exact, genera = families
    for header, sequence in fasta(path):
        fields = dict(word.split(':', 1) for word in header.split()[1:] if ':' in word)
        if header.split()[0] != locus or not all(fields.get(k) for k in ('Species', 'Repository', 'Sequence_ID')):
            raise ValueError(f'unexpected FASTA header in Kew locus {locus}: {header}')
        species = fields['Species']
        accession = fields['Sequence_ID']
        genus = species.split('_', 1)[0]
        family = exact.get(f'{species}_{accession}', genera.get(genus))
        ranks = taxonomy.reference_ranks(genus, family)
        if ranks[0] is not None and sequence:
            yield species, fields['Repository'], accession, ranks, sequence
