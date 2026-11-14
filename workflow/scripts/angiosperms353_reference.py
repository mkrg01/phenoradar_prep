"""Cache Kew Angiosperms353 markers and resolve reference taxonomy offline."""
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from contextlib import closing, nullcontext
from functools import lru_cache
import fcntl
import gzip
import re
from pathlib import Path
from urllib.request import urlopen

from common import atomic_writer, read_tsv, text_open, write_json
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


    @lru_cache(maxsize=None)
    def species_ids(self, name):
        candidates = set()
        for table in ('species', 'synonym'):
            for (taxid,) in self.db.execute(f'SELECT taxid FROM {table} WHERE spname=? COLLATE NOCASE', (name,)):
                node = self.node(taxid)
                if node is not None and node[2] == 'species':
                    candidates.add(node[0])
        return tuple(sorted(candidates))

    def resolve_reference(self, species, family):
        name = ' '.join(species.replace('_', ' ').split())
        candidates = self.species_ids(name)
        if len(candidates) == 1:
            taxid = candidates[0]
            lineage = self.lineage(taxid)
            family_ids = self.name_ids(family, 'family') if family else set()
            family_conflict = bool(family_ids and lineage.get('family', (None,))[0] not in family_ids)
            return dict(method='species', taxid=taxid, family_conflict=family_conflict,
                        ranks=[lineage.get(rank, (None,))[0] for rank in RANKS])
        if candidates:
            lineages = [self.lineage(taxid) for taxid in candidates]
            family_ids = {lineage.get('family', (None,))[0] for lineage in lineages}
            expected = self.name_ids(family, 'family') if family else family_ids
            # Retain a common family only if the references agree on it.
            family_id = next(iter(family_ids)) if len(family_ids) == 1 and family_ids <= expected else None
            lineage = self.lineage(family_id) if family_id is not None else {}
            return dict(method='ambiguous_species', candidates=candidates,
                        ranks=[lineage.get(rank, (None,))[0] for rank in RANKS])
        ranks = self.reference_ranks(species.split('_', 1)[0], family)
        method = 'genus' if ranks[RANKS.index('genus')] is not None else 'family' if ranks[0] else 'unresolved'
        return dict(method=method, ranks=ranks)


def reference_label(header, locus):
    parts = header.split()
    gene, tokens = (parts[0], parts[1:]) if parts else ('', [])
    fields = dict(word.split(':', 1) for word in tokens if ':' in word)
    keys = ('Species', 'Repository', 'Sequence_ID')
    if gene != locus or not all(fields.get(key) for key in keys):
        raise ValueError(f'unexpected FASTA header in Kew locus {locus}: {header}')
    return tuple(fields[key] for key in keys)


def reference_key(species, repository, accession):
    return '|'.join((species, repository, accession))


def prepare_reference_taxonomy(paths, taxonomy_path, tree_path, output):
    """Resolve the union of FASTA reference labels once, before running any batch."""
    labels = set()
    for locus, path in paths.items():
        with text_open(path) as handle:
            labels.update(reference_label(line[1:], locus) for line in handle if line.startswith('>'))
    if not labels:
        raise ValueError('Kew reference files contain no sequence labels')
    exact, genera = kew_families(tree_path)
    tax = ReferenceTaxonomy(taxonomy_path)
    with closing(tax.db):
        resolved = {}
        for species, repository, accession in sorted(labels):
            genus = species.split('_', 1)[0]
            family = exact.get(f'{species}_{accession}', genera.get(genus))
            resolved[reference_key(species, repository, accession)] = dict(
                kew_family=family, **tax.resolve_reference(species, family))
    write_json(output, resolved, compact=True)


def reference_rows(path, locus, assignments):
    for header, sequence in fasta(path):
        species, repository, accession = reference_label(header, locus)
        ranks = assignments[reference_key(species, repository, accession)]['ranks']
        if ranks[0] is not None and sequence:
            yield species, repository, accession, ranks, sequence
