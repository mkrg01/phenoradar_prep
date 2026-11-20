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
LAND_PLANTS_TAXID = 3193
LEGACY_RANKS = ('family', 'genus')
INTERMEDIATE_RANKS = ('family', 'subfamily', 'supertribe', 'tribe', 'subtribe', 'genus')
RANKS = (*INTERMEDIATE_RANKS, 'superfamily', 'parvorder', 'infraorder', 'suborder',
         'order', 'superorder', 'infraclass', 'subclass', 'class', 'superclass',
         'subphylum', 'phylum', 'superphylum', 'subkingdom', 'kingdom', 'superkingdom', 'domain')
ASSESSMENT_RANKS = frozenset(RANKS[:RANKS.index('order') + 1])


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
    def __init__(self, path):
        from ete4 import NCBITaxa

        path = Path(path).resolve()
        if not path.is_file():
            raise ValueError(f'taxonomy database missing: {path}')
        self.ncbi = NCBITaxa(dbfile=str(path), update=False)
        self.db = self.ncbi.db
        self.db.execute('PRAGMA query_only=ON')

    @lru_cache(maxsize=None)
    def ancestor_ids(self, taxid):
        node = self.node(taxid)
        return tuple(self.ncbi.get_lineage(node[0])) if node is not None else ()

    @lru_cache(maxsize=None)
    def name_ids(self, name, rank):
        name = ' '.join((name or '').replace('_', ' ').split())
        if not name:
            return ()
        candidates = set()
        for (taxid,) in self.db.execute('''
            SELECT taxid FROM species WHERE spname=? COLLATE NOCASE
            UNION
            SELECT taxid FROM synonym WHERE spname=? COLLATE NOCASE
        ''', (name, name)):
            node = self.node(taxid)
            if node is not None and node[2] == rank and LAND_PLANTS_TAXID in self.ancestor_ids(node[0]):
                candidates.add(node[0])
        return tuple(sorted(candidates))

    def resolve_reference(self, species, family):
        name = ' '.join(species.replace('_', ' ').split())
        genus = name.split()[0] if name else ''
        for rank, query in [('species', name), ('genus', genus), ('family', family)]:
            candidates = self.name_ids(query, rank)
            if not candidates:
                continue
            paths = [self.ancestor_ids(taxid) for taxid in candidates]
            anchor = next(taxid for taxid in reversed(paths[0]) if all(taxid in path for path in paths[1:]))
            lineage = self.lineage(anchor)
            return dict(method=rank if len(candidates) == 1 else f'ambiguous_{rank}',
                        matched_rank=rank, taxid=anchor, assigned_rank=self.node(anchor)[2],
                        candidates=candidates, ranks=[lineage.get(level, (None,))[0] for level in RANKS])
        return dict(method='unresolved', matched_rank=None, taxid=None, assigned_rank=None,
                    candidates=(), ranks=[None] * len(RANKS))


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
        if any(taxid is not None for taxid in ranks) and sequence:
            yield species, repository, accession, ranks, sequence
