"""Small species-table fixtures, including deliberate damage for validation tests."""
import csv
import gzip
import json
from pathlib import Path
from common import write_json
from dataset_assets import record
from mapping_tables import write_tables


def make_mapping(path, genes, pairs):
    path = Path(path); staging = path.parent / '.fixtures'; staging.mkdir(parents=True, exist_ok=True)
    entries = {}
    for species in sorted({s for g,s in genes}):
        values = {g:[] for g,s in genes if s == species}
        for gene,og in pairs:
            if gene in values: values[gene].append(og)
        table = staging / (species + '.tsv.gz')
        # Replace rather than truncate links held by older snapshots.
        temporary = table.with_suffix('.tmp')
        with gzip.open(temporary, 'wt') as f:
            f.write('gene_id\torthogroup\n')
            for gene, groups in values.items():
                for group in groups or ['']: f.write(f'{gene}\t{group}\n')
        temporary.replace(table)
        n = sum(map(len, values.values()))
        entries[species] = {'odb_species':species.replace('-','_'), 'protein_sha256':'0'*64,
                            'table':record(table), 'qc':{'protein_genes':len(values),
                            'unique_gene_og_pairs':n,'input_annotation_rows':n,
                            'duplicate_pairs_removed':0,'ambiguous_genes':sum(len(v)>1 for v in values.values())}}
    result = write_tables(path.parent, entries)
    if path.name != result.name: result.replace(path)
    return path


def mapping_rows(path):
    path=Path(path); data=json.loads(path.read_text()); genes=[]; pairs=[]
    for species,entry in data['tables'].items():
        seen=set()
        with gzip.open(path.parent/entry['table']['path'],'rt') as f:
            for row in csv.DictReader(f,delimiter='\t'):
                gene=row['gene_id']
                if gene not in seen: genes.append((gene,species)); seen.add(gene)
                if row['orthogroup']:pairs.append((gene,row['orthogroup']))
    return genes,pairs


def edit_mapping(path, genes=None, pairs=None):
    old_genes,old_pairs=mapping_rows(path)
    return make_mapping(path, genes(old_genes) if genes else old_genes, pairs(old_pairs) if pairs else old_pairs)
