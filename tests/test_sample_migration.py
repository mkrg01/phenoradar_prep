"""A migrated registration is reusable without invoking any scientific tools."""
import gzip
import json
from pathlib import Path

import pytest

from common import read_tsv, write_json, write_tsv
from dataset_assets import COUNTS, digest, identities, record, resolve
from mapping_tables import read_species, write_tables
from migrate_sample_ids import migrate
from protein_cache import cached_translate


@pytest.fixture
def legacy_dataset(tmp_path):
    source=tmp_path/'old';bundle=source/'products';store=source/'store'
    old='Plant_alpha';run='SRR1';row=dict(scientific_name='Plant alpha',run=run,taxid='42')
    bundle.mkdir(parents=True)
    cds=bundle/'cds'/f'{old}_longestCDS.fa.gz';cds.parent.mkdir()
    with gzip.open(cds,'wt') as f:
        for i in range(3):f.write(f'>{old}_g{i}\nATGAAATAA\n')
    full=bundle/'busco/full'/f'{old}.busco.full.tsv';full.parent.mkdir(parents=True)
    full.write_text(f'# The lineage dataset is: embryophyta_odb12\nM1\tComplete\t{old}_g0:1-9\t100\t3\n')
    abundance=bundle/'quant'/old/run/f'{run}_abundance.tsv'
    write_tsv(abundance,['target_id','tpm'],[dict(target_id=f'{old}_g{i}',tpm=v) for i,v in enumerate(['0.0000012300','20','1e+06'])])
    protein=bundle/'proteins'/f'{old}_protein.fa';protein.parent.mkdir()
    protein.write_text(''.join(f'>{old}_g{i}\nMK*\n' for i in range(3)))
    translation=dict(cds=record(cds),protein=record(protein),seqkit='synthetic',translation_table=1,sequences=3)
    write_json(protein.with_suffix('.json'),translation)
    ref=dict(schema_version=1,species=old,taxid='42',reference_id=record(cds)['sha256'],cds=record(cds),provenance={'source':'legacy'})
    root=store/old/ref['reference_id'];write_json(root/'reference.json',ref)
    counts=dict(zip(COUNTS,[1,0,0,0,1]))
    write_json(root/'busco.json',dict(schema_version=1,reference_id=ref['reference_id'],lineage='embryophyta_odb12',counts=counts,full=record(full),short=None))
    write_json(root/'quant'/f'{run}.json',dict(schema_version=1,reference_id=ref['reference_id'],run=run,abundance=record(abundance),sample={}))
    table=source/'table.tsv.gz'
    with gzip.open(table,'wt') as f:f.write(f'gene_id\torthogroup\n{old}_g0\tOG1\n{old}_g1\tOG2\n{old}_g2\t\n')
    entry=dict(odb_species=old,protein_sha256=record(protein)['sha256'],table=record(table),
               qc=dict(protein_genes=3,unique_gene_og_pairs=2,ambiguous_genes=0),annotation_sha256='a'*64)
    mapping=write_tables(bundle/'odb',{old:entry})
    def relative(path):return dict(record(path),path=str(path.relative_to(bundle)))
    product=dict(row=row,reference_id=ref['reference_id'],protein=relative(protein),translation=relative(protein.with_suffix('.json')))
    data=dict(schema_version=3,kind='completed_build',products={old:product},mapping=relative(mapping),
              lineage='embryophyta_odb12',translation={'table':1},odb={'version':'v12','node':3193})
    data['sha256']=digest(data);write_json(bundle/'manifest.json',data)
    metadata=source/'metadata.tsv';write_tsv(metadata,list(row),[row])
    return dict(metadata=metadata,store=store,products=bundle,legacy_results=source/'unused',old_cache=source/'cache',output=tmp_path/'new')


def test_one_time_migration_preserves_data_and_seeds_caches(legacy_dataset):
    paths=legacy_dataset
    before={str(p):p.read_bytes() for p in paths['metadata'].parent.rglob('*') if p.is_file()}
    result=migrate(**paths,workers=1)
    out=paths['output'];name='Plant_alpha_SRR1'
    assert set(result)=={name}
    assert json.loads((out/'completed.json').read_text())['scientific_recomputation'] is False
    item=identities(out/'input/metadata.tsv')[1][0]
    products=resolve(out/'dataset_assets',item,'embryophyta_odb12',need_full=True)
    assert all(products[k] for k in ('reference','busco','quant'))
    assert [r['tpm'] for r in read_tsv(products['quant']['abundance']['path'])]==['0.0000012300','20','1e+06']
    genes,_=read_species(out/'mapping/snapshot.json',name)
    assert genes=={name+'_g0':['OG1'],name+'_g1':['OG2'],name+'_g2':[]}
    # An invalid executable proves that migration seeded the translation cache.
    assert cached_translate(products['reference']['cds']['path'],out/'reused.fa',out/'reused.json',
                            out/'dataset_assets/.proteins',seqkit='must-not-be-executed')
    second=migrate(**paths,workers=1)
    assert result.keys()==second.keys()
    assert before=={str(p):p.read_bytes() for p in paths['metadata'].parent.rglob('*') if p.is_file()}


def test_pilot_does_not_publish_completion_and_changed_sources_fail(legacy_dataset):
    paths=legacy_dataset
    migrate(**paths,workers=1,limit=1)
    assert not (paths['output']/'completed.json').exists()
    metadata=read_tsv(paths['metadata']);metadata[0]['run']='DIFFERENT'
    write_tsv(paths['metadata'],list(metadata[0]),metadata)
    with pytest.raises(ValueError,match='association unavailable'):
        migrate(**paths,workers=1)


def test_publish_migration_produces_reusable_schema4_build_without_tools(legacy_dataset,tiny_inputs,monkeypatch):
    import shutil
    import subprocess
    import yaml
    from publish_sample_migration import publish
    from build_products import load_complete
    from dataset import plan
    paths=legacy_dataset
    root=paths['output'].parent
    repo=Path(__file__).resolve().parents[1]
    for name in ('workflow','config','profiles'):
        shutil.copytree(repo/name,root/name,ignore=shutil.ignore_patterns('__pycache__'))
    shutil.copy2(repo/'run_pipeline.sh',root/'run_pipeline.sh')
    cfg=yaml.safe_load((root/'config/build.yaml').read_text())
    cfg['metadata']='old/metadata.tsv';cfg['excluded_accessions']=None
    (root/'config/build.yaml').write_text(yaml.safe_dump(cfg))
    database=root/'resources/taxonomy/taxa.sqlite';database.parent.mkdir(parents=True)
    shutil.copy2(tiny_inputs['taxonomy_db'],database)
    migrate(**paths,workers=1)
    monkeypatch.setattr(subprocess,'run',lambda *a,**k:pytest.fail('scientific command executed'))
    build=publish(root,paths['output'],'converted')
    done=load_complete(build)
    assert done['schema_version']==4 and set(done['products'])=={'Plant_alpha_SRR1'}
    report=plan(root,paths['output']/'build.yaml')[-1]
    assert [report[0][s] for s in ('assembly','busco','quant')]==['reuse']*3
    assert json.loads((paths['output']/'build_validation.json').read_text())['mapping_reused']==1
