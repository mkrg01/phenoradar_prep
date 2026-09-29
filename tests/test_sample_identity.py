"""Sample independence, one-time relabelling, and biological species joins."""
import gzip
import json
from pathlib import Path

import pytest

from common import read_tsv, write_tsv
from dataset_assets import identities
from sample_identity import sample_id, select_samples, traits_for_samples
from relabel_sample import relabel, text_open


def test_repeated_species_have_distinct_sample_keys(tmp_path):
    path=tmp_path/'metadata.tsv'
    rows=[dict(scientific_name='Plant alpha',run=r,taxid='42') for r in ('SRR1','LOCAL_read_2')]
    write_tsv(path,list(rows[0]),rows)
    fields,items=identities(path)
    assert {i['species'] for i in items}=={'Plant_alpha_SRR1','Plant_alpha_LOCAL_read_2'}
    assert {i['row']['species_id'] for i in items}=={'Plant_alpha'}
    assert 'analysis_sample_id' in fields
    samples=[dict(i['row'],species=i['species']) for i in items]
    assert select_samples(samples,['Plant_alpha'])=={i['species'] for i in items}
    assert select_samples(samples,['Plant_alpha_SRR1'])=={'Plant_alpha_SRR1'}
    assert traits_for_samples(samples,{'Plant_alpha':'1'})=={i['species']:'1' for i in items}


def test_concatenation_and_odb_collisions_are_rejected(tmp_path):
    path=tmp_path/'metadata.tsv'
    for rows in ([dict(scientific_name='Plant alpha',run='x_R1',taxid='1'),dict(scientific_name='Plant alpha x',run='R1',taxid='1')],
                 [dict(scientific_name='Plant alpha-x',run='R1',taxid='1'),dict(scientific_name='Plant alpha',run='x_R1',taxid='1')]):
        write_tsv(path,list(rows[0]),rows)
        with pytest.raises(ValueError,match='unique.*identities'): identities(path)


@pytest.mark.parametrize('kind,content',[
 ('fasta','>Plant_alpha_g1 description\nATGAAATAG\n>Plant_alpha_g22\nATGCCC\n'),
 ('busco','# The lineage dataset is: embryophyta_odb12\nM1\tComplete\tPlant_alpha_g1:1-9\t10\t3\nM2\tMissing\n'),
 ('quant','target_id\ttpm\nPlant_alpha_g1\t0.0000012300\nPlant_alpha_g22\t1e+06\n'),
 ('mapping','gene_id\torthogroup\nPlant_alpha_g1\tOG1\nPlant_alpha_g22\t\n')])
def test_relabel_is_lossless_and_resumable(tmp_path,kind,content):
    source=tmp_path/'source.gz';target=tmp_path/'target.gz'
    with gzip.open(source,'wt') as f:f.write(content)
    original=source.read_bytes()
    r=relabel(source,target,'Plant_alpha','Plant_alpha_SRR1',kind)
    with text_open(target) as f:assert f.read()==content.replace('Plant_alpha_g','Plant_alpha_SRR1_g')
    first=target.read_bytes();relabel(source,target,'Plant_alpha','Plant_alpha_SRR1',kind)
    assert target.read_bytes()==first and source.read_bytes()==original
    assert r['decoded_sha256']


def test_relabel_refuses_foreign_ids_and_existing_conflicts(tmp_path):
    source=tmp_path/'source.fa';target=tmp_path/'target.fa'
    source.write_text('>Other_g1\nATG\n');target.write_text('keep')
    with pytest.raises(ValueError,match='does not belong'):relabel(source,target,'Plant_alpha','Plant_alpha_SRR1')
    assert target.read_text()=='keep'
    source.write_text('>Plant_alpha_g1\nATG\n')
    with pytest.raises(ValueError,match='conflicting'):relabel(source,target,'Plant_alpha','Plant_alpha_SRR1')
    with pytest.raises(ValueError,match='separate'):relabel(source,source,'Plant_alpha','Plant_alpha_SRR1')


from test_datasets import dataset_project, new_dataset


def test_two_runs_of_one_species_execute_independent_native_stages(dataset_project):
    from dataset import submit,worker,status,materialize
    root=dataset_project
    build=new_dataset(root,('New plant','New plant'))
    submit(build,until='quant',dry_run=True)
    for stage in ('assembly','busco','quant'):worker(build,stage,1)
    states=status(build)
    assert states[0]['quant']=='reuse' and states[1]['assembly']=='pending'
    for stage in ('assembly','busco','quant'):worker(build,stage,2)
    assert all(r['quant']=='reuse' for r in status(build))
    inputs=materialize(build)
    rows=read_tsv(inputs/'metadata.tsv')
    assert {r['scientific_name'] for r in rows}=={'New plant'}
    assert {r['analysis_sample_id'] for r in rows}=={'New_plant_SRR1','New_plant_SRR2'}
    for row in rows:
        name=row['analysis_sample_id']
        events=(build/'work/genegalleon'/name/'events.jsonl').read_text().splitlines()
        assert len(events)==3
        with gzip.open(inputs/'cds'/f'{name}_longestCDS.fa.gz','rt') as f:assert f.readline().startswith('>'+name+'_g')
        quant=read_tsv(inputs/'quant'/name/row['run']/f'{row["run"]}_abundance.tsv')
        assert all(r['target_id'].startswith(name+'_g') for r in quant)


def test_alignment_min_taxa_counts_species_not_sample_replicates(tmp_path):
    from infer_phylogeny import alignment_qc
    manifest=tmp_path/'samples.tsv'
    rows=[dict(species=f'Plant_alpha_R{i}',species_id='Plant_alpha') for i in range(4)]
    write_tsv(manifest,list(rows[0]),rows)
    records=[(r['species'],'AAAAC' if i else 'AAAAD') for i,r in enumerate(rows)]
    output,qc=alignment_qc(records,dict(min_taxa=4,min_protein_length=1,sample_manifest=str(manifest)))
    assert qc['taxa']==4 and qc['biological_species']==1
    assert qc['status']=='too_few_taxa'


def test_timetree_keeps_same_species_samples_and_deduplicates_queries(tmp_path,tiny_inputs):
    from timetree_calibrations import species_taxids,interpret_response
    from test_timetree_calibrations import payload
    metadata=tmp_path/'samples.tsv'
    rows=[dict(species='Alpha_plant_A1',species_id='Alpha_plant',taxid='42'),
          dict(species='Alpha_plant_A2',species_id='Alpha_plant',taxid='42'),
          dict(species='Beta_sp-X_B1',species_id='Beta_sp-X',taxid='43')]
    write_tsv(metadata,list(rows[0]),rows)
    mapping,report=species_taxids(metadata,tiny_inputs['taxonomy_db'],{r['species'] for r in rows})
    assert len(mapping)==3 and set(mapping.values())=={42,43}
    query={'node':'N','clade_taxa':3,'query_taxa':sorted(mapping),
           'children':[['Alpha_plant_A1','Alpha_plant_A2'],['Beta_sp-X_B1']]}
    result=interpret_response(payload([42,43]),query,mapping)
    assert result['query_taxids']==[42,43] and len(result['used_taxa'])==3


def test_ncbi_guide_expands_duplicate_taxids(tmp_path,tiny_inputs,monkeypatch):
    import phylogeny_root
    from ete4 import Tree
    # Test the actual installed constrain implementation; version pinning is separately tested.
    import nwkit.constrain
    monkeypatch.setattr(phylogeny_root,'nwkit_backend',lambda:None)
    samples=tmp_path/'samples.tsv'
    rows=[dict(species='Alpha_plant_A1',taxid='42',cds='a'),dict(species='Alpha_plant_A2',taxid='42',cds='b'),
          dict(species='Beta_sp-X_B1',taxid='43',cds='c'),dict(species='Gamma_plant_G1',taxid='44',cds='d')]
    write_tsv(samples,list(rows[0]),rows)
    tree=tmp_path/'guide.nwk'
    phylogeny_root.ncbi_tree(samples,tiny_inputs['taxonomy_db'],tree,tmp_path/'taxids.tsv')
    t=Tree(tree.read_text(),parser=9)
    assert set(t.leaf_names())=={r['species'] for r in rows}
    assert set(t.common_ancestor(['Alpha_plant_A1','Alpha_plant_A2']).leaf_names())=={'Alpha_plant_A1','Alpha_plant_A2'}


def test_phenotyped_selection_and_tpm_keep_replicates_separate(tmp_path):
    from species_traits import select_phenotyped
    from export_species_tpm import export
    rows = [dict(species=f"Plant_{name}_R{i}", species_id=f"Plant_{name}", run=f"R{i}")
            for i, name in enumerate("AABCD")]
    samples, traits = tmp_path / "samples.tsv", tmp_path / "traits.tsv"
    write_tsv(samples, list(rows[0]), rows)
    write_tsv(traits, ["species", "C4"],
              [dict(species=f"Plant_{n}", C4=str(i % 2)) for i, n in enumerate("ABCD")])
    select_phenotyped(samples, traits, tmp_path / "selected")
    report = json.loads((tmp_path / "selected/selection.json").read_text())
    assert report["inference_species"] == 4 and report["inference_samples"] == 5
    assert len(read_tsv(tmp_path / "selected/samples.tsv")) == 5
    tpm, output = tmp_path / "tpm.tsv", tmp_path / "export.tsv"
    values = ["0.0000012300", "9e+05", "10", "20", "0"]
    write_tsv(tpm, ["species", "run", "orthogroup", "tpm"],
              [dict(species=r["species"], run=r["run"], orthogroup="OG1", tpm=v)
               for r, v in zip(rows, values)])
    export(samples, tpm, output)
    assert [r["tpm"] for r in read_tsv(output)] == values
    assert len({r["species"] for r in read_tsv(output)}) == 5


def test_real_phylogeny_retains_independent_samples(tmp_path, command_environment, workflow_project, seed_taxonomy):
    import os
    import shutil
    import subprocess
    import sys
    import yaml
    from ete4 import Tree
    from test_phylogeny import phylogeny_inputs, trimal_binary, ROOT
    from sample_identity import annotate

    famsa = os.environ.get("FAMSA_BIN") or shutil.which("famsa")
    vft = os.environ.get("VERYFASTTREE_BIN") or shutil.which("VeryFastTree")
    snakemake = shutil.which("snakemake")
    if not all([famsa, vft, snakemake]) or not (ROOT / "resources/phylogeny_tools/bin/astral4_int128").exists():
        pytest.skip("real inference tools required")
    source, biological = phylogeny_inputs(tmp_path)
    old_rows = read_tsv(source / "metadata.tsv")
    summaries = read_tsv(source / "busco/summary.tsv")
    # Copy a second independent assembly for one species, using distinct gene IDs.
    old_rows.append(dict(old_rows[1], run="R6"))
    metadata, busco = [], []
    for row in old_rows:
        item = annotate(row)
        old, new = item["species_id"], item["analysis_sample_id"]
        origin_run = "R1" if row["run"] == "R6" else row["run"]
        relabel(source / "cds" / f"{old}_longestCDS.fa.gz",
                source / "cds" / f"{new}_longestCDS.fa.gz", old, new)
        relabel(source / "busco/full" / f"{old}.busco.full.tsv",
                source / "busco/full" / f"{new}.busco.full.tsv", old, new, "busco")
        relabel(source / "quant" / old / origin_run / f"{origin_run}_abundance.tsv",
                source / "quant" / new / row["run"] / f"{row['run']}_abundance.tsv", old, new, "quant")
        metadata.append(item)
        busco.append(dict(summaries[biological.index(old)], Species=new))
    write_tsv(source / "metadata.tsv", list(metadata[0]), metadata)
    write_tsv(source / "busco/summary.tsv", list(busco[0]), busco)
    seed_taxonomy(source / "taxa.sqlite")
    config = tmp_path / "config.yaml"
    config.write_text(yaml.safe_dump(dict(run_name="test", phylogeny=dict(
        trees=["all"], outgroup=metadata[0]["analysis_sample_id"], max_markers=3))))
    env = command_environment(dict(python=sys.executable, famsa=famsa, trimal=trimal_binary(), VeryFastTree=vft))
    result = subprocess.run([snakemake, "--snakefile", str(ROOT / "workflow/Snakefile"),
        "--configfile", str(config), "--cores", "2", "--resources", "mem_mb=8000",
        "--set-threads", "align_busco_marker=1", "infer_busco_gene_tree=1", "infer_busco_species_tree=2",
        "--set-resources", "infer_busco_species_tree:mem_mb=4000", "--", "phylogeny"],
        cwd=workflow_project, env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    assert result.returncode == 0, result.stdout + "\n" + "\n".join(
        p.read_text()[-3000:] for p in (tmp_path / "logs").rglob("*.log"))
    out = tmp_path / "results/test"
    assert set(Tree((out / "phylogeny/all/species_tree.nwk").read_text()).leaf_names()) == {
        r["analysis_sample_id"] for r in metadata}
    selected = json.loads((out / "metadata/selection.json").read_text())
    assert selected["selected_species"] == 6 and selected["selected_samples"] == 7
    assert len(read_tsv(out / "phylogeny/all/plan/species.tsv")) == 7


def test_failed_registration_can_retry_with_new_output(dataset_project, monkeypatch):
    import dataset
    build = new_dataset(dataset_project)
    dataset.submit(build, until="quant", dry_run=True)
    dataset.worker(build, "assembly", 1)
    dataset.worker(build, "busco", 1)
    register = dataset.register_quant
    def fail(*args, **kwargs):
        raise ValueError("interrupted registration")
    monkeypatch.setattr(dataset, "register_quant", fail)
    with pytest.raises(ValueError, match="interrupted registration"):
        dataset.worker(build, "quant", 1)
    staged = build / "work/genegalleon/New_plant_SRR1/products/SRR1_abundance.tsv"
    staged.write_text("previous failed attempt\n")
    monkeypatch.setattr(dataset, "register_quant", register)
    dataset.worker(build, "quant", 1)
    assert dataset.status(build)[0]["quant"] == "reuse"
    assert read_tsv(staged)[0]["target_id"] == "New_plant_SRR1_g1"
    assert any(p.read_text() == "previous failed attempt\n" for p in
               (build / "jobs/incomplete/New_plant_SRR1/quant").rglob("products-*.tsv"))


def test_posthoc_filter_retains_species_traits_for_remaining_samples(tmp_path):
    from filter_species import export
    source = tmp_path / "analysis"
    rows = [dict(species=f"Plant_{n}_R{i}", species_id=f"Plant_{n}",
                 scientific_name=f"Plant {n}", run=f"R{i}", taxid="42" if n == "A" else "43",
                 odb_species=f"Plant_{n}_R{i}", family="Plantaceae") for i, n in enumerate("AAB")]
    write_tsv(source / "metadata/samples.tsv", list(rows[0]), rows)
    write_tsv(source / "metadata/metadata_high_busco.tsv", list(rows[0]), rows)
    traits = tmp_path / "traits.tsv"
    write_tsv(traits, ["species", "C4"], [dict(species="Plant_A", C4="0"), dict(species="Plant_B", C4="1")])
    export(source, ["Plant_A_R0"], traits=traits)
    filtered = source / "filtered/metadata"
    assert {r["species"] for r in read_tsv(filtered / "species_trait.tsv")} == {"Plant_A", "Plant_B"}
    assert {r["species"]: r["C4"] for r in read_tsv(filtered / "species_metadata.tsv")} == {
        "Plant_A_R1": "0", "Plant_B_R2": "1"}
    export(source, ["Plant_A_R0", "Plant_A_R1"], traits=traits)
    assert [r["species"] for r in read_tsv(filtered / "species_trait.tsv")] == ["Plant_B"]


def test_normalized_biological_names_require_consistent_taxids(tmp_path):
    path = tmp_path / "metadata.tsv"
    rows = [dict(scientific_name="Plant alpha", run="R1", taxid="42"),
            dict(scientific_name="Plant_alpha", run="R2", taxid="43")]
    write_tsv(path, list(rows[0]), rows)
    with pytest.raises(ValueError, match="conflicting taxids"):
        identities(path)
