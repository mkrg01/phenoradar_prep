from xml.etree import ElementTree

import pytest

from ncbi_sample import NCBISamples


@pytest.mark.parametrize("bio_id", ["11", "22"])
def test_confirm_matches_run_and_biosample_species(monkeypatch, bio_id):
    from types import SimpleNamespace

    xml = {
        "sra": """<EXPERIMENT_PACKAGE_SET><EXPERIMENT_PACKAGE>
          <SAMPLE><IDENTIFIERS><EXTERNAL_ID>SAMN1</EXTERNAL_ID></IDENTIFIERS>
          <SAMPLE_NAME><TAXON_ID>12</TAXON_ID></SAMPLE_NAME></SAMPLE>
          <RUN_SET><RUN accession="SRR1"/></RUN_SET>
          </EXPERIMENT_PACKAGE></EXPERIMENT_PACKAGE_SET>""",
        "biosample": f"""<BioSampleSet><BioSample accession="SAMN1">
          <Description><Organism taxonomy_id="{bio_id}"/></Description>
          </BioSample></BioSampleSet>""",
    }
    ncbi = NCBISamples()
    monkeypatch.setattr(ncbi, "fetch", lambda database, accession: ElementTree.fromstring(xml[database]))
    taxonomy = SimpleNamespace(resolve=lambda tid: (
        tid, {"species": ("11" if tid == "12" else tid, "Species")}))
    if bio_id == "11":
        evidence = ncbi.confirm("SRR1", "11", taxonomy)
        assert (evidence["run"], evidence["biosample"], evidence["sra_taxid"]) == ("SRR1", "SAMN1", "12")
    else:
        with pytest.raises(ValueError, match="BioSample organism"):
            ncbi.confirm("SRR1", "11", taxonomy)
