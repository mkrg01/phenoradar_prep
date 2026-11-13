"""Confirm the organism registered for a public SRA Run and its BioSample."""
import os
import re
import time
from urllib.parse import urlencode
from urllib.request import urlopen
from xml.etree import ElementTree

from common import now

ENDPOINT = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'


class NCBISamples:
    def __init__(self):
        self.next_request = 0

    def fetch(self, database, accession):
        params = dict(db=database, id=accession, retmode='xml', tool='phenoradar_prep')
        if email := os.environ.get('NCBI_EMAIL'):
            params['email'] = email
        url = ENDPOINT + '?' + urlencode(params)
        for attempt in range(2):
            time.sleep(max(0, self.next_request - time.monotonic()))
            self.next_request = time.monotonic() + 0.4
            try:
                with urlopen(url, timeout=15) as response:
                    content = response.read(4 * 1024 * 1024 + 1)
                if len(content) > 4 * 1024 * 1024:
                    raise ValueError('NCBI response exceeded the metadata size limit')
                return ElementTree.fromstring(content)
            except (OSError, ElementTree.ParseError) as error:
                if attempt:
                    raise ValueError(f'could not verify {database} metadata for {accession}: {error}') from error
                self.next_request = time.monotonic() + 1

    def confirm(self, run, expected_species, taxonomy):
        if not re.fullmatch(r'[SED]RR[0-9]+', run):
            raise ValueError('reassignment requires a public SRA Run accession')
        sra = self.fetch('sra', run)
        packages = [package for package in sra.iter('EXPERIMENT_PACKAGE')
                    if any(entry.get('accession') == run for entry in package.findall('./RUN_SET/RUN'))]
        if len(packages) != 1:
            raise ValueError('SRA Run does not resolve to one experiment package')
        sample = packages[0].find('./SAMPLE')
        if sample is None:
            raise ValueError('SRA Run has no linked sample')
        accessions = {entry.text for entry in sample.findall('./IDENTIFIERS/EXTERNAL_ID')
                      if entry.text and re.fullmatch(r'SAM[END][A-Z]?[0-9]+', entry.text)}
        if len(accessions) != 1:
            raise ValueError('SRA Run does not identify one BioSample')
        accession, = accessions
        biosample = self.fetch('biosample', accession)
        records = [entry for entry in biosample.iter('BioSample') if entry.get('accession') == accession]
        organism = records[0].find('./Description/Organism') if len(records) == 1 else None
        sra_id = sample.findtext('./SAMPLE_NAME/TAXON_ID')
        bio_id = organism.get('taxonomy_id') if organism is not None else None
        for source, taxid in (('SRA', sra_id), ('BioSample', bio_id)):
            resolved = taxonomy.resolve(taxid) if taxid and taxid.isdigit() else None
            if not resolved or resolved[1].get('species', ('', ''))[0] != expected_species:
                raise ValueError(f'{source} organism does not resolve to the requested species')
        return {'run': run, 'biosample': accession, 'sra_taxid': sra_id,
                'biosample_taxid': bio_id, 'species_taxid': expected_species, 'checked_at': now(),
                'sra_url': f'https://www.ncbi.nlm.nih.gov/sra/{run}',
                'biosample_url': f'https://www.ncbi.nlm.nih.gov/biosample/{accession}'}
