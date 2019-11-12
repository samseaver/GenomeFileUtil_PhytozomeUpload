import unittest
import os
import time
from uuid import uuid4
from configparser import ConfigParser

from GenomeFileUtil.GenomeFileUtilServer import MethodContext
from GenomeFileUtil.GenomeFileUtilImpl import GenomeFileUtil
from installed_clients.WorkspaceClient import Workspace

_DATA_PATH = "/kb/module/test/data"
_TEST_GENOME = {
    "assembly_ref": "44970/1/1",
    "cdss": [
        {
            "aliases": [["protein_id", "BAC24147.1"], ["gene", "gidA"]],
            "db_xrefs": [["NCBI_GP", "BAC24147.1"]],
            "dna_sequence": "ATCG",
            "dna_sequence_length": 1887,
            "functions": ["BAC24147.1"],
            "id": "gene-gidA.CDS",
            "location": [["BA000021.3", 182, "+", 1887]],
            "md5": "5fda22c7efeab80748d2b91270d9ba10",
            "note": "Unknown function",
            "parent_gene": "gene-gidA",
            "protein_md5": "5569b50f7ff67531a4db47af662a44c8",
            "protein_translation": "PPP",
            "protein_translation_length": 628
        },
    ],
    "contig_ids": ["BA000021.3"],
    "contig_lengths": [697724],
    "dna_size": 697724,
    "feature_counts": {
        "CDS": 11,
        "gene": 11,
        "non_coding_features": 0,
        "protein_encoding_gene": 11
    },
    "features": [
        {
            "aliases": [["gene", "gidA"], ["protein_id", "BAC24147.1"]],
            "cdss": ["gene-gidA.CDS"],
            "dna_sequence": "ATCG",
            "dna_sequence_length": 1887,
            "functions": ["BAC24147.1"],
            "id": "gene-gidA",
            "location": [["BA000021.3", 182, "+", 1887]],
            "md5": "5fda22c7efeab80748d2b91270d9ba10",
            "protein_translation": "PPP",
            "protein_translation_length": 628
        },
    ],
    "gc_content": 0.22479,
    "genome_tiers": ["User"],
    "id": "a643f789-bf5d-4e08-a57f-52e853b15f75",
    "md5": "325ef1bc3886187384b1c483c4d3b11d",
    "molecule_type": "SingleLetterAlphabet",
    "mrnas": [],
    "non_coding_features": [],
    "num_contigs": 1,
    "ontologies_present": {},
    "ontology_events": [
        {
            "id": "GO",
            "method": "GenomeFileUtils Genbank uploader from annotations",
            "method_version": "0.9.0",
            "ontology_ref": "6308/3/2",
            "timestamp": "2019_11_11_23_11_46"
        }
    ],
    "source": "User",
    "source_id": "unknown",
    "taxon_assignments": {"ncbi": "3702"},
}


class IntegrationTest(unittest.TestCase):
    """
    Basic integration tests for save_one_genome.
    """

    @classmethod
    def setUpClass(cls):
        token = os.environ.get('KB_AUTH_TOKEN', None)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'provenance': [
                            {'service': 'GenomeFileUtil',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)  # type: ignore
        for nameval in config.items('GenomeFileUtil'):
            cls.cfg[nameval[0]] = nameval[1]
        cls.wsURL = cls.cfg['workspace-url']
        cls.ws = Workspace(cls.wsURL, token=token)
        cls.gfu = GenomeFileUtil(cls.cfg)
        # create one WS for all tests
        suffix = int(time.time() * 1000)
        cls.ws_name = "test_GenomeAnnotationAPI_" + str(suffix)
        ws_info = cls.ws.create_workspace({'workspace': cls.ws_name})
        cls.ws_id = ws_info[0]

    def test_save_one_genome_ok(self):
        """
        Test the `genbank_to_genome` method
        Taxonomy ID should be pulled from the genbank.
        """
        result = self.gfu.save_one_genome(self.ctx, {
            'workspace_name': self.ws_name,
            'data': _TEST_GENOME,
            'name': str(uuid4()),
        })
        print('xyz result', result)
        # ref = result[0]['genome_ref']
        # self.assertTrue(ref, 'Genome ref exists')
        # info = result[0]['genome_info']
        # typ = info[2]
        # self.assertTrue(typ.startswith('KBaseGenomes.Genome'))
        # info_details = info[-1]
        # self.assertEqual(info_details['Taxonomy'], (
        #     'cellular organisms;Bacteria;Proteobacteria;'
        #     'Gammaproteobacteria;Enterobacterales;Erwiniaceae;'
        #     'Wigglesworthia;Wigglesworthia glossinidia'
        # ))
        # self.assertEqual(info_details['Size'], '697724')
        # self.assertEqual(info_details['Source'], 'Genbank')
        # self.assertEqual(info_details['Name'], 'Wigglesworthia glossinidia endosymbiont of Glossina brevipalpis')
        # self.assertEqual(info_details['GC content'], '0.22479')
        # self.assertEqual(info_details['Genetic code'], '11')
        # self.assertEqual(info_details['Number of Genome Level Warnings'], '0')
        # self.assertEqual(info_details['Source ID'], 'BA000021')
        # self.assertEqual(info_details['Number of Protein Encoding Genes'], '20')
        # self.assertEqual(info_details['Domain'], 'Bacteria')
        # self.assertTrue(info_details['Assembly Object'])
        # self.assertEqual(info_details['Number contigs'], '1')
        # self.assertEqual(info_details['Number of CDS'], '20')
        # self.assertTrue(info_details['MD5'])
