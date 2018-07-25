import unittest
from configparser import ConfigParser
from os import environ
import logging

import mock

from GenomeFileUtil.GenomeFileUtilImpl import GenomeFileUtil, SDKConfig
from GenomeFileUtil.GenomeFileUtilServer import MethodContext
from GenomeFileUtil.core.GenomeInterface import GenomeInterface
from GenomeFileUtil.core import GenomeUtils
from Workspace.WorkspaceClient import Workspace as workspaceService


class GenomeFileUtilTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        logging.basicConfig(level=logging.info)
        token = environ.get('KB_AUTH_TOKEN', None)
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
        config_file = environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('GenomeFileUtil'):
            cls.cfg[nameval[0]] = nameval[1]
        cls.wsURL = cls.cfg['workspace-url']
        cls.ws = workspaceService(cls.wsURL, token=token)
        cls.serviceImpl = GenomeFileUtil(cls.cfg)
        gi_config = SDKConfig(cls.cfg)
        cls.genome_interface = GenomeInterface(gi_config)

    @classmethod
    def tearDownClass(cls):
        pass

    def test_retreve_taxon(self):
        self.assertEqual(self.genome_interface.retrieve_taxon("meh", "Arabidopsis thaliana"),
                         ('cellular organisms; Eukaryota; Viridiplantae; Streptophyta; Streptophytina; Embryophyta; Tracheophyta; Euphyllophyta; Spermatophyta; Magnoliophyta; Mesangiospermae; eudicotyledons; Gunneridae; Pentapetalae; rosids; malvids; Brassicales; Brassicaceae; Camelineae; Arabidopsis',
                          'meh/3702_taxon', 'Eukaryota', 11))
        self.assertEqual(self.genome_interface.retrieve_taxon("meh", "Escherichia coli"),
                         ('cellular organisms; Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacterales; Enterobacteriaceae; Escherichia',
                          'meh/562_taxon', 'Bacteria', 11))
        self.assertEqual(self.genome_interface.retrieve_taxon("meh", "rhodobacter"),
                         ('Unconfirmed Organism: rhodobacter',
                          'ReferenceTaxons/unknown_taxon', 'Unknown', 11)
                         )
        self.assertEqual(self.genome_interface.retrieve_taxon("meh", "foo"),
                         ('Unconfirmed Organism: foo',
                          'ReferenceTaxons/unknown_taxon', 'Unknown', 11))

    @mock.patch("GenomeFileUtil.core.GenomeInterface.MAX_GENOME_SIZE", 1)
    def test_max_genome_size(self):
        with self.assertRaisesRegex(ValueError, "genome exceeds the maximum permitted size"):
            GenomeInterface.validate_genome({"taxon_ref": "", "domain": ""})

    def test_user(self):
        self.assertEqual(GenomeInterface.determine_tier('RefSeq user'),
                         ('RefSeq', ['ExternalDB', 'User']))

    def test_is_parent(self):
        gene_1 = {"type": "gene", "location": [["A", 100, "+", 400]]}
        gene_2 = {"type": "gene", "location": [["A", 500, "-", 400]]}
        gene_3 = {"type": "gene", "location": [["B", 100, "+", 400]]}
        mrna_1 = {"type": "mRNA", "location": [["A", 100, "+", 50], ["A", 175, "+", 75],
                                               ["A", 300, "+", 75], ["A", 400, "+", 100]]}
        mrna_2 = {"type": "mRNA", "location": [["A", 500, "-", 50], ["A", 425, "-", 75],
                                               ["A", 300, "-", 75], ["A", 200, "-", 100]]}
        mrna_3 = {"type": "mRNA", "location": [["B", 100, "+", 100]]}
        cds_1 = {"type": "mRNA", "location": [["A", 100, "+", 50], ["A", 175, "+", 75],
                                              ["A", 300, "+", 75], ["A", 400, "+", 50]]}
        cds_1a = {"type": "CDS", "location": [["A", 200, "+", 50], ["A", 300, "+", 75],
                                              ["A", 400, "+", 50]]}
        cds_2 = {"type": "CDS", "location": [["A", 500, "-", 50], ["A", 425, "-", 75],
                                             ["A", 300, "-", 75], ["A", 200, "-", 50]]}
        cds_2a = {"type": "CDS", "location": [["A", 400, "-", 50], ["A", 300, "-", 75],
                                              ["A", 200, "-", 50]]}
        cds_3 = {"type": "CDS", "location": [["B", 125, "+", 50]]}
        cds_3a = {"type": "CDS", "location": [["B", 150, "+", 50], ["B", 175, "+", 75],
                                             ["B", 300, "+", 75], ["B", 400, "+", 50]]}
        cds_4 = {"type": "mRNA", "location": [["A", 100, "+", 50], ["A", 175, "+", 50],
                                              ["A", 300, "+", 75], ["A", 400, "+", 50]]}
        cds_5 = {"type": "mRNA", "location": [["A", 100, "+", 50], ["A", 400, "+", 50]]}

        # Test gene parentage
        self.assertTrue(GenomeUtils.is_parent(gene_1, mrna_1))
        self.assertTrue(GenomeUtils.is_parent(gene_1, cds_1))

        self.assertTrue(GenomeUtils.is_parent(gene_2, mrna_2))
        self.assertTrue(GenomeUtils.is_parent(gene_2, cds_2))

        self.assertTrue(GenomeUtils.is_parent(gene_3, mrna_3))
        self.assertTrue(GenomeUtils.is_parent(gene_3, cds_3))

        self.assertFalse(GenomeUtils.is_parent(gene_1, mrna_2))
        self.assertFalse(GenomeUtils.is_parent(gene_1, cds_2))
        self.assertFalse(GenomeUtils.is_parent(gene_1, mrna_3))
        self.assertFalse(GenomeUtils.is_parent(gene_1, cds_3))

        # test mrna parentage
        self.assertTrue(GenomeUtils.is_parent(mrna_1, cds_1))
        self.assertTrue(GenomeUtils.is_parent(mrna_1, cds_1a))

        self.assertTrue(GenomeUtils.is_parent(mrna_2, cds_2))
        self.assertTrue(GenomeUtils.is_parent(mrna_2, cds_2a))

        self.assertTrue(GenomeUtils.is_parent(mrna_3, cds_3))

        self.assertFalse(GenomeUtils.is_parent(mrna_1, cds_2))
        self.assertFalse(GenomeUtils.is_parent(mrna_1, cds_3))
        self.assertFalse(GenomeUtils.is_parent(mrna_1, cds_3a))
        self.assertFalse(GenomeUtils.is_parent(mrna_1, cds_4))
        self.assertFalse(GenomeUtils.is_parent(mrna_1, cds_5))
