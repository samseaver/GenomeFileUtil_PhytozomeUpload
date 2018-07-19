import unittest
from configparser import ConfigParser
from os import environ

import mock

from GenomeFileUtil.GenomeFileUtilImpl import GenomeFileUtil, SDKConfig
from GenomeFileUtil.GenomeFileUtilServer import MethodContext
from GenomeFileUtil.core.GenomeInterface import GenomeInterface
from Workspace.WorkspaceClient import Workspace as workspaceService


class GenomeFileUtilTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
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
