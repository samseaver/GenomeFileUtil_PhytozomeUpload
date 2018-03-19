import unittest
import os
import json
import time
import shutil
import filecmp

from os import environ
try:
    from ConfigParser import ConfigParser  # py2
except:
    from configparser import ConfigParser  # py3

from Workspace.WorkspaceClient import Workspace as workspaceService
from GenomeFileUtil.GenomeFileUtilImpl import GenomeFileUtil, SDKConfig
from GenomeFileUtil.GenomeFileUtilServer import MethodContext
from GenomeFileUtil.core.GenomeInterface import GenomeInterface


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
        self.assertEqual(self.genome_interface.retrieve_taxon("meh", "Escherichia coli"),
                         ('Escherichia coli', u'meh/562_taxon', u'Bacteria', 11))
        print(self.genome_interface.old_retrieve_taxon("meh", "Escherichia coli"))
        self.assertEqual(self.genome_interface.retrieve_taxon("meh", "rhodobacter"),
                         ('Unconfirmed Organism: rhodobacter',
                          'ReferenceTaxons/unknown_taxon', 'Unknown', 11)
                         )
        self.assertEqual(self.genome_interface.retrieve_taxon("meh", "foo"),
                         ('Unconfirmed Organism: foo',
                          'ReferenceTaxons/unknown_taxon', 'Unknown', 11))
