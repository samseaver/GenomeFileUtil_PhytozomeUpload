# -*- coding: utf-8 -*-
import unittest
import os  # noqa: F401
import json  # noqa: F401
import time
import shutil
import re

from os import environ
try:
    from ConfigParser import ConfigParser  # py2
except:
    from configparser import ConfigParser  # py3

from pprint import pprint  # noqa: F401

from biokbase.workspace.client import Workspace as workspaceService
from GenomeFileUtil.GenomeFileUtilImpl import GenomeFileUtil
from GenomeFileUtil.GenomeFileUtilImpl import SDKConfig
from GenomeFileUtil.GenomeFileUtilServer import MethodContext
from GenomeFileUtil.authclient import KBaseAuth as _KBaseAuth
from GenomeFileUtil.core.FastaGFFToGenome import FastaGFFToGenome
from DataFileUtil.DataFileUtilClient import DataFileUtil


class FastaGFFToGenomeUploadTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        print('setting up class')
        cls.token = environ.get('KB_AUTH_TOKEN', None)
        config_file = environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('GenomeFileUtil'):
            cls.cfg[nameval[0]] = nameval[1]
        authServiceUrl = cls.cfg.get('auth-service-url',
                                     "https://kbase.us/services/authorization/Sessions/Login")
        auth_client = _KBaseAuth(authServiceUrl)
        cls.user_id = auth_client.get_user(cls.token)
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': cls.token,
                        'user_id': cls.user_id,
                        'provenance': [
                            {'service': 'GenomeFileUtil',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = workspaceService(cls.wsURL, token=cls.token)
        cls.serviceImpl = GenomeFileUtil(cls.cfg)
        cls.dfu = DataFileUtil(os.environ['SDK_CALLBACK_URL'], token=cls.token)
        cls.scratch = cls.cfg['scratch']
        cls.shockURL = cls.cfg['shock-url']
        cls.gfu_cfg = SDKConfig(cls.cfg)
        cls.wsName = "Athaliana_Test"
        cls.prepare_data()

#    @classmethod
#    def tearDownClass(cls):
#        if hasattr(cls, 'wsName'):
#            cls.wsClient.delete_workspace({'workspace': cls.wsName})
#            print('Test workspace was deleted')

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        if hasattr(self.__class__, 'wsName'):
            return self.__class__.wsName
        suffix = int(time.time() * 1000)
        wsName = "test_GenomeFileUtil_" + str(suffix)
        ret = self.getWsClient().create_workspace({'workspace': wsName})  # noqa
        self.__class__.wsName = wsName
        return wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    @classmethod
    def prepare_data(cls):
#        cls.wsClient.create_workspace({'workspace': cls.wsName})
        cls.importer = FastaGFFToGenome(cls.gfu_cfg)

        cls.dtn_phytozome_root = "/kb/module/data/PhytozomeV11_[PhytozomeV11]/"

        cls.gff_filename = 'Test_v1.0.gene.gff3.gz'
        cls.gff_path = os.path.join(cls.scratch, cls.gff_filename)
        shutil.copy(os.path.join("data", "fasta_gff", "JGI", "Plant_Data", cls.gff_filename), cls.gff_path)

        cls.fa_filename = 'Test_v1.0.fa.gz'
        cls.fa_path = os.path.join(cls.scratch, cls.fa_filename)
        shutil.copy(os.path.join("data", "fasta_gff", "JGI", "Plant_Data", cls.fa_filename), cls.fa_path)

    def check_minimal_items_exist(self, result):

        self.assertTrue('genome_info' in result)
        self.assertTrue('genome_ref' in result)
        self.assertTrue('report_name' in result)
        self.assertTrue('report_ref' in result)

        genome_info = result['genome_info']
        self.assertEquals(genome_info[10]['Domain'], 'Eukaryota')
        self.assertEquals(genome_info[10]['Genetic code'], '1')
        self.assertEquals(genome_info[10]['Name'], 'Arabidopsis thaliana')
        self.assertEquals(genome_info[10]['Source'], 'JGI')
        self.assertTrue('GC content' in genome_info[10])
        self.assertTrue(re.match("^\d+?\.\d+?$", genome_info[10]['GC content']) is not None)
        self.assertTrue('Number features' in genome_info[10])
        self.assertTrue(genome_info[10]['Number features'].isdigit())
        self.assertTrue('Size' in genome_info[10])
        self.assertTrue(genome_info[10]['Size'].isdigit())

    def test_athaliana_fasta_gff_to_genome(self):
        species = "Athaliana"
        gene_model = "TAIR10"

        species_root = os.path.join(self.dtn_phytozome_root,species)

        fa_path = os.path.join(species_root,"assembly/Athaliana_167_TAIR9.fa.gz");
        gff_path = os.path.join(species_root,"annotation/Athaliana_167_TAIR10.gene.gff3.gz")

        input_params = {
            'fasta_file': {'path': fa_path},
            'gff_file': {'path': gff_path},
            'genome_name': 'Athaliana',
            'workspace_name': self.getWsName(),
            'source': 'JGI',
            'type': 'Reference',
            'scientific_name': 'Arabidopsis thaliana'
        }

        result = self.getImpl().fasta_gff_to_genome(self.getContext(), input_params)[0]
        self.check_minimal_items_exist(result)
