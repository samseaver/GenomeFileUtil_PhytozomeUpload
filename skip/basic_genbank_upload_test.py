import unittest
import os
import time
import shutil
import json

from os import environ
try:
    from ConfigParser import ConfigParser  # py2
except:
    from configparser import ConfigParser  # py3

from pprint import pprint

from Workspace.WorkspaceClient import Workspace as workspaceService
from GenomeFileUtil.GenomeFileUtilImpl import GenomeFileUtil
from GenomeFileUtil.GenomeFileUtilServer import MethodContext

from DataFileUtil.DataFileUtilClient import DataFileUtil

from GenomeFileUtil.JsonIOHelper import (download_genome_to_json_files, 
                                         compare_genome_json_files)


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
        cls.wsClient = workspaceService(cls.wsURL, token=token)
        cls.serviceImpl = GenomeFileUtil(cls.cfg)

        cls.TEST_ECOLI_FILE_FTP = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Escherichia_coli/reference/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gbff.gz'

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        if hasattr(self.__class__, 'wsName'):
            return self.__class__.wsName
        suffix = int(time.time() * 1000)
        wsName = "test_GenomeFileUtil_" + str(suffix)
        ret = self.getWsClient().create_workspace({'workspace': wsName})
        self.__class__.wsName = wsName
        return wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    def test_simple_local_upload(self):
        # fetch the test files and set things up
        gbk_path = "data/e_coli/GCF_000005845.2_ASM584v2_genomic.gbff"

        ### Test for a Local Function Call
        print('attempting upload via local function directly')
        ws_obj_name = 'MyGenome'
        result = self.getImpl().genbank_to_genome(self.getContext(),
            {
                'file': {'path': gbk_path },
                'workspace_name': self.getWsName(),
                'genome_name': ws_obj_name
            })[0]
        pprint(result)
        self.assertIsNotNone(result['genome_ref'])
        target_dir = os.path.join("/kb/module/work/tmp", "GCF_000005845")
        download_genome_to_json_files(self.getContext()['token'], result['genome_ref'],
                                      target_dir)
        #self.assertEqual(0, len(compare_genome_json_files(target_dir, 
        #                                                  os.path.join("/kb/module/test/data", 
        #                                                               "GCF_000005845"))))
        # todo: add test that result is correct

    def test_simple_shock_upload(self):
        ### Test for upload from SHOCK - upload the file to shock first
        print('attempting upload through shock')
        gbk_path = "data/e_coli/GCF_000005845.2_ASM584v2_genomic.gbff"
        data_file_cli = DataFileUtil(os.environ['SDK_CALLBACK_URL'])
        shutil.copy(gbk_path, self.__class__.cfg['scratch'])
        shock_id = data_file_cli.file_to_shock({
            'file_path': os.path.join(self.__class__.cfg['scratch'], gbk_path.split("/")[-1])
        })['shock_id']
        print("Running test")
        ws_obj_name2 = 'MyGenome.2'
        result = self.getImpl().genbank_to_genome(self.getContext(), {
                'file': {'shock_id': shock_id},
                'workspace_name': self.getWsName(),
                'genome_name': ws_obj_name2,
            })[0]
        pprint(result)
        self.assertIsNotNone(result['genome_ref'])
        # todo: add test that result is correct

    def test_simple_ftp_upload(self):
        ### Test for upload via FTP- use something from genbank
        print('attempting upload through ftp url')
        ws_obj_name3 = 'MyGenome.3'
        result = self.getImpl().genbank_to_genome(self.getContext(), {
                'file': {'ftp_url': self.__class__.TEST_ECOLI_FILE_FTP},
                'workspace_name': self.getWsName(),
                'genome_name': ws_obj_name3,
            })[0]
        pprint(result)
        self.assertIsNotNone(result['genome_ref'])



