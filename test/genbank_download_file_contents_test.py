import unittest
import os
import time

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

    def test_simple_download(self):
        # fetch the test files and set things up
        genomeFileUtil = self.getImpl()
        gbk_path = "data/e_coli/GCF_000005845.2_ASM584v2_altered_genomic.gbff"

        ### Test for a Local Function Call
        print('attempting upload via local function directly')
        ws_obj_name = 'TestEcoliAltered'
        result = genomeFileUtil.genbank_to_genome(self.getContext(), 
            {
                'file' : { 'path': gbk_path },
                'workspace_name':self.getWsName(),
                'genome_name':ws_obj_name
            })[0]
        pprint(result)
        self.assertIsNotNone(result['genome_ref'])
        print('testing Genbank download by building the file')
        genomeFileUtil.export_genome_as_genbank(
            self.getContext(), {'input_ref': result['genome_ref']})
        file_path = "/kb/module/work/tmp/TestEcoliAltered/KBase_derived_TestEcoliAltered.gbff"
        gb_file = open(file_path, 'r').read()
        gene_any_count = gb_file.count("gene")  #  4555
        rev_feature_count = gb_file.count("complement")  #  2283
        gene_true_count = gb_file.count("     gene            ")  #  4314
        contig_count = gb_file.count("LOCUS       ")  #  2
        #self.assertEqual(gene_any_count, 4555)
        #self.assertEqual(rev_feature_count, 2283)
        #self.assertEqual(gene_true_count, 4314)
        self.assertEqual(contig_count, 2)




