import os
import shutil
import time
import unittest
import mock
from configparser import ConfigParser

from DataFileUtil.DataFileUtilClient import DataFileUtil
from GenomeFileUtil.GenomeFileUtilImpl import GenomeFileUtil
from GenomeFileUtil.GenomeFileUtilServer import MethodContext
from GenomeFileUtil.core.GenomeInterface import GenomeInterface
from Workspace.WorkspaceClient import Workspace as workspaceService

class GenomeFileUtilTest(unittest.TestCase):
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
        config.read(config_file)
        for nameval in config.items('GenomeFileUtil'):
            cls.cfg[nameval[0]] = nameval[1]
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = workspaceService(cls.wsURL, token=token)
        cls.serviceImpl = GenomeFileUtil(cls.cfg)
        cls.token = token

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
        self.getWsClient().create_workspace({'workspace': wsName})
        self.__class__.wsName = wsName
        return wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    def test_full_sequence(self):
        # features should not have sequences in it. But both non_coding_features and CDSs should have sequences.
        print("test_full_sequence")
        gbk_path = "data/e_coli/GCF_000005845.2_ASM584v2_genomic.gbff"
        ws_obj_name = 'full_sequence'
        result = self.getImpl().genbank_to_genome(
            self.getContext(),
            {
              'file': {
                  'path': gbk_path},
              'workspace_name': self.getWsName(),
              'genome_name': ws_obj_name,
              'generate_ids_if_needed': 1
            })[0]
        data_file_cli = DataFileUtil(os.environ['SDK_CALLBACK_URL'], 
                                    token=self.__class__.token,
                                    service_ver='dev')
        genome = data_file_cli.get_objects({'object_refs': [result['genome_ref']]})['data'][0]['data']
        count_features_without_dna_sequence = 0
        for feature in genome['features']:
            if "dna_sequence" not in feature:
                count_features_without_dna_sequence += 1
        count_non_coding_features_without_sequence = 0
        for feature in genome['non_coding_features']:
            if "dna_sequence" not in feature:
                if feature["dna_sequence_length"] <= 10000:
                    count_non_coding_features_without_sequence += 1
                    print("non_coding_feature_without_sequence: " + str(feature))
        count_cdss_without_sequence = 0
        for feature in genome['cdss']:
            if "dna_sequence" not in feature:
                count_cdss_without_sequence += 1        
        self.assertTrue(count_features_without_dna_sequence == 0,"All features should have DNA sequences.")
        self.assertTrue(count_non_coding_features_without_sequence == 0,
                        "All non_coding_features should have DNA sequences.")
        self.assertTrue(count_cdss_without_sequence == 0,"All CDSs should have DNA sequences.")

    @mock.patch("GenomeFileUtil.core.GenomeInterface.MAX_GENOME_SIZE", 14000000)
    def test_partial_sequence(self):
        # features should not have sequences in it. But both non_coding_features and CDSs should have sequences.
        print("test_partial_sequence")
        gbk_path = "data/e_coli/GCF_000005845.2_ASM584v2_genomic.gbff"
        ws_obj_name = 'partial_sequence'
        result = self.getImpl().genbank_to_genome(
            self.getContext(),
            {
              'file': {
                  'path': gbk_path},
              'workspace_name': self.getWsName(),
              'genome_name': ws_obj_name,
              'generate_ids_if_needed': 1
            })[0]
        data_file_cli = DataFileUtil(os.environ['SDK_CALLBACK_URL'], 
                                    token=self.__class__.token,
                                    service_ver='dev')
        genome = data_file_cli.get_objects({'object_refs': [result['genome_ref']]})['data'][0]['data']
        count_features_with_dna_sequence = 0
        for feature in genome['features']:
            if "dna_sequence" in feature:
                count_features_with_dna_sequence += 1
        count_non_coding_features_without_sequence = 0
        for feature in genome['non_coding_features']:
            if "dna_sequence" not in feature:
                if feature["dna_sequence_length"] <= 10000:
                    count_non_coding_features_without_sequence += 1
                    print("non_coding_feature_without_sequence: " + str(feature))
        count_cdss_without_sequence = 0
        for feature in genome['cdss']:
            if "dna_sequence" not in feature:
                count_cdss_without_sequence += 1        
        self.assertTrue(count_features_with_dna_sequence == 0,"All features should not have DNA sequences.")
        self.assertTrue(count_non_coding_features_without_sequence == 0,
                        "All non_coding_features should have DNA sequences.")
        self.assertTrue(count_cdss_without_sequence == 0,"All CDSs should have DNA sequences.")

    @mock.patch("GenomeFileUtil.core.GenomeInterface.MAX_GENOME_SIZE", 9000000)
    def test_no_sequence_kept(self):
        # features, cds, and non_coding_features should not have sequences in it.
        print("test_no_sequence_kept")
        gbk_path = "data/e_coli/GCF_000005845.2_ASM584v2_genomic.gbff"
        ws_obj_name = 'no_sequence'
        result = self.getImpl().genbank_to_genome(
            self.getContext(),
            {
              'file': {
                  'path': gbk_path},
              'workspace_name': self.getWsName(),
              'genome_name': ws_obj_name,
              'generate_ids_if_needed': 1
            })[0]
        data_file_cli = DataFileUtil(os.environ['SDK_CALLBACK_URL'], 
                                    token=self.__class__.token,
                                    service_ver='dev')
        genome = data_file_cli.get_objects({'object_refs': [result['genome_ref']]})['data'][0]['data']
        count_features_with_dna_sequence = 0
        for feature in genome['features']:
            if "dna_sequence" in feature:
                count_features_with_dna_sequence += 1
        count_non_coding_features_with_sequence = 0
        for feature in genome['non_coding_features']:
            if "dna_sequence" in feature:
                count_non_coding_features_with_sequence += 1
        count_cdss_with_sequence = 0
        for feature in genome['cdss']:
            if "dna_sequence" in feature:
                count_cdss_with_sequence += 1        
        self.assertTrue(count_features_with_dna_sequence == 0,"All features should not have DNA sequences.")
        self.assertTrue(count_non_coding_features_with_sequence == 0,
                        "All non_coding_features should not have DNA sequences.")
        self.assertTrue(count_cdss_with_sequence == 0,"All CDSs should not have DNA sequences.")

    @mock.patch("GenomeFileUtil.core.GenomeInterface.MAX_GENOME_SIZE", 1)
    def test_max_genome_size(self):
        with self.assertRaisesRegex(ValueError, "This genome size of "):
            GenomeInterface.validate_genome({"taxon_ref": "", "domain": ""})