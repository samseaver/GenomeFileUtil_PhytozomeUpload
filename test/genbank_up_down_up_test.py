import unittest
import time
import os
import shutil

try:
    from ConfigParser import ConfigParser  # py2
except:
    from configparser import ConfigParser  # py3

from Workspace.WorkspaceClient import Workspace as workspaceService
from GenomeFileUtil.GenomeFileUtilImpl import GenomeFileUtil
from GenomeFileUtil.GenomeFileUtilServer import MethodContext
from DataFileUtil.DataFileUtilClient import DataFileUtil
from pprint import pprint


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
        gbk_path = "data/e_coli/GCF_000005845.2_ASM584v2_altered_genomic.gbff"
        ws_obj_name = 'ecoli_2contigs_orig'
        suffix = int(time.time() * 1000)
        cls.wsName = "test_GenomeFileUtil_" + str(suffix)
        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})
        result = cls.serviceImpl.genbank_to_genome(
            cls.ctx,
            {
              'file': {
                  'path': gbk_path},
              'workspace_name': cls.wsName,
              'genome_name': ws_obj_name,
              'generate_ids_if_needed': 1
            })[0]
#        print("HERE IS THE RESULT:")
        data_file_cli = DataFileUtil(os.environ['SDK_CALLBACK_URL'], 
                                token=cls.ctx['token'],
                                service_ver='dev')
        cls.genome_orig = data_file_cli.get_objects({'object_refs': [result['genome_ref']]})['data'][0]['data']
#        genome_data = self.wsClient.get_objects2({'objects': [
#            {'ref': result['genome_ref']}]})['data'][0]['data']
#        json.dump(genome_data, open('/kb/module/work/tmp/{}.json'.format(
#            genome_data['id']), 'w'))
        print('testing Genbank download by building the file')
        #genomeFileUtil.export_genome_as_genbank(cls.ctx,
        cls.serviceImpl.export_genome_as_genbank(cls.ctx,
                                {'input_ref': result['genome_ref']})
#        old_file_path = "/kb/module/test/data/e_coli/KBase_derived_TestEcoliAltered.gbff"
        new_gbk_path = "/kb/module/work/tmp/ecoli_2contigs_orig/KBase_derived_ecoli_2contigs_orig.gbff"
        new_ws_obj_name = 'ecoli_2contigs_new'
        new_result = cls.serviceImpl.genbank_to_genome(
            cls.ctx,
            {
              'file': {
                  'path': new_gbk_path},
              'workspace_name': cls.wsName,
              'genome_name': new_ws_obj_name,
              'generate_ids_if_needed': 1
            })[0]
        cls.genome_new = data_file_cli.get_objects({'object_refs': [new_result['genome_ref']]})['data'][0]['data']

#        print("GENE 1: ")
#        pprint(cls.genome['features'][0])
#        pprint(result)



    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def test_gene_count(self):
        genome_orig = self.__class__.genome_orig
        genome_new = self.__class__.genome_new
        print "Len GO: " + str(len(genome_orig['features']))
        print "Len GN: " + str(len(genome_new['features']))
        self.assertTrue(len(genome_orig['features']) == len(genome_new['features']))
        self.maxDiff = None
        self.assertEqual(genome_orig, genome_new)





