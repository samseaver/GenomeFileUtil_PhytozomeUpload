import unittest
import os
import json
import time
import shutil

from os import environ
try:
    from ConfigParser import ConfigParser  # py2
except:
    from configparser import ConfigParser  # py3

from pprint import pprint

from Workspace.WorkspaceClient import Workspace as workspaceService
from GenomeAnnotationAPI.GenomeAnnotationAPIClient import GenomeAnnotationAPI
from GenomeFileUtil.GenomeFileUtilImpl import GenomeFileUtil
from GenomeFileUtil.GenomeFileUtilServer import MethodContext

from DataFileUtil.DataFileUtilClient import DataFileUtil
from AssemblyUtil.AssemblyUtilClient import AssemblyUtil


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
        cls.gaapi = GenomeAnnotationAPI(os.environ['SDK_CALLBACK_URL'])
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('GenomeFileUtil'):
            cls.cfg[nameval[0]] = nameval[1]
        cls.wsURL = cls.cfg['workspace-url']
        cls.ws = workspaceService(cls.wsURL, token=token)
        cls.serviceImpl = GenomeFileUtil(cls.cfg)

        # create one WS for all tests
        suffix = int(time.time() * 1000)
        wsName = "test_GenomeAnnotationAPI_" + str(suffix)
        ret = cls.ws.create_workspace({'workspace': wsName})
        cls.wsName = wsName

        

       

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.ws.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def getWsClient(self):
        return self.__class__.ws

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


    def load_genome_direct(self, filename, assembly_filename, obj_name):
        au = AssemblyUtil(os.environ['SDK_CALLBACK_URL'])
        assembly_ref = au.save_assembly_from_fasta({
            'workspace_name': self.getWsName(),
            'assembly_name': obj_name + '.assembly',
            'file': {'path': assembly_filename}
        })
        pprint('created test assembly: ' + assembly_ref)

        with open(filename, 'r') as file:
            data_str = file.read()
        data = json.loads(data_str)
        data['assembly_ref'] = assembly_ref
        # save to ws
        save_info = {
                'workspace': self.getWsName(),
                'data': data,
                'name': obj_name + '.genome'
            }
        result = self.gaapi.save_one_genome_v1(save_info)
        info = result['info']
        ref = str(info[6]) + '/' + str(info[0]) + '/' + str(info[4])
        print('created test genome: ' + ref + ' from file ' + filename)
        return ref

    def load_genome_with_cache(self, filename, gbff_cache_filename):
        """ cache filename needs to in scratch space """
        with open(filename, 'r') as file:
            data_str = file.read()
        data = json.loads(data_str)
        # save to ws
        save_info = {
                'workspace': self.getWsName(),
                'objects': [{
                    'type': 'KBaseGenomes.Genome',
                    'data': data,
                    'name': 'e_coli'
                }]
            }
        dfu = DataFileUtil(os.environ['SDK_CALLBACK_URL'])
        shock_file = dfu.file_to_shock({
                            'file_path': gbff_cache_filename,
                            'make_handle': 1
                        })
        data['genbank_handle_ref'] = shock_file['handle']['hid']
        # save to ws
        save_info['objects'][0]['name'] = 'e_coli_with_genbank'
        result = self.ws.save_objects(save_info)
        info = result[0]
        ref = str(info[6]) + '/' + str(info[0]) + '/' + str(info[4])
        print('created test genome with gbff cache: ' + ref + ' from file ' + filename)
        return ref

    def test_simple_genbank_download(self):
        # load test data data
        assembly_file_path = os.path.join(self.cfg['scratch'], 'e_coli_assembly.fasta')
        shutil.copy('data/e_coli/e_coli_assembly.fasta', assembly_file_path)
        e_coli_ref = self.load_genome_direct('data/e_coli/e_coli.json',
                                             assembly_file_path, 'tax_bug_test')

        # run the test
        genomeFileUtil = self.getImpl()
        print('testing Genbank download by building the file')
        res1 = genomeFileUtil.genome_to_genbank(self.getContext(),
            {'genome_ref': e_coli_ref})[0]
        self.assertEqual(res1['from_cache'], 0)

    def test_check_for_taxonomy_bug(self):
        # load test data data
        assembly_file_path = os.path.join(self.cfg['scratch'], 'e_coli_assembly.fasta')
        shutil.copy('data/e_coli/e_coli_assembly.fasta', assembly_file_path)
        e_coli_ref = self.load_genome_direct(
            'data/taxonomy_bug_test_genome.json', assembly_file_path, 'tax_bug_test')
        # run the test
        genomeFileUtil = self.getImpl()
        print('testing Genbank download by building the file')
        res1 = genomeFileUtil.genome_to_genbank(self.getContext(),
            {'genome_ref': e_coli_ref})[0]
        self.assertEqual(res1['from_cache'], 0)



