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

from biokbase.workspace.client import Workspace as workspaceService
from GenomeFileUtil.GenomeFileUtilImpl import GenomeFileUtil
from GenomeFileUtil.GenomeFileUtilServer import MethodContext

from DataFileUtil.DataFileUtilClient import DataFileUtil


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

        # create one WS for all tests
        suffix = int(time.time() * 1000)
        wsName = "test_GenomeAnnotationAPI_" + str(suffix)
        ret = cls.ws.create_workspace({'workspace': wsName})
        cls.wsName = wsName

        # preload with reference data
        with open ('data/e_coli.json', 'r') as file:
            data_str=file.read()
        data = json.loads(data_str)
        # save to ws
        save_info = {
                'workspace':wsName,
                'objects': [{
                    'type':'KBaseGenomes.Genome',
                    'data':data,
                    'name':'e_coli'
                }]
            }
        result = cls.ws.save_objects(save_info)
        info = result[0]
        cls.e_coli_ref = str(info[6]) +'/' + str(info[0]) + '/' + str(info[4])
        print('created e.coli test genome: ' + cls.e_coli_ref)

        # save a genbank file to shock, preload a genome pointing to it
        dfu = DataFileUtil(os.environ['SDK_CALLBACK_URL'])
        shutil.copy('data/GCF_000005845.2_ASM584v2_genomic.gbff', cls.cfg['scratch'])
        shock_file = dfu.file_to_shock({
                            'file_path': os.path.join(cls.cfg['scratch'], 'GCF_000005845.2_ASM584v2_genomic.gbff'),
                            'make_handle': 1
                        })
        data['genbank_handle_ref'] = shock_file['handle']['hid']

        # save to ws
        save_info['objects'][0]['name'] = 'e_coli_with_genbank'
        result = cls.ws.save_objects(save_info)
        info = result[0]
        cls.e_coli_ref_with_genbank = str(info[6]) +'/' + str(info[0]) + '/' + str(info[4])
        print('created e.coli test genome with handle: ' + cls.e_coli_ref_with_genbank)

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

    def test_simple_genbank_download(self):
        genomeFileUtil = self.getImpl()
        print('testing Genbank download by building the file')
        res1 = genomeFileUtil.genome_to_genbank(self.getContext(),
            {'genome_ref': self.e_coli_ref})[0]
        self.assertEqual(res1['from_cache'], 0)

    def test_simple_genbank_download_from_cache(self):
        genomeFileUtil = self.getImpl()
        print('testing Genbank download from cached file')
        res2 = genomeFileUtil.genome_to_genbank(self.getContext(),
            {'genome_ref': self.e_coli_ref_with_genbank})[0]
        self.assertEqual(res2['from_cache'], 1)
        self.assertTrue('genbank_file' in res2)
        self.assertTrue('file_path' in res2['genbank_file'])
        self.assertTrue(bool(res2['genbank_file']['file_path']))

    def test_simple_genbank_export_from_cache(self):
        genomeFileUtil = self.getImpl()
        print('testing Genbank export from cached file')
        res = genomeFileUtil.export_genome_as_genbank(self.getContext(),
            {'input_ref': self.e_coli_ref_with_genbank})[0]
        self.assertTrue('shock_id' in res)
        self.assertTrue(bool(res['shock_id']))


