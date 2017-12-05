import difflib
import time
from os import environ
import unittest
import json
try:
    from ConfigParser import ConfigParser  # py2
except:
    from configparser import ConfigParser  # py3

from pprint import pprint

from Workspace.WorkspaceClient import Workspace as workspaceService
from GenomeFileUtil.GenomeFileUtilImpl import GenomeFileUtil
from GenomeFileUtil.GenomeFileUtilServer import MethodContext


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

    def test_genbank_roundtrip(self):
        # fetch the test files and set things up
        genomeFileUtil = self.getImpl()
        gbk_path = "data/e_coli/GCF_000005845.2_ASM584v2_altered_genomic.gbff"

        ### Test for a Local Function Call
        print('attempting upload via local function directly')
        ws_obj_name = 'TestEcoliAltered'
        result = genomeFileUtil.genbank_to_genome(
            self.getContext(),
            {
                'file': {'path': gbk_path},
                'workspace_name': self.getWsName(),
                'genome_name': ws_obj_name
            })[0]
        pprint(result)
        self.assertIsNotNone(result['genome_ref'])
        genome_data = self.wsClient.get_objects2({'objects': [
            {'ref': result['genome_ref']}]})['data'][0]['data']
        json.dump(genome_data, open('/kb/module/work/tmp/{}.json'.format(
            genome_data['id']), 'w'))
        print('testing Genbank download by building the file')
        genomeFileUtil.export_genome_as_genbank(
            self.getContext(), {'input_ref': result['genome_ref']})
        old_file_path = "/kb/module/test/data/e_coli/KBase_derived_TestEcoliAltered.gbff"
        new_file_path = "/kb/module/work/tmp/TestEcoliAltered/KBase_derived_TestEcoliAltered.gbff"
        differ = difflib.ndiff(open(old_file_path).readlines(),
                               open(new_file_path).readlines())
        diffs = [x for x in differ if x[0] == "+" or x[0] == "-"]
        if len(diffs) > 4:  # expect date line in each contig differ
            raise AssertionError("Output file has changed {}".format(diffs))
