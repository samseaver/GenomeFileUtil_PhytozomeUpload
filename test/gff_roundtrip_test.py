import difflib
import time
from os import environ
import unittest
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

    def test_genbank_to_gff_roundtrip(self):
        # fetch the test files and set things up
        genomeFileUtil = self.getImpl()
        gbk_path = "data/drosophila/small_test.tar.gz"

        ### Test for a Local Function Call
        print('attempting upload via local function directly')
        ws_obj_name = 'DrosophilaGenome.1'
        result = genomeFileUtil.genbank_to_genome(
            self.getContext(),
            {
                'file': {'path': gbk_path},
                'workspace_name': self.getWsName(),
                'genome_name': ws_obj_name
            })[0]
        pprint(result)
        self.assertIsNotNone(result['genome_ref'])
        print('testing gff download by building the file')
        down_result = genomeFileUtil.genome_to_gff(
            self.getContext(), {'genome_ref': result['genome_ref']})[0]
        old_file_path = "/kb/module/test/data/drosophila/DrosophilaGenome.1.gff"
        new_file_path = down_result['file_path']
        differ = difflib.ndiff(open(old_file_path).readlines(),
                               open(new_file_path).readlines())
        diffs = [x for x in differ if x[0] == "+" or x[0] == "-"]
        if len(diffs) > 4:  # expect date line in each contig differ
            raise AssertionError("Output file has changed {}".format(diffs))
