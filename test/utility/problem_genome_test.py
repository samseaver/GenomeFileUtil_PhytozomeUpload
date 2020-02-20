import os
import time
import unittest
from configparser import ConfigParser

from GenomeFileUtil.GenomeFileUtilImpl import GenomeFileUtil
from GenomeFileUtil.GenomeFileUtilServer import MethodContext
from installed_clients.WorkspaceClient import Workspace as workspaceService


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

    def test_problem_genome_for_json(self):
        # TO RUN THIS TEST:
        # 1) Add your gbff fil to the test/data directory
        # 2) Change the gbk_path below appropriately
        # 3) Move this test into the test directory
        # 4) Remove the comments from the json dump printing command that is just before
        #       the call to 
        #      "result = self.gi.save_one_genome({
        #                    'workspace': params['workspace_name'],
        #                    'name': params['genome_name'],
        #                    'data': genome,
        #                    "meta": params['metadata'],
        #                    })"
        # 5) Run kb-sdk test
        # 6) Look for the json file for the object at test_local/workdir/tmp/ProblemGenome.json
        #gbk_path = "data/YOUR_GENOME_HERE.gbff"
        gbk_path = "data/mouse_chr1.gbff"
        ws_obj_name = 'problem_genome'
        result = self.getImpl().genbank_to_genome(
            self.getContext(),
            {
              'file': {
                  'path': gbk_path},
              'workspace_name': self.getWsName(),
              'genome_name': ws_obj_name,
              'generate_ids_if_needed': 1
            })[0]
        self.assertTrue(int(
            result['genome_info'][10]['Number of Protein Encoding Genes']) > 0)

