import unittest
import os
import json
import time
import shutil
import urllib2
from contextlib import closing

from os import environ
try:
    from ConfigParser import ConfigParser  # py2
except:
    from configparser import ConfigParser  # py3

from pprint import pprint

from biokbase.workspace.client import Workspace as workspaceService
from GenomeFileUtil.GenomeFileUtilImpl import GenomeFileUtil
from GenomeFileUtil.GenomeFileUtilServer import MethodContext

from DataFileUtil.DataFileUtilClient import DataFileUtil


class MinimalGenbankUploadTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        print('setting up class')
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
        cls.impl = GenomeFileUtil(cls.cfg)

        cls.MINIMAL_TEST_FILE = os.path.join( cls.cfg['scratch'], 'minimal.gbff')
        shutil.copy('data/minimal.gbff', cls.MINIMAL_TEST_FILE )

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.ws.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def getWsClient(self):
        return self.ws

    def getWsName(self):
        if hasattr(self.__class__, 'wsName'):
            return self.__class__.wsName
        suffix = int(time.time() * 1000)
        wsName = "test_GenomeFileUtil_" + str(suffix)
        ret = self.getWsClient().create_workspace({'workspace': wsName})
        self.__class__.wsName = wsName
        return wsName

    def getImpl(self):
        return self.__class__.impl

    def getContext(self):
        return self.__class__.ctx

    def test_upload(self):

        # fetch the test files and set things up
        genomeFileUtil = self.getImpl()
        gbk_path = self.MINIMAL_TEST_FILE


        genomeFileUtil = self.getImpl()

        result = genomeFileUtil.genbank_to_genome(self.getContext(),
                                    {
                                        'file':{'path': gbk_path},
                                        'workspace_name': self.getWsName(),
                                        'genome_name': 'something',
                                    })[0]

        pprint(result)




