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

    def getTempGenbank(self):
        tmp_dir = self.__class__.cfg['scratch']
        file_name = "GCF_000005845.2_ASM584v2_genomic.gbff.gz"
        gbk_path = os.path.join(tmp_dir, file_name)
        if not os.path.exists(gbk_path):
            ftp_url = self.TEST_ECOLI_FILE_FTP
            print('downloading test data from:' + ftp_url)
            with closing(urllib2.urlopen(ftp_url)) as r:
                with open(gbk_path, 'wb') as f:
                    shutil.copyfileobj(r, f)
        return gbk_path

    def test_simple_upload(self):
        # fetch the test files and set things up
        genomeFileUtil = self.getImpl()
        tmp_dir = self.__class__.cfg['scratch']
        gbk_path = self.getTempGenbank()

        ### Test for a Local Function Call
        print('attempting upload via local function directly')
        ws_obj_name = 'MyGenome'
        result = genomeFileUtil.genbank_to_genome(self.getContext(), 
            {
                'file' : { 'path':gbk_path },
                'workspace_name':self.getWsName(),
                'genome_name':ws_obj_name
            });
        pprint(result)
        self.assertIsNotNone(result['genome_ref'])
        # todo: add test that result is correct

        ### Test for upload from SHOCK - upload the file to shock first
        print('attempting upload through shock')
        data_file_cli = DataFileUtil(os.environ['SDK_CALLBACK_URL'], 
                                token=self.__class__.ctx['token'],
                                service_ver='dev')
        shock_id = data_file_cli.file_to_shock({'file_path': gbk_path})['shock_id']
        ws_obj_name2 = 'MyGenome.2'
        result2 = genomeFileUtil.genbank_to_genome(self.getContext(), 
            {
                'file': {'shock_id':shock_id},
                'workspace_name':self.getWsName(),
                'genome_name':ws_obj_name2,
            });
        pprint(result2)
        self.assertIsNotNone(result['genome_ref'])
        # todo: add test that result is correct

        ### Test for upload via FTP- use something from genbank
        print('attempting upload through ftp url')
        ws_obj_name3 = 'MyGenome.3'
        result3 = genomeFileUtil.genbank_to_genome(self.getContext(), 
            {
                'file':{'ftp_url':'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Escherichia_coli/reference/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gbff.gz'},
                'workspace_name':self.getWsName(),
                'genome_name':ws_obj_name3,
            });
        pprint(result3)
        self.assertIsNotNone(result3['genome_ref'])



