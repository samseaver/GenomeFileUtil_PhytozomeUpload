import json
import os
import shutil
import time
import unittest
from os import environ
from pprint import pprint

from configparser import ConfigParser  # py3

from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.GenomeAnnotationAPIClient import GenomeAnnotationAPI
from GenomeFileUtil.GenomeFileUtilImpl import GenomeFileUtil
from GenomeFileUtil.GenomeFileUtilServer import MethodContext
from installed_clients.WorkspaceClient import Workspace as workspaceService


class GenomeFileUtilTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        print('xyzxyzxyz')
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
        cls.gaa = GenomeAnnotationAPI(os.environ['SDK_CALLBACK_URL'])
        cls.serviceImpl = GenomeFileUtil(cls.cfg)

        # create one WS for all tests
        suffix = int(time.time() * 1000)
        wsName = "test_GenomeAnnotationAPI_" + str(suffix)
        cls.ws.create_workspace({'workspace': wsName})
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
        self.getWsClient().create_workspace({'workspace': wsName})
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
        save_info = {
            'workspace': self.getWsName(),
            'data': data,
            'name': obj_name + '.genome'
        }
        info = self.gaa.save_one_genome_v1(save_info)['info']
        ref = "{}/{}/{}".format(info[6], info[0], info[4])
        print('created test genome: ' + ref + ' from file ' + filename)
        return ref

    def test_simple_genbank_download(self):
        # load test data data
        assembly_file_path = os.path.join(self.cfg['scratch'], 'e_coli_assembly.fasta')
        shutil.copy('data/e_coli/e_coli_assembly.fasta', assembly_file_path)
        e_coli_ref = self.load_genome_direct(
            'data/e_coli/e_coli.json',
            assembly_file_path,
            'tax_bug_test'
        )
        # run the test
        genomeFileUtil = self.getImpl()
        print('testing Genbank download by building the file')
        res1 = genomeFileUtil.genome_to_genbank(self.getContext(), {
            'genome_ref': e_coli_ref
        })[0]
        self.assertEqual(res1['from_cache'], 0)

    def test_check_for_taxonomy_bug(self):
        # load test data data
        assembly_file_path = os.path.join(self.cfg['scratch'], 'e_coli_assembly.fasta')
        shutil.copy('data/e_coli/e_coli_assembly.fasta', assembly_file_path)
        e_coli_ref = self.load_genome_direct(
            'data/taxonomy_bug_test_genome.json',
            assembly_file_path,
            'tax_bug_test'
        )
        # run the test
        genomeFileUtil = self.getImpl()
        print('testing Genbank download by building the file')
        res1 = genomeFileUtil.genome_to_genbank(self.getContext(), {'genome_ref': e_coli_ref})[0]
        self.assertEqual(res1['from_cache'], 0)

    def test_contig_set_download(self):
        contig_data = json.load(open('data/rhodobacter_contigs.json'))
        save_info = {
            'workspace': self.getWsName(),
            'objects': [{
                'type': 'KBaseGenomes.ContigSet',
                'data': contig_data,
                'name': 'rhodobacter_contigs'
            }]
        }
        info = self.ws.save_objects(save_info)[0]
        cs_ref = "{}/{}/{}".format(info[6], info[0], info[4])
        genome_data = json.load(open('data/rhodobacter.json'))
        genome_data['contigset_ref'] = cs_ref
        save_info = {
            'workspace': self.getWsName(),
            'data': genome_data,
            'name': 'rhodobacter_genome'
        }
        info = self.gaa.save_one_genome_v1(save_info)['info']
        ref = "{}/{}/{}".format(info[6], info[0], info[4])
        # run the test
        genomeFileUtil = self.getImpl()
        print('testing Genbank download by building the file')
        res1 = genomeFileUtil.genome_to_genbank(self.getContext(), {'genome_ref': ref})[0]
        self.assertEqual(res1['from_cache'], 0)
