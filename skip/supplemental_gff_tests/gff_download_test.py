import json
import os
import shutil
import time
import unittest
from configparser import ConfigParser
from os import environ

from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.GenomeAnnotationAPIClient import GenomeAnnotationAPI
from GenomeFileUtil.GenomeFileUtilImpl import GenomeFileUtil
from GenomeFileUtil.GenomeFileUtilServer import MethodContext
from installed_clients.WorkspaceClient import Workspace as workspaceService


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
        cls.gaa = GenomeAnnotationAPI(os.environ['SDK_CALLBACK_URL'])
        cls.serviceImpl = GenomeFileUtil(cls.cfg)

        # create one WS for all tests
        suffix = int(time.time() * 1000)
        wsName = "test_GenomeAnnotationAPI_" + str(suffix)
        ret = cls.ws.create_workspace({'workspace': wsName})
        cls.wsName = wsName

        cls.ara_ref = cls.serviceImpl.genbank_to_genome(
            cls.ctx,
            {
                'file': {
                    'path': 'data/Arabidopsis_gbff/A_thaliana_Ensembl_TAIR10_38_chr4_minus_xref.gbff'},
                'workspace_name': cls.wsName,
                'genome_name': "arab",
                'source': 'Ensembl',
                'generate_ids_if_needed': 1,
            })[0]['genome_ref']

        # preload with reference data
        data = json.load(open('data/rhodobacter.json'))
        # save to ws
        save_info = {
                'workspace': wsName,
                'data': data,
                'name': 'rhodobacter'
            }
        info = cls.gaa.save_one_genome_v1(save_info)['info']
        cls.rhodobacter_ref = str(info[6]) +'/' + str(info[0]) + '/' + str(info[4])
        print(('created rhodobacter test genome: ' + cls.rhodobacter_ref))

        # save new genome
        assembly_file_path = os.path.join(cls.cfg['scratch'],
                                          'e_coli_assembly.fasta')
        shutil.copy('data/e_coli/e_coli_assembly.fasta', assembly_file_path)
        au = AssemblyUtil(os.environ['SDK_CALLBACK_URL'])
        assembly_ref = au.save_assembly_from_fasta({
            'workspace_name': cls.wsName,
            'assembly_name': 'ecoli.assembly',
            'file': {'path': assembly_file_path}
        })
        data = json.load(open('data/e_coli/new_ecoli_genome.json'))
        data['assembly_ref'] = assembly_ref
        # save to ws
        save_info = {
                'workspace': wsName,
                'objects': [{
                    'type': 'KBaseGenomes.Genome',
                    'data': data,
                    'name': 'new_ecoli'
                }]
            }
        result = cls.ws.save_objects(save_info)
        info = result[0]
        cls.ecoli_ref = str(info[6]) +'/' + str(info[0]) + '/' + str(info[4])
        print(('created ecoli test genome: ' + cls.rhodobacter_ref))

        # save a GFF file to shock, preload a genome pointing to it
        dfu = DataFileUtil(os.environ['SDK_CALLBACK_URL'])
        shutil.copy('data/rhodobacter.gtf',cls.cfg['scratch'])
        shock_file = dfu.file_to_shock({
                            'file_path': os.path.join(cls.cfg['scratch'], 'rhodobacter.gtf'),
                            'make_handle': 1
                        })
        data['gff_handle_ref']=shock_file['handle']['hid']

        # save to ws
        save_info['objects'][0]['name'] = 'rhodobacter_with_gff'
        result = cls.ws.save_objects(save_info)
        info = result[0]
        cls.rhodobacter_ref_with_gff = str(info[6]) +'/' + str(info[0]) + '/' + str(info[4])
        print(('created rhodobacter test genome with handle: ' + cls.rhodobacter_ref_with_gff))

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

    def test_old_genome_gff_download(self):
        genomeFileUtil = self.getImpl()
        print('testing GFF download by building the file')
        res = genomeFileUtil.genome_to_gff(
            self.getContext(),
            {'genome_ref': self.rhodobacter_ref})[0]
        self.assertEqual(res['from_cache'], 0)

    def test_old_genome_gtf_download(self):
        genomeFileUtil = self.getImpl()
        print('testing GTF download by building the file')
        res2 = genomeFileUtil.genome_to_gff(
            self.getContext(),
            {'genome_ref': self.rhodobacter_ref,
             'target_dir': '/kb/module/work/tmp/rhodo_gft', 'is_gtf': 1})[0]
        self.assertEqual(res2['from_cache'], 0)

    def test_old_genome_gff_export(self):
        # fetch the test files and set things up
        genomeFileUtil = self.getImpl()
        print('testing GTF export')
        res = genomeFileUtil.export_genome_as_gff(
            self.getContext(), {'input_ref': self.rhodobacter_ref_with_gff})[0]
        assert 'shock_id' in res

    def test_new_genome_gff_download(self):
        # fetch the test files and set things up
        genomeFileUtil = self.getImpl()
        print('testing GFF download by building the file')
        res = genomeFileUtil.genome_to_gff(
            self.getContext(), {'genome_ref': self.ecoli_ref})[0]
        self.assertEqual(res['from_cache'], 0)

    def test_new_genome_gtf_download(self):
        # fetch the test files and set things up
        genomeFileUtil = self.getImpl()
        print('testing GTF download by building the file')
        res = genomeFileUtil.genome_to_gff(
            self.getContext(), {'genome_ref': self.ecoli_ref, 'is_gtf': 1})[0]
        self.assertEqual(res['from_cache'], 0)

    def test_ara_gff_download(self):
        # fetch the test files and set things up
        genomeFileUtil = self.getImpl()
        print('testing GFF download by building the file')
        res = genomeFileUtil.genome_to_gff(
            self.getContext(), {'genome_ref': self.ara_ref})[0]
        self.assertEqual(res['from_cache'], 0)

    def test_ara_gtf_download(self):
        # fetch the test files and set things up
        genomeFileUtil = self.getImpl()
        print('testing GTF download by building the file')
        res = genomeFileUtil.genome_to_gff(
            self.getContext(), {'genome_ref': self.ara_ref, 'is_gtf': 1})[0]
        self.assertEqual(res['from_cache'], 0)