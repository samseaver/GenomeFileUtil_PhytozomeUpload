import os
import time
import unittest
import json
import shutil
from configparser import ConfigParser

from GenomeFileUtil.GenomeFileUtilImpl import GenomeFileUtil
from GenomeFileUtil.GenomeFileUtilServer import MethodContext
from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.WorkspaceClient import Workspace as workspaceService


class UtilTest(unittest.TestCase):
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
        suffix = int(time.time() * 1000)
        cls.wsName = "test_GenomeFileUtil_" + str(suffix)
        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})
        #cls.genome_ref = 'KBaseExampleData/Escherichia_coli_K-12_MG1655'
        cls.genome_ref, cls.assembly_ref = cls.load_genome_direct(
            'data/e_coli/new_ecoli_genome.json', 'data/e_coli/e_coli_assembly.fasta', 'TestGenome')

    @classmethod
    def load_genome_direct(cls, filename, assembly_filename, obj_name):
        au = AssemblyUtil(os.environ['SDK_CALLBACK_URL'])
        assembly_path = os.path.join(cls.cfg['scratch'], os.path.basename(assembly_filename))
        shutil.copy(assembly_filename, assembly_path)
        assembly_ref = au.save_assembly_from_fasta({
            'workspace_name': cls.wsName,
            'assembly_name': obj_name + '.assembly',
            'file': {'path': assembly_path}
        })

        data = json.load(open(filename))
        data['assembly_ref'] = assembly_ref
        save_info = {
            'workspace': cls.wsName,
            'objects': [{
                'data': data,
                'name': obj_name + '.genome',
                'type': 'KBaseGenomes.Genome',
            }],
        }
        info = cls.wsClient.save_objects(save_info)[0]
        ref = f"{info[6]}/{info[0]}/{info[4]}"
        print('created test genome: ' + ref + ' from file ' + filename)
        return ref, assembly_ref

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def test_genome_protein_to_fasta(self):
        ret = self.serviceImpl.genome_proteins_to_fasta(self.ctx,
                                                        {'genome_ref': self.genome_ref,
                                                         'include_functions': False,
                                                         'include_aliases': False,
                                                         })[0]
        self.assertIn('file_path', ret)
        with open(ret['file_path']) as infile:
            header1 = '>b0001_CDS_1\n'
            self.assertEqual(infile.readline(), header1)
            self.assertEqual(infile.readline(), 'MKRISTTITTTITITTGNGAG\n')
            infile.readline()
            self.assertEqual(len(infile.readline()), 71)

    def test_genome_features_to_fasta(self):
        ret = self.serviceImpl.genome_features_to_fasta(self.ctx,
                                                        {'genome_ref': self.genome_ref,
                                                         'filter_ids': ['b0001', 'b0001_CDS_1'],
                                                         'feature_lists': ['features', 'cdss']
                                                         })[0]
        self.assertIn('file_path', ret)
        with open(ret['file_path']) as infile:
            header1 = infile.readline()
            self.assertIn('>b0001 functions=leader,Amino acid biosynthesis:', header1)
            self.assertIn('GeneID:944742', header1)
            self.assertIn('thrL', header1)
            self.assertEqual(infile.readline(), 'ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGA\n')
            self.assertIn('>b0001_CDS_1', infile.readline())

    def test_bad_inputs(self):
        with self.assertRaisesRegexp(ValueError, 'required field "genome_ref"'):
            ret = self.serviceImpl.genome_features_to_fasta(self.ctx, {})
        with self.assertRaisesRegexp(ValueError, 'Unknown parameter'):
            ret = self.serviceImpl.genome_features_to_fasta(self.ctx,
                                                            {'genome_ref': self.genome_ref,
                                                             'foo': 'bar'})
        with self.assertRaisesRegexp(ValueError, 'Unknown feature_lists'):
            ret = self.serviceImpl.genome_features_to_fasta(self.ctx,
                                                            {'genome_ref': self.genome_ref,
                                                             'feature_lists': ['foo']})
        with self.assertRaisesRegexp(ValueError, 'Object is not a Genome'):
            ret = self.serviceImpl.genome_features_to_fasta(self.ctx, {
                'genome_ref': self.assembly_ref})

