import os
import time
import unittest
from configparser import ConfigParser

from installed_clients.DataFileUtilClient import DataFileUtil
from GenomeFileUtil.GenomeFileUtilImpl import GenomeFileUtil
from GenomeFileUtil.GenomeFileUtilServer import MethodContext
from installed_clients.WorkspaceClient import Workspace as workspaceService


@unittest.skip('x')
class GenomeFileUtilTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        print('xyzxyzxyz')
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
        gff_path = "data/fasta_gff/RefSeq/Bacterial_Data/NC_021490.gff.gz"
        fasta_path = "data/fasta_gff/RefSeq/Bacterial_Data/NC_021490.fasta.gz"
        suffix = int(time.time() * 1000)
        cls.wsName = "test_GenomeFileUtil_" + str(suffix)
        cls.wsClient.create_workspace({'workspace': cls.wsName})

        print('Uploading GFF file')
        result = cls.serviceImpl.fasta_gff_to_genome(cls.ctx, {
            'workspace_name': cls.wsName,
            'genome_name': 'MyGenome',
            'taxon_id': '3702',
            'fasta_file': {'path': fasta_path},
            'gff_file': {'path': gff_path},
            'source': 'GFF',
            'type': 'Reference'
        })[0]
        data_file_cli = DataFileUtil(os.environ['SDK_CALLBACK_URL'])
        cls.genome_orig = data_file_cli.get_objects(
            {'object_refs': [result['genome_ref']]})['data'][0]['data']

        print('testing GFF download by building the file')
        down_result = cls.serviceImpl.genome_to_gff(
            cls.ctx, {'genome_ref': result['genome_ref']})[0]

        print('Reuploading GFF file')
        new_result = cls.serviceImpl.fasta_gff_to_genome(
            cls.ctx,
            {
                'workspace_name': cls.wsName,
                'genome_name': 'MyGenome',
                'fasta_file': {'path': fasta_path},
                'gff_file': {'path': down_result['file_path']},
                'source': 'GFF',
                'type': 'Reference'
            })[0]
        cls.genome_new = data_file_cli.get_objects({'object_refs': [new_result['genome_ref']]})['data'][0]['data']

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def feature_list_comparison(self, feature_list_name):
        genome_orig = self.genome_orig
        genome_new = self.genome_new
        self.assertTrue(len(genome_orig[feature_list_name]) == len(genome_new[feature_list_name]),
                        feature_list_name + " list is not of equal length in Original and New Genomes.")
        print("\n\n" + feature_list_name + " TOTAL NUMBER:" + str(len(genome_orig[feature_list_name])))
        orig_dict = dict([(x['id'], x) for x in genome_orig[feature_list_name]])
        new_dict = dict([(x['id'], x) for x in genome_new[feature_list_name]])

        first_pass_matches = 0
        first_pass_non_match = 0
        second_pass_matches = 0

        for key in orig_dict:
            orig_feature = orig_dict[key]
            new_feature = new_dict[key]
            if "aliases" in orig_feature:
                orig_feature['aliases'] = sorted(orig_feature.get('aliases', []))
                new_feature['aliases'] = sorted(new_feature.get('aliases', []))
            if "db_xrefs" in orig_feature:
                orig_feature['db_xrefs'] = sorted(orig_feature.get('db_xrefs', []))
                new_feature['db_xrefs'] = sorted(new_feature.get('db_xrefs', []))
            if "functions" in orig_feature:
                orig_feature["functions"] = sorted(orig_feature.get('functions', []))
                new_feature["functions"] = sorted(new_feature.get('functions', []))
            if orig_feature == new_feature:
                first_pass_matches += 1
            else:
                first_pass_non_match += 1
                orig_feature.pop("note", None)
                new_feature.pop("note", None)
                orig_feature.pop('inference_data', None)
                new_feature.pop('inference_data', None)
                if "warnings" in orig_feature and "warnings" not in new_feature:
                    del(orig_feature["warnings"])
                # THESE ARE TEMPORARY TO FIND OTHER ISSUES:
                # if feature_list_name == "features":
                #     if 'protein_translation' in orig_feature:
                #         del orig_feature['protein_translation']
                #         del new_feature['protein_translation']
                # if 'protein_translation_length' in orig_feature:
                #     del orig_feature['protein_translation_length']
                #     del new_feature['protein_translation_length']
                # if 'protein_md5' in orig_feature:
                #     del orig_feature['protein_md5']
                #     del new_feature['protein_md5']
                # if 'functions' in orig_feature:
                #     del orig_feature['functions']
                #     del new_feature['functions']
                if orig_feature == new_feature:
                    second_pass_matches += 1
                else:
                    self.maxDiff = None
                    self.assertEqual(orig_feature, new_feature)
        self.assertEqual(len(orig_dict), (first_pass_matches + second_pass_matches),
                         "There were %d first pass matches and %d second pass matches out of %d items in %s" %
                         (first_pass_matches, second_pass_matches, len(orig_dict), feature_list_name))

    def test_gene_features(self):
        self.feature_list_comparison("features")

    def test_cds_features(self):
        self.feature_list_comparison("cdss")

    def test_mrna_features(self):
        self.feature_list_comparison("mrnas")

    def test_ncf_features(self):
        self.feature_list_comparison("non_coding_features")
