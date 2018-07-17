import os
import time
import unittest
from configparser import ConfigParser

from DataFileUtil.DataFileUtilClient import DataFileUtil
from GenomeFileUtil.GenomeFileUtilImpl import GenomeFileUtil
from GenomeFileUtil.GenomeFileUtilServer import MethodContext
from Workspace.WorkspaceClient import Workspace as workspaceService


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
        gbk_path = "data/Arabidopsis_gbff/Arab_Chloro_Modified.gbff"
        ws_obj_name = 'oddities_orig'
        suffix = int(time.time() * 1000)
        cls.wsName = "test_GenomeFileUtil_" + str(suffix)
        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})
        result = cls.serviceImpl.genbank_to_genome(
            cls.ctx,
            {
              'file': {
                  'path': gbk_path},
              'workspace_name': cls.wsName,
              'genome_name': ws_obj_name,
              'generate_ids_if_needed': 1
            })[0]
#        print("HERE IS THE RESULT:")
        data_file_cli = DataFileUtil(os.environ['SDK_CALLBACK_URL'], 
                                token=cls.ctx['token'],
                                service_ver='dev')
        cls.genome_orig = data_file_cli.get_objects({'object_refs': [result['genome_ref']]})['data'][0]['data']
        print('testing Genbank download by building the file')
        #genomeFileUtil.export_genome_as_genbank(cls.ctx,
        cls.serviceImpl.export_genome_as_genbank(cls.ctx,
                                {'input_ref': result['genome_ref']})
        new_gbk_path = "/kb/module/work/tmp/oddities_orig/KBase_derived_oddities_orig.gbff"
        new_ws_obj_name = 'oddities_new'
        new_result = cls.serviceImpl.genbank_to_genome(
            cls.ctx,
            {
              'file': {
                  'path': new_gbk_path},
              'workspace_name': cls.wsName,
              'genome_name': new_ws_obj_name,
              'generate_ids_if_needed': 1
            })[0]
        cls.genome_new = data_file_cli.get_objects({'object_refs': [new_result['genome_ref']]})['data'][0]['data']


    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def feature_list_comparison(self, feature_list_name):
        genome_orig = self.__class__.genome_orig
        genome_new = self.__class__.genome_new
        i = 0 #counter loop over feature
        self.assertTrue(len(genome_orig[feature_list_name]) == len(genome_new[feature_list_name]),
                    feature_list_name + " list is not of equal length in Original and New Genomes.")
        print("\n\n" + feature_list_name + " TOTAL NUMBER:" + str(len(genome_orig[feature_list_name])))
        orig_dict = dict([(x['id'],x) for x in genome_orig[feature_list_name]])
        new_dict = dict([(x['id'],x) for x in genome_new[feature_list_name]])

        first_pass_matches = 0
        first_pass_non_match = 0
        second_pass_matches = 0

        for key in orig_dict:
            orig_feature = orig_dict[key]
            new_feature = new_dict[key]
            if "aliases" in orig_feature:
                orig_feature['aliases'] = sorted(orig_feature['aliases'])
                new_feature['aliases'] = sorted(new_feature['aliases'])
            if "db_xrefs" in orig_feature:
                orig_feature['db_xrefs'] = sorted(orig_feature['db_xrefs'])
                new_feature['db_xrefs'] = sorted(new_feature['db_xrefs'])
            if "functions" in orig_feature:
                orig_feature["functions"] = sorted(orig_feature["functions"])
                new_feature["functions"] = sorted(new_feature["functions"])
            if orig_feature == new_feature:
                first_pass_matches += 1
            else:
                first_pass_non_match += 1
                note_orig = orig_feature.pop("note",None)
                note_new = new_feature.pop("note",None)
                inference_orig = orig_feature.pop('inference_data',None)
                inference_new = new_feature.pop('inference_data',None)
                if "warnings" in orig_feature:
                    if 'The coordinates supplied for this feature are non-exact. DNA or protein translations are approximate.' in  orig_feature["warnings"]:
                        orig_feature["warnings"] = [warning for warning in orig_feature["warnings"] if warning != 
                            'The coordinates supplied for this feature are non-exact. DNA or protein translations are approximate.']
                if "warnings" in orig_feature and "warnings" not in new_feature:
                    del(orig_feature["warnings"])

                if orig_feature == new_feature:
                    second_pass_matches += 1
                else:
                    self.maxDiff = None
                    self.assertEqual(orig_feature,new_feature)        
        self.assertEqual(len(orig_dict),(first_pass_matches + second_pass_matches),
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
    
    def test_genome_level_attributes(self):
        genome_orig = self.__class__.genome_orig
        genome_new = self.__class__.genome_new 
        self.maxDiff = None       
        self.assertEqual(genome_orig["scientific_name"],genome_new["scientific_name"])
        self.assertEqual(genome_orig["domain"],genome_new["domain"])
        self.assertEqual(genome_orig["genome_tiers"],genome_new["genome_tiers"])
        self.assertEqual(genome_orig["genetic_code"],genome_new["genetic_code"])
        self.assertEqual(genome_orig["dna_size"],genome_new["dna_size"])
        self.assertEqual(genome_orig["num_contigs"],genome_new["num_contigs"])
        self.assertEqual(genome_orig["contig_lengths"],genome_new["contig_lengths"])
        self.assertEqual(genome_orig["contig_ids"],genome_new["contig_ids"])
        self.assertEqual(genome_orig["source"],genome_new["source"])
        self.assertEqual(genome_orig["md5"],genome_new["md5"])
        self.assertEqual(genome_orig["taxonomy"],genome_new["taxonomy"])
        self.assertEqual(genome_orig["gc_content"],genome_new["gc_content"])
        for publication in genome_orig["publications"]:
            self.assertTrue(publication in genome_new["publications"])
        self.assertEqual(genome_orig["ontologies_present"],genome_new["ontologies_present"])
