import unittest
from configparser import ConfigParser
from os import environ
import logging
import json

import mock

from GenomeFileUtil.GenomeFileUtilImpl import GenomeFileUtil, SDKConfig
from GenomeFileUtil.GenomeFileUtilServer import MethodContext
from GenomeFileUtil.core.GenomeInterface import GenomeInterface
from GenomeFileUtil.core import GenomeUtils
from installed_clients.WorkspaceClient import Workspace as workspaceService


class GenomeFileUtilTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
#        logging.basicConfig(level=logging.info)
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
        cls.serviceImpl = GenomeFileUtil(cls.cfg)
        gi_config = SDKConfig(cls.cfg)
        cls.genome_interface = GenomeInterface(gi_config)

    @classmethod
    def tearDownClass(cls):
        pass


    def test_non_unique_ids(self):
        # Non unique ids
        with open('data/non_unique_ids.json') as json_file:  
            non_unique_ids_genome = json.load(json_file)
        uniqueness_results = GenomeUtils.check_feature_ids_uniqueness(non_unique_ids_genome)
        self.assertTrue(len(uniqueness_results) > 0)
        self.assertTrue(uniqueness_results == {'RL742_CDS_1': 2, 'RL4742_mRNA_2': 2, 'RL4742': 3})

    def test_unique_ids(self):
        # Unique ids.
        with open('data/unique_ids.json') as json_file:  
            non_unique_ids_genome = json.load(json_file)
        uniqueness_results = GenomeUtils.check_feature_ids_uniqueness(non_unique_ids_genome)
        self.assertTrue(len(uniqueness_results) == 0)

    def test_single_feature_relationships_good(self):
        # Tests confirm_feature_relationships called separately for valid relationships.
        test_count = 0
        with open('data/good_relationships.json') as json_file:  
            good_relationships_genome = json.load(json_file)
        feature_id_sets_dict = dict()
        feature_lists = ["features","cdss","mrnas","non_coding_features"]
        for feature_list in feature_lists:
            if feature_list in good_relationships_genome:
                feature_id_sets_dict[feature_list] = GenomeUtils.make_id_set(good_relationships_genome[feature_list])
            else:
                feature_id_sets_dict[feature_list] = set()
        for feature in good_relationships_genome["features"]:
            if feature["id"] == "RL4742":
                gene_results = GenomeUtils.confirm_feature_relationships(feature,"features",feature_id_sets_dict)
                self.assertTrue(len(gene_results)==0)
                test_count += 1
        for feature in good_relationships_genome["cdss"]:
            if feature["id"] == "RL742_CDS_1":        
                cds_results = GenomeUtils.confirm_feature_relationships(feature,"cdss",feature_id_sets_dict)
                self.assertTrue(len(cds_results)==0)
                test_count += 1
        for feature in good_relationships_genome["mrnas"]:
            if feature["id"] == "RL742_mRNA_1":           
                mrna_results = GenomeUtils.confirm_feature_relationships(feature,"mrnas",feature_id_sets_dict)
                self.assertTrue(len(mrna_results)==0)
                test_count += 1
        for feature in good_relationships_genome["non_coding_features"]:
            if feature["id"] == "test_gene_1":         
                non_coding_parent_results = GenomeUtils.confirm_feature_relationships(feature,"non_coding_features",feature_id_sets_dict)
                self.assertTrue(len(non_coding_parent_results)==0) 
                test_count += 1
            elif feature["id"] == "test_gene_1_tRNA_2":         
                non_coding_nc_child_results = GenomeUtils.confirm_feature_relationships(feature,"non_coding_features",feature_id_sets_dict)
                self.assertTrue(len(non_coding_nc_child_results)==0) 
                test_count += 1
            elif feature["id"] == "RL4742_promoter":         
                non_coding_gene_child_results = GenomeUtils.confirm_feature_relationships(feature,"non_coding_features",feature_id_sets_dict)
                self.assertTrue(len(non_coding_gene_child_results)==0) 
                test_count += 1         
        self.assertTrue(test_count == 6,"Not all good inividual feature relationships tested. " + str(test_count))

    def test_single_feature_relationships_bad(self):
        # Tests confirm_feature_relationships called separately for missing relationship targets.
        with open('data/bad_relationships.json') as json_file:  
            bad_relationships_genome = json.load(json_file)
        feature_id_sets_dict = dict()
        feature_lists = ["features","cdss","mrnas","non_coding_features"]
        for feature_list in feature_lists:
            if feature_list in bad_relationships_genome:
                feature_id_sets_dict[feature_list] = GenomeUtils.make_id_set(bad_relationships_genome[feature_list])
            else:
                feature_id_sets_dict[feature_list] = set()
        test_count = 0
        for feature in bad_relationships_genome["features"]:
            if feature["id"] == "no_real_cds":
                gene_results1 = GenomeUtils.confirm_feature_relationships(feature,"features",feature_id_sets_dict)
                print("Bad Gene Results 1:" + str(gene_results1))
                self.assertTrue(len(gene_results1)>0)
                test_count += 1
            elif feature["id"] == "no_real_mrna":
                gene_results2 = GenomeUtils.confirm_feature_relationships(feature,"features",feature_id_sets_dict)
                print("Bad Gene Results 2:" + str(gene_results2))
                self.assertTrue(len(gene_results2)>0)
                test_count += 1
        for feature in bad_relationships_genome["cdss"]:
            if feature["id"] == "no_real_parent_gene_cds":        
                cds_results1 = GenomeUtils.confirm_feature_relationships(feature,"cdss",feature_id_sets_dict)
                print("Bad CDS Results 1:" + str(cds_results1))
                self.assertTrue(len(cds_results1)>0)
                test_count += 1
            elif feature["id"] == "no_real_parent_mrna_cds":  
                cds_results2 = GenomeUtils.confirm_feature_relationships(feature,"cdss",feature_id_sets_dict)
                print("Bad CDS Results 2:" + str(cds_results2))
                self.assertTrue(len(cds_results2)>0)
                test_count += 1
        for feature in bad_relationships_genome["mrnas"]:
            if feature["id"] == "no_real_parent_gene_mrna":       
                mrna_results1 = GenomeUtils.confirm_feature_relationships(feature,"mrnas",feature_id_sets_dict)
                print("Bad mRNA Results 1:" + str(mrna_results1))
                self.assertTrue(len(mrna_results1)>0)
                test_count += 1
            elif feature["id"] == "no_child_cds_mrna":  
                mrna_results2 = GenomeUtils.confirm_feature_relationships(feature,"mrnas",feature_id_sets_dict)
                print("Bad mRNA Results 2:" + str(mrna_results2))
                self.assertTrue(len(mrna_results2)>0)
                test_count += 1
        for feature in bad_relationships_genome["non_coding_features"]:
            if feature["id"] == "no_real_parent_gene_nc":   
                NC_results1 = GenomeUtils.confirm_feature_relationships(feature,"non_coding_features",feature_id_sets_dict)
                print("Bad NC Results 1:" + str(NC_results1))
                self.assertTrue(len(NC_results1)>0)
                test_count += 1
            elif feature["id"] == "no_real_nc_child":   
                NC_results2 = GenomeUtils.confirm_feature_relationships(feature,"non_coding_features",feature_id_sets_dict)
                print("Bad NC Results 2:" + str(NC_results2))
                self.assertTrue(len(NC_results2)>0)
                test_count += 1
        self.assertTrue(test_count == 8,"Not all bad inividual feature relationships tested." + str(test_count))

    def test_genome_relationships_good(self):
        # Tests confirm_feature_relationships for Genome with only valid relationships.
        with open('data/good_relationships.json') as json_file:  
            good_relationships_genome = json.load(json_file)
        genome_results = GenomeUtils.confirm_genomes_feature_relationships(good_relationships_genome)
        #print("Good Genome_results: " + str(genome_results))
        self.assertTrue(len(genome_results) == 0)

    def test_genome_relationships_bad(self):
        # Tests confirm_feature_relationships for Genome with only valid relationships.
        with open('data/bad_relationships.json') as json_file:  
            bad_relationships_genome = json.load(json_file)
        genome_results = GenomeUtils.confirm_genomes_feature_relationships(bad_relationships_genome)
        print("Bad Genome_results: " + str(genome_results))
        self.assertTrue(len(genome_results) > 0)

