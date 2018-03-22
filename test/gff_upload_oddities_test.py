import unittest
import time
import os
import shutil

try:
    from ConfigParser import ConfigParser  # py2
except:
    from configparser import ConfigParser  # py3

from Workspace.WorkspaceClient import Workspace as workspaceService
from GenomeFileUtil.GenomeFileUtilImpl import GenomeFileUtil
from GenomeFileUtil.GenomeFileUtilServer import MethodContext
from DataFileUtil.DataFileUtilClient import DataFileUtil
from GenomeFileUtil.core.GenomeUtils import warnings
from pprint import pprint


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
        gff_path = "data/fasta_gff/RefSeq/Red_Algae/RedAlgaeModified.gff"
        fna_path = "data/fasta_gff/RefSeq/Red_Algae/RedAlgaeModified.fna"
        ws_obj_name = 'red_algae_contigs'
        suffix = int(time.time() * 1000)
        cls.wsName = "test_GenomeFileUtil_" + str(suffix)
        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})
        result = cls.serviceImpl.fasta_gff_to_genome(
            cls.ctx,
            {
              'gff_file': {
                  'path': gff_path},
              'fasta_file': {
                  'path': fna_path},
              'workspace_name': cls.wsName,
              'genome_name': ws_obj_name,
              'generate_missing_genes' : 1,
              'generate_ids_if_needed': 1
            })[0]
#        print("HERE IS THE RESULT:")
        data_file_cli = DataFileUtil(os.environ['SDK_CALLBACK_URL'], 
                                token=cls.ctx['token'],
                                service_ver='dev')
        cls.genome = data_file_cli.get_objects({'object_refs': [result['genome_ref']]})['data'][0]['data']
#        print("GENE 1: ")
#        pprint(cls.genome['features'][0])
#        pprint(result)



    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

#    def test_incorrect(self):
#        self.assertTrue( 1 == 0, "1 ne 0")

    def test_feature_type_with_counts(self):
        genome = self.__class__.genome        
        list_names = ["features","cdss","mrnas","non_coding_features"]
        for list_name in list_names:
            if list_name in genome:
                self.assertTrue(len(genome[list_name]) > 0, list_name + " is empty.")

    def test_for_empty_functions(self):
        genome = self.__class__.genome
        empty_function_count = 0
        found_function_count = 0
        features_with_functions_count = 0
        features_without_functions_count = 0
        list_names = ["features","cdss","mrnas","non_coding_features"]
        for list_name in list_names:
            if list_name in genome:
                for feature in genome[list_name]:
                    if "functions" in feature:
                        features_with_functions_count += 1
                        for function in feature["functions"]:
                            if function.strip() == '':
                                empty_function_count += 1
                            else:
                                found_function_count += 1
                    else:
                        features_without_functions_count += 1
        print "EMPTY FUNCTION COUNT: " + str(empty_function_count)
        print "FOUND FUNCTION COUNT: " + str(found_function_count)
        print "FEATURES WITH FUNCTIONS COUNT: " + str(features_with_functions_count)
        print "FEATURES WITHOUT FUNCTIONS COUNT: " + str(features_without_functions_count)
        self.assertTrue(empty_function_count == 0, str(empty_function_count) + " features had empty functions.")     
        self.assertTrue(found_function_count > 0, "No features had functions.")    

    def test_for_not_multiple_of_3(self):
        #Tests for not a multiple of 3 warning. TEST CASE 
        genome = self.__class__.genome
        found_not_3multiple_warning = False
        found_making_upper_case = False
        for feature in genome['cdss']:
            if feature['id'] == 'rna0.CDS':
                #print "rna0.CDS WARNINGS: " + str(feature["warnings"])
                #print "rna0.CDS Sequence: " + str(feature["dna_sequence"])
                #print "rna0.CDS SubSequence: " + str(feature["dna_sequence"][:3])
                #print "rna0.CDS SubSequence Upper: " + warnings["non_standard_start_codon"].format(feature["dna_sequence"][:3].upper())
                #NEEDS TO BE MADE UPPER CASE.
                if feature["dna_sequence"][:5] == "ATGTT":
                    found_making_upper_case = True
                if warnings['not_multiple_of_3CDS'].format(len(feature["dna_sequence"])) in feature["warnings"]:
                    found_not_3multiple_warning = True   
                self.assertTrue(found_not_3multiple_warning, "Not multiple of 3 warning not found. Warnings: " + str(feature["warnings"]))  
                self.assertTrue(found_making_upper_case, "NOT MAKING THE SEQUENCE UPPER CASE : " + feature["dna_sequence"][:5])

    def test_for_out_of_order(self):
        #Tests for  out of order exons. Also test not a start codon
        genome = self.__class__.genome
        found_out_of_order_warning = False
        found_not_start_codon_warning = False
        for feature in genome['cdss']:
            if feature['id'] == 'rna1.CDS':
                if warnings["non_standard_start_codon"].format(feature["dna_sequence"][:3]) in feature["warnings"]:
                    found_not_start_codon_warning = True
                if warnings['out_of_order'] in feature["warnings"]:
                    found_out_of_order_warning = True
                self.assertTrue(found_out_of_order_warning, "Out of order warning not found. Warnings were: "
                            + str(feature["warnings"])) 
                self.assertTrue(found_not_start_codon_warning, "Not start codon warning not found. Warnings: " 
                            + str(feature["warnings"]))


#CHECK SEQUENCE

#JGI STRANDS

        
    def test_mrna_has_2_plus_locations(self):
        genome = self.__class__.genome
        for feature in genome['mrnas']:
            if feature['id'] == 'rna7':
                #print "RNA7 Feature: " + str(feature)
                #print "Feature locations : " + str(feature["location"])
                self.assertTrue(len(feature["location"]) == 2, "mRNA is not 2 locations as it should be.")               

    def test_mrna_has_2_minus_locations(self):
        genome = self.__class__.genome
        for feature in genome['mrnas']:
            if feature['id'] == 'rna3':
                #print "RNA3 Feature: " + str(feature)
                self.assertTrue(len(feature["location"]) == 2, "mRNA is not 2 locations as it should be.")    

    def test_CDS_has_2_plus_locations(self):
        genome = self.__class__.genome
        for feature in genome['cdss']:
            if feature['id'] == 'rna7.CDS':
                #print "RNA7.CDS Feature: " + str(feature)
                self.assertTrue(len(feature["location"]) == 2, "mRNA is not 2 locations as it should be.")                        

    def test_trans_splicing_both_strands(self):
        genome = self.__class__.genome
        #CHECK ALL THREE LEVELS. CHECK LOCATIONS. CHECK CATCHING TRANS_SPLICING. CHECK OUT OF ORDER ERROR.
        #NOTE CAN CATCH "exception=trans-splicing" in the last column
        gene_trans_splicing_flag = False
        mRNA_trans_splicing_flag = False
        CDS_trans_splicing_flag = False
        CDS_premature_stop_warning = False
        for feature in genome['features']:
            if feature['id'] == 'gene4ts':
                #print "GENE4ts Feature: " + str(feature))
                self.assertTrue(len(feature["location"]) == 2, "Gene is not 2 locations as it should be.")
                self.assertTrue(feature['location'][0] == [u'NC_010127.1', 69724, u'-', 114])
                self.assertTrue(feature['location'][1] == [u'NC_010127.1', 139856, u'+', 795]) 
                if "flags" in feature:
                    if "trans_splicing" in feature["flags"]:
                        gene_trans_splicing_flag = False           
        for feature in genome['mrnas']:
            if feature['id'] == 'rna4ts':
                #print "RNA4ts Feature: " + str(feature)
                self.assertTrue(len(feature["location"]) == 3, "mRNA is not 3 locations as it should be.")
                self.assertTrue(feature['location'][0] == ['NC_010127.1', 69724, '-', 114])
                self.assertTrue(feature['location'][1] == ['NC_010127.1', 139856, '+', 232])   
                self.assertTrue(feature['location'][2] == ['NC_010127.1', 140625, '+', 26]) 
                if "flags" in feature:
                    if "trans_splicing" in feature["flags"]:
                        mRNA_trans_splicing_flag = False   
        for feature in genome['cdss']:
            if feature['id'] == 'rna4ts.CDS':
                #print "RNA4ts.CDS Feature: " + str(feature)
                self.assertTrue(len(feature["location"]) == 3, "CDS is not 3 locations as it should be.")
                if "flags" in feature:
                    if "trans_splicing" in feature["flags"]:
                        CDS_trans_splicing_flag = False
                if "warnings" in feature:
                    if warnings["premature_stop_codon"] in feature["warnings"]:
                        CDS_premature_stop_warning = True
        self.assertTrue(CDS_premature_stop_warning, "Premature stop codon not found.")            
        self.assertTrue(gene_trans_splicing_flag, "gene trans_splicing flag not set.")
        self.assertTrue(mRNA_trans_splicing_flag, "gene trans_splicing flag not set.")   
        self.assertTrue(CDS_trans_splicing_flag, "gene trans_splicing flag not set.")

    def test_for_empty_feature_warnings(self):
        genome = self.__class__.genome
        empty_warning_count = 0
        found_warning_count = 0
        features_with_warnings_count = 0
        features_without_warnings_count = 0
        list_names = ["features","cdss","mrnas","non_coding_features"]
        for list_name in list_names:
            if list_name in genome:
                for feature in genome[list_name]:
                    if "warnings" in feature:
                        features_with_warnings_count += 1
                        print "FEATURE WITH A WARNING: " + str(feature)
                        for warning in feature["warnings"]:
                            if warning.strip() == '':
                                empty_warning_count += 1
                            else:
                                found_warning_count += 1
                    else:
                        features_without_warnings_count += 1
        print "EMPTY FEATURE WARNING COUNT: " + str(empty_warning_count)
        print "FOUND FEATURE WARNING COUNT: " + str(found_warning_count)
        print "FEATURES WITH WARNINGS COUNT: " + str(features_with_warnings_count)
        print "FEATURES WITHOUT WARNINGS COUNT: " + str(features_without_warnings_count)
        self.assertTrue(empty_warning_count == 0, str(empty_warning_count) + " features had empty warnings.")     
        self.assertTrue(found_warning_count > 0, "No features had warnings.")

    def test_non_coding_feature_ids(self):
        genome = self.__class__.genome
        underscore_start_count = 0
        overall_count = 0
        for feature in genome["non_coding_features"]:
            overall_count += 1
            if feature["id"].startswith("_"):
                underscore_start_count += 1 
        print "Starts with underscore count : " + str(underscore_start_count)
        print "Overall noncoding count : " + str(overall_count)
        self.assertTrue(underscore_start_count == 0, "Non coding features are starting with an underscore.")




