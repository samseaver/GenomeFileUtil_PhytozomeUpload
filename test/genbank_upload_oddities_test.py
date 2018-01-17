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
from GenomeFileUtil.core.GenomeUtils import is_parent
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
        gbk_path = "data/Arabidopsis_gbff/Arab_Chloro_Modified.gbff"
        ws_obj_name = 'ArabidopsisChloro'
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
        data_file_cli = DataFileUtil(os.environ['SDK_CALLBACK_URL'], 
                                token=cls.ctx['token'],
                                service_ver='dev')
        cls.genome = data_file_cli.get_objects({'object_refs': [result['genome_ref']]})['data'][0]['data']


#        print("GENE 1: ")
#        pprint(cls.genome['features'][0])
#        print("HERE ARE THE CDSs:")
#        pprint(cls.genome['cdss'])
#        pprint(result)

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def test_incorrect(self):
        self.assertTrue( 1 == 0, "1 ne 0")

    def test_for_alias_colon(self):
        genome = self.__class__.genome
        colon_included = False
        for feature in genome["cdss"]:
            if feature['id'] == "ArthCp001_CDS_1":
                print "Found ArthCp001_CDS_1"
                for alias_tuple in feature["aliases"]:
                    if alias_tuple[0] == "gene_synonym":
                        if alias_tuple[1] == "TEST:COLON":
                            colon_included = True 
        self.assertTrue(colon_included, "The synonym TEST:COLON was not found.")

    def test_for_trans_splicing(self):
        genome = self.__class__.genome
        gene_flag_found = False
        cds_flag_found = False
        found_mRNA = False
        for feature in genome["features"]:
            if feature['id'] == "ArthCp001":
                print "Found ArthCp001"
                if "flags" in feature:
                    for flag in feature["flags"]:
                        if flag == "trans_splicing":
                            gene_flag_found = True
                gene_sequence = "ATGCCAACCATTAAACAACTTATTAGAAATACAAGACAGCCAATCCGAAACGTCACGAAATCCCCAGCGCTTCGGGGATGCCCTCAGCGACGAGGAACATGTACTCGGGTGTATACTATCACCCCCAAAAAACCAAACTCTGCTTTACGTAAAGTTGCCAGAGTACGATTAACCTCGGGATTTGAAATCACTGCTTATATACCTGGTATTGGCCATAATTTACAAGAACATTCTGTAGTCTTAGTAAGAGGGGGAAGGGTTAAGGATTTACCCGGTGTGAGATATCACATTGTTCGAGGAACCCTAGATGCTGTCGGAGTAAAGGATCGTCAACAAGGGCGTTCTAGTGCGTTGTAGATTCTTATCCAAGACTTGTATCATTTGATGATGCCATGTGAATCGCTAGAAACATGTGAAGTGTATGGCTAACCCAATAACGAAAGTTTCGTAAGGGGACTGGAGCAGGCTACCACGAGACAAAAGATCTTCTTTCAAAAGAGATTCGATTCGGAACTCTTATATGTCCAAGGTTCAATATTGAAATAATTTCAGAGGTTTTCCCTGACTTTGTCCGTGTCAACAAACAATTCGAAATGCCTCGACTTTTTTAGAACAGGTCCGGGTCAAATAGCAATGATTCGAAGCACTTATTTTTACACTATTTCGGAAACCCAAGGACTCAATCGTATGGATATGTAAAATACAGGATTTCCAATCCTAGCAGGAAAAGGAGGGAAACGGATACTCAATTTAAAAGTGAGTAAACAGAATTCCATACTCGATTTCAGAGATACATATATAATTCTGTGGAAAGCCGTATTCGATGAAAGTCGTATGTACGGTTTGGAGGGAGATCTTTCATATCTTTCGAGATCCACCCTACAATATGGGGTCAAAAAGCCAAAATAA"
                self.assertTrue(feature['dna_sequence'] == gene_sequence, "The DNA sequence for the gene ArthCp001 was not as expected. It contained the following sequence : " + str(feature['dna_sequence']))         
                if 'cdss' in feature:
                    self.assertTrue("ArthCp001_CDS_1" in feature['cdss'], "The child CDS of ArthCp001_CDS_1 was not found")
                else:
                    self.assertTrue('cdss' in feature, "There was no child CDS for ArthCp001")    
                if 'mrnas' in feature:
                    self.assertTrue("ArthCp001_mRNA_1" in feature['mrnas'], "The child mRNA of ArthCp001_mRNA_1 was not found")
                else:
                    self.assertTrue('mrnas' in feature, "There was no child mRNA for ArthCp001")             
        self.assertTrue(gene_flag_found, "The trans_splicing flag for the gene ArthCp001 was not found.")
        for feature in genome["mrnas"]:
            if feature['id'] == "ArthCp001_mRNA_1":
                print "Found ArthCp001_mRNA_1"
                found_mRNA = True
                if 'parent_gene' in feature:
                    self.assertTrue(feature['parent_gene'] == 'ArthCp001',"The parent gene for ArthCp001_CDS_1 was not as expected")
                else:
                    self.assertTrue('parent_gene' in feature, "The parent gene for ArthCp001_CDS_1 was not populated")
        self.assertTrue(found_mRNA, "The mRNA ArthCp001_mRNA_1 was not found.")
        for feature in genome["cdss"]:
            if feature['id'] == "ArthCp001_CDS_1":
                print "Found ArthCp001_CDS_1"
                if "flags" in feature:
                    for flag in feature["flags"]:
                        if flag == "trans_splicing":
                            cds_flag_found = True
                cds_sequence = "ATGCCAACCATTAAACAACTTATTAGAAATACAAGACAGCCAATCCGAAACGTCACGAAATCCCCAGCGCTTCGGGGATGCCCTCAGCGACGAGGAACATGTACTCGGGTGTATACTATCACCCCCAAAAAACCAAACTCTGCTTTACGTAAAGTTGCCAGAGTACGATTAACCTCGGGATTTGAAATCACTGCTTATATACCTGGTATTGGCCATAATTTACAAGAACATTCTGTAGTCTTAGTAAGAGGGGGAAGGGTTAAGGATTTACCCGGTGTGAGATATCACATTGTTCGAGGAACCCTAGATGCTGTCGGAGTAAAGGATCGTCAACAAGGGCGTTCTAAATATGGGGTCAAAAAGCCAAAATAA"
                self.assertTrue(feature['dna_sequence'] == cds_sequence, "The DNA sequence for the cds ArthCp001_CDS_1 was not as expected. It contained the following sequence : " + str(feature['dna_sequence']))         
                cds_translation = "MPTIKQLIRNTRQPIRNVTKSPALRGCPQRRGTCTRVYTITPKKPNSALRKVARVRLTSGFEITAYIPGIGHNLQEHSVVLVRGGRVKDLPGVRYHIVRGTLDAVGVKDRQQGRSKYGVKKPK"
#                print "FEATURE::::" + str(feature)
                self.assertTrue(feature['protein_translation'] == cds_translation, "The AA sequence for the cds ArthCp001_CDS_1 was not as expected. It contained the following sequence : " + str(feature['protein_translation']))         
                if 'parent_gene' in feature:
                    self.assertTrue(feature['parent_gene'] == 'ArthCp001',"The parent gene for ArthCp001_CDS_1 was not as expected")
                else:
                    self.assertTrue('parent_gene' in feature, "The parent gene for ArthCp001_CDS_1 was not populated")
        self.assertTrue(cds_flag_found, "The trans_splicing flag for the CDS ArthCp001 was not found.")

    
    def test_for_trans_splicing_invalid_parentage(self):
        genome = self.__class__.genome
        found_gene = False
        found_CDS = False
        found_mRNA = False
        for feature in genome["features"]:
            if feature['id'] == "ArthCp001A":
                print "FEATURE::::" + str(feature)
                print "Found ArthCp001A"
                found_gene = True
                #self.assertFalse('mrnas' in feature, "Their should be no child mRNAs for ArthCp001A, should have failed on coordinates.")
                #self.assertFalse('cdss' in feature, "Their should be no child CDSs for ArthCp001A, should have failed on coordinates.")
        for feature in genome["mrnas"]:
            if feature['id'] == "ArthCp001A_mRNA_1":
                print "Found ArthCp001A_mRNA_1"
                found_mRNA = True
                print(feature['warnings'])
                self.assertFalse('parent_gene' in feature, "Their should be no parent_gene for ArthCp001A_mRNA_1, should have failed on coordinates.")
        for feature in genome["features"]:
            if feature['id'] == "ArthCp001A_CDS_1":
                print "Found ArthCp001A_CDS_1"
                found_CDS = True
                self.assertFalse('parent_gene' in feature, "Their should be no parent_gene for ArthCp001A_CDS_1, should have failed on coordinates.")
        self.assertTrue(found_gene, "The gene ArthCp001A was not found.")
        self.assertTrue(found_mRNA, "The mRNA ArthCp001A_mRNA_1 was not found.")
        self.assertTrue(found_CDS, "The CDS ArthCp001A_CDS_1 was not found.")   

    def test_is_parent(self):
        mrna = {'location': [
            ["Supercontig_5", 5204180, "-", 1901],
            ["Supercontig_5", 5202226, "-", 1412]
        ]}
        cds1 = {'location': [
            ["Supercontig_5", 5203923, "-", 1644],
            ["Supercontig_5", 5202226, "-", 1071]
        ]}
        cds2 = {'location': [
            ["Supercontig_5", 5203923, "-", 1644],
            ["Supercontig_5", 5202226, "+", 1071]
        ]}
        cds3 = {'location': [
            ["Supercontig_5", 5203923, "-", 1644],
            ["Supercontig_5", 5202228, "-", 1071]
        ]}
        self.assertTrue(is_parent(mrna, cds1))
        self.assertFalse(is_parent(mrna, cds2), "Mismatched strands")
        self.assertFalse(is_parent(mrna, cds3), "Mismatched start")


'''
#################################################
    def test_unknown_molecule(self):
        genome = self.__class__.genome
        if 'warnings' in genome:
            for genome_warning in genome['warnings']:
                print("WARNING")
                print(str(genome_warning)) 
                self.assertNotIn("Genome molecule_type Unknown is not expected for domain Bacteria.", genome_warning)

    def test_empty_publications(self):
        genome = self.__class__.genome
        if "publications" in genome:
            for publication in genome["publications"]:
                self.assertFalse((publication[0] == 0) and 
                                 (publication[1] == '') and
                                 (publication[2] == '') and
                                 (publication[3] == '') and
                                 (publication[4] == '') and
                                 (publication[5] == '') and
                                 (publication[6] == ''),
                                 "Stored an Empty Publication")

    def test_accurately_stored_publication(self):
        genome = self.__class__.genome
        found_publication = False
        if "publications" in genome:
            for publication in genome["publications"]:
#                print publication
                if ((publication[0] == 0) and 
                    (publication[1] == '') and
                    (publication[2] == 'Escherichia coli K-12 MG1655 yqiK-rfaE intergenic region, genomic sequence correction') and
                    (publication[3] == '') and
                    (publication[4] == '') and
                    (publication[5] == 'Perna,N.T.') and
                    (publication[6] == 'Unpublished')):
#                    print "FOUND THE PUBLICATION!!!!!!!!!!!!!!!!!"
                    found_publication = True
        self.assertTrue(found_publication,"Expected stored publication was not found.") 



    def test_neg_strand_off_by_one_issue(self):
        genome = self.__class__.genome
        found_synonyms = False
        for feature in genome["features"]:
            if feature['id'] == "b0618":
                self.assertTrue(feature["location"][0][1] == 651856, 
                                "The negative strand location start is off; " +
                                "It is " + str(feature["location"][0][1]) +
                                " when it should be 651856.")

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


    def test_getting_all_go_ontologies(self):
        genome = self.__class__.genome    
        all_ontologies_accounted_for = True
        check_all_go_ontologies = {
                                    "GO:0005737":0,
                                    "GO:0016563":0,
                                    "GO:0016564":0,
                                    "GO:0006350":0
        }
        for cds in genome["cdss"]:
            if cds['id'] == "b3357_CDS_1":
                print "Found b3357_CDS_1"
                for ontology in cds["ontology_terms"]["GO"]:
                    print "Ontology : " + str(ontology)
                    if ontology in check_all_go_ontologies:
                        check_all_go_ontologies[ontology] = 1
        for ontology in check_all_go_ontologies:
            if check_all_go_ontologies[ontology] == 0:
                all_ontologies_accounted_for = False            
        self.assertTrue(all_ontologies_accounted_for, "Not all expected ontologies were accounted for : " + str(check_all_go_ontologies))

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
                        for warning in feature["warnings"]:
                            if warning.strip() == '':
                                empty_warning_count += 1
                            else:
                                found_warning_count += 1
                    else:
                        features_without_warnings_count += 1
        print "EMPTY WARNING COUNT: " + str(empty_warning_count)
        print "FOUND WARNING COUNT: " + str(found_warning_count)
        print "FEATURES WITH WARNINGS COUNT: " + str(features_with_warnings_count)
        print "FEATURES WITHOUT WARNINGS COUNT: " + str(features_without_warnings_count)
        self.assertTrue(empty_warning_count == 0, str(empty_warning_count) + " features had empty warnings.")     
        self.assertTrue(found_warning_count > 0, "No features had warnings.")

#    def test_no_empty_mRNAs(self):
#        genome = self.__class__.genome
#        if "mrnas" in genome:
#            self.assertTrue(len(genome["mrnas"]) > 0, "The mRNA list is empty and is still present.")

#    def test_no_empty_genome_level_warnings(self):
#        genome = self.__class__.genome
#        if "warnings" in genome:
#            if len(genome["warnings"]) > 0:
#                for warning in genome["warnings"]:
#                    self.assertTrue(genome["warnings"][0] != '', "The Genome level warnings list is empty and is still present.")
#            self.assertTrue(len(genome["warnings"]) > 0, "The Genome level warnings is empty and is still present.")

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


    def test_flags_being_caught(self):
        genome = self.__class__.genome
        found_synonyms = False
        for feature in genome["features"]:
            if feature['id'] == "b4659":
#                print "Found b4659"
                self.assertTrue("flags" in feature, "This is a pseudo gene and should have the flags list.")
                found_pseudo = False
                for flag in feature["flags"]:
                    if flag == "pseudo":
                        found_pseudo = True
                self.assertTrue(found_pseudo, "This is a pseudo gene and should have a flag for it.")
'''

