import unittest
import time
import os
import shutil
import re

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
        found_noncoding = False
        for feature in genome["features"]:
            if feature['id'] == "ArthCp001A":
#                print "FEATURE::::" + str(feature)
                print "Found ArthCp001A"
                found_gene = True
                self.assertFalse('mrnas' in feature, "Their should be no child mRNAs for ArthCp001A, should have failed on coordinates.")
                self.assertFalse('cdss' in feature, "Their should be no child CDSs for ArthCp001A, should have failed on coordinates.")
##TODO
#ADD CHECKS FOR WARNINGS                
        for feature in genome["mrnas"]:
            if feature['id'] == "ArthCp001A_mRNA_1":
#                print "MRNA::::" + str(feature)                
                print "Found ArthCp001A_mRNA_1"
                found_mRNA = True
                self.assertFalse('parent_gene' in feature, "There should be no parent_gene for ArthCp001A_mRNA_1, should have failed on coordinates.")
##TODO
#ADD CHECKS FOR WARNINGS
        for feature in genome["cdss"]:
            if feature['id'] == "ArthCp001A_CDS_1":
                print "Found ArthCp001A_CDS_1"
#                print "CDS::::" + str(feature)
                found_CDS = True
                self.assertFalse('parent_gene' in feature, "There should be no parent_gene for ArthCp001A_CDS_1, should have failed on coordinates.")
##TODO
#ADD CHECKS FOR WARNINGS
        self.assertTrue(found_gene, "The gene ArthCp001A was not found.")
        self.assertTrue(found_mRNA, "The mRNA ArthCp001A_mRNA_1 was not found.")
        self.assertTrue(found_CDS, "The CDS ArthCp001A_CDS_1 was not found.")   

    def test_both_strand_trans_splicing(self):
        genome = self.__class__.genome
        found_gene = False
        found_cds = False
        for feature in genome["features"]:
            if feature['id'] == "ArthCp047":
                print "Found ArthCp047 :: " + str(feature)
                found_gene = True
                if "flags" in feature:
                    for flag in feature["flags"]:
                        if flag == "trans_splicing":
                            gene_flag_found = True
                gene_sequence = "ATGCCAACCATTAAACAACTTATTAGAAATACAAGACAGCCAATCCGAAACGTCACGAAATCCCCAGCGCTTCGGGGATGCCCTCAGCGACGAGGAACATGTACTCGGGTGTATACTATCACCCCCAAAAAACCAAACTCTGCTTTACGTAAAGTTGCCAGAGTACGATTAACCTCGGGATTTGAAATCACTGCTTATATACCTGGTATTGGCCATAATTTACAAGAACATTCTGTAGTCTTAGTAAGAGGGGGAAGGGTTAAGGATTTACCCGGTGTGAGATATCACATTGTTCGAGGAACCCTAGATGCTGTCGGAGTAAAGGATCGTCAACAAGGGCGTTCTAGTGCGTTGTAGATTCTTATCCAAGACTTGTATCATTTGATGATGCCATGTGAATCGCTAGAAACATGTGAAGTGTATGGCTAACCCAATAACGAAAGTTTCGTAAGGGGACTGGAGCAGGCTACCACGAGACAAAAGATCTTCTTTCAAAAGAGATTCGATTCGGAACTCTTATATGTCCAAGGTTCAATATTGAAATAATTTCAGAGGTTTTCCCTGACTTTGTCCGTGTCAACAAACAATTCGAAATGCCTCGACTTTTTTAGAACAGGTCCGGGTCAAATAGCAATGATTCGAAGCACTTATTTTTACACTATTTCGGAAACCCAAGGACTCAATCGTATGGATATGTAAAATACAGGATTTCCAATCCTAGCAGGAAAAGGAGGGAAACGGATACTCAATTTAAAAGTGAGTAAACAGAATTCCATACTCGATTTCAGAGATACATATATAATTCTGTGGAAAGCCGTATTCGATGAAAGTCGTATGTACGGTTTGGAGGGAGATCTTTCATATCTTTCGAGATCCACCCTACAATATGGGGTCAAAAAGCCAAAATAA"                
                self.assertTrue(feature['dna_sequence'] == gene_sequence, "The DNA sequence for the gene ArthCp047 was not as expected. It contained the following sequence : " + str(feature['dna_sequence']))         
                if 'cdss' in feature:
                    self.assertTrue("ArthCp001_CDS_1" in feature['cdss'], "The child CDS of ArthCp047_CDS_1 was not found")
                else:
                    self.assertTrue('cdss' in feature, "There was no child CDS for ArthCp047")    
                self.assertTrue('mrnas' in feature in feature['mrnas'], "Their should not be children mrnas for ArthCp047") 
        self.assertTrue(found_gene, "The gene ArthCp047 was not found.")      
        self.assertTrue(gene_flag_found, "The trans_splicing flag for the gene ArthCp047 was not found.")
        for feature in genome["cdss"]:
            if feature['id'] == "ArthCp047_CDS_1":
                found_cds = True
                print "Found ArthCp047_CDS_1"
                if "flags" in feature:
                    for flag in feature["flags"]:
                        if flag == "trans_splicing":
                            cds_flag_found = True
                cds_sequence = "ATGCCAACCATTAAACAACTTATTAGAAATACAAGACAGCCAATCCGAAACGTCACGAAATCCCCAGCGCTTCGGGGATGCCCTCAGCGACGAGGAACATGTACTCGGGTGTATACTATCACCCCCAAAAAACCAAACTCTGCTTTACGTAAAGTTGCCAGAGTACGATTAACCTCGGGATTTGAAATCACTGCTTATATACCTGGTATTGGCCATAATTTACAAGAACATTCTGTAGTCTTAGTAAGAGGGGGAAGGGTTAAGGATTTACCCGGTGTGAGATATCACATTGTTCGAGGAACCCTAGATGCTGTCGGAGTAAAGGATCGTCAACAAGGGCGTTCTAAATATGGGGTCAAAAAGCCAAAATAA"
                self.assertTrue(feature['dna_sequence'] == cds_sequence, "The DNA sequence for the cds ArthCp047_CDS_1 was not as expected. It contained the following sequence : " + str(feature['dna_sequence']))         
                cds_translation = "MPTIKQLIRNTRQPIRNVTKSPALRGCPQRRGTCTRVYTITPKKPNSALRKVARVRLTSGFEITAYIPGIGHNLQEHSVVLVRGGRVKDLPGVRYHIVRGTLDAVGVKDRQQGRSKYGVKKPK"
#                print "FEATURE::::" + str(feature)
                self.assertTrue(feature['protein_translation'] == cds_translation, "The AA sequence for the cds ArthCp001_CDS_1 was not as expected. It contained the following sequence : " + str(feature['protein_translation']))         
                if 'parent_gene' in feature:
                    self.assertTrue(feature['parent_gene'] == 'ArthCp047',"The parent gene for ArthCp047_CDS_1 was not as expected")
                else:
                    self.assertTrue('parent_gene' in feature, "The parent gene for ArthCp047_CDS_1 was not populated")
        self.assertTrue(found_cds, "The cds ArthCp047_CDS_1 was not found.")  
        self.assertTrue(cds_flag_found, "The trans_splicing flag for the CDS ArthCp047_CDS_1 was not found.")

    def test_for_trans_splicing_multicontig(self):
        genome = self.__class__.genome
        gene_flag_found = False
        cds_flag_found = False
        for feature in genome["features"]:
            if feature['id'] == "MultiContigTransSpliced":
                print "Found MultiContigTransSpliced"
                if "flags" in feature:
                    for flag in feature["flags"]:
                        if flag == "trans_splicing":
                            gene_flag_found = True
                gene_sequence = "ATGCCAACCATTAAACAACTTATTAGAAATACAAGACAGCCAATCCGAAACGTCACGAAATCCCCAGCGCTTCGGGGATGCCCTCAGCGACGAGGAACATGTACTCGGGTGTATACTATCACCCCCAAAAAACCAAACTCTGCTTTACGTAAAGTTGCCAGAGTACGATTAACCTCGGGATTTGAAATCACTGCTTATATACCTGGTATTGGCCATAATTTACAAGAACATTCTGTAGTCTTAGTAAGAGGGGGAAGGGTTAAGGATTTACCCGGTGTGAGATATCACATTGTTCGAGGAACCCTAGATGCTGTCGGAGTAAAGGATCGTCAACAAGGGCGTTCTAGTGCGTTGTAGATTCTTATCCAAGACTTGTATCATTTGATGATGCCATGTGAATCGCTAGAAACATGTGAAGTGTATGGCTAACCCAATAACGAAAGTTTCGTAAGGGGACTGGAGCAGGCTACCACGAGACAAAAGATCTTCTTTCAAAAGAGATTCGATTCGGAACTCTTATATGTCCAAGGTTCAATATTGAAATAATTTCAGAGGTTTTCCCTGACTTTGTCCGTGTCAACAAACAATTCGAAATGCCTCGACTTTTTTAGAACAGGTCCGGGTCAAATAGCAATGATTCGAAGCACTTATTTTTACACTATTTCGGAAACCCAAGGACTCAATCGTATGGATATGTAAAATACAGGATTTCCAATCCTAGCAGGAAAAGGAGGGAAACGGATACTCAATTTAAAAGTGAGTAAACAGAATTCCATACTCGATTTCAGAGATACATATATAATTCTGTGGAAAGCCGTATTCGATGAAAGTCGTATGTACGGTTTGGAGGGAGATCTTTCATATCTTTCGAGATCCACCCTACAATATGGGGTCAAAAAGCCAAAAATGTAG"
                self.assertTrue(feature['dna_sequence'] == gene_sequence, "The DNA sequence for the gene MultiContigTransSpliced was not as expected. It contained the following sequence : " + str(feature['dna_sequence']))         
                if 'cdss' in feature:
                    self.assertTrue("MultiContigTransSpliced_CDS_1" in feature['cdss'], "The child CDS of MultiContigTransSpliced was not found")
                else:
                    self.assertTrue('cdss' in feature, "There was no child CDS for MultiContigTransSpliced")           
        self.assertTrue(gene_flag_found, "The trans_splicing flag for the gene MultiContigTransSpliced was not found.")
        for feature in genome["cdss"]:
            if feature['id'] == "MultiContigTransSpliced_CDS_1":
                print "Found MultiContigTransSpliced_CDS_1"
                if "flags" in feature:
                    for flag in feature["flags"]:
                        if flag == "trans_splicing":
                            cds_flag_found = True
                cds_sequence = "ATGCCAACCATTAAACAACTTATTAGAAATACAAGACAGCCAATCCGAAACGTCACGAAATCCCCAGCGCTTCGGGGATGCCCTCAGCGACGAGGAACATGTACTCGGGTGTATACTATCACCCCCAAAAAACCAAACTCTGCTTTACGTAAAGTTGCCAGAGTACGATTAACCTCGGGATTTGAAATCACTGCTTATATACCTGGTATTGGCCATAATTTACAAGAACATTCTGTAGTCTTAGTAAGAGGGGGAAGGGTTAAGGATTTACCCGGTGTGAGATATCACATTGTTCGAGGAACCCTAGATGCTGTCGGAGTAAAGGATCGTCAACAAGGGCGTTCTAAATATGGGGTCAAAAAGCCAAAAATGTAG"
                self.assertTrue(feature['dna_sequence'] == cds_sequence, "The DNA sequence for the cds MultiContigTransSpliced_CDS_1 was not as expected. It contained the following sequence : " + str(feature['dna_sequence']))         
                cds_translation = "MPTIKQLIRNTRQPIRNVTKSPALRGCPQRRGTCTRVYTITPKKPNSALRKVARVRLTSGFEITAYIPGIGHNLQEHSVVLVRGGRVKDLPGVRYHIVRGTLDAVGVKDRQQGRSKYGVKKPKM"
#                print "FEATURE::::" + str(feature)
                self.assertTrue(feature['protein_translation'] == cds_translation, "The AA sequence for the cds MultiContigTransSpliced_CDS_1 was not as expected. It contained the following sequence : " + str(feature['protein_translation']))         
                if 'parent_gene' in feature:
                    self.assertTrue(feature['parent_gene'] == 'MultiContigTransSpliced',"The parent gene for MultiContigTransSpliced_CDS_1 was not as expected")
                else:
                    self.assertTrue('parent_gene' in feature, "The parent gene for MultiContigTransSpliced_CDS_1 was not populated")
        self.assertTrue(cds_flag_found, "The trans_splicing flag for the CDS MultiContigTransSpliced_CDS_1 was not found.")


    def test_for_invalid_order(self):
        genome = self.__class__.genome
        found_gene = False
        found_cds = False
        gene_transpliced_flag = False
        cds_transpliced_flag = False
        has_gene_warning = False
        has_cds_warning = False
        genome_suspect = False
        genome_warning = False
        p = re.compile("SUSPECT: This Genome has \d+ features with coordinates that are out of order and are not trans_splicing")
        for feature in genome["features"]:
            if feature['id'] == "InvalidOrder":
#                print "FEATURE::::" + str(feature)
                print "Found InvalidOrder"
                found_gene = True
                if "flags" in feature:
                    for flag in feature["flags"]:
                        if flag == "trans_splicing":
                            gene_transpliced_flag = True
                if "warnings" in feature:
                    for warning in feature["warnings"]:
                        if warning == "The feature coordinates order are suspect and the feature is not listed as being trans_splicing":
                            has_gene_warning = True
        self.assertTrue(has_gene_warning, "The position coordinates for gene 'InvalidOrder' are out of order and this is not listed as a transpliced gene.  It should have a warning.")                
        self.assertFalse(gene_transpliced_flag, "The trans_splicing flag for the gene 'InvalidOrder' was set, technically it appears it may be transpliced, but the file does not state it to be.")
        for feature in genome["cdss"]:
            if feature['id'] == "InvalidOrder_CDS_1":
                print "Found InvalidOrder_CDS_1"
                found_cds = True
                if "flags" in feature:
                    for flag in feature["flags"]:
                        if flag == "trans_splicing":
                            cds_transpliced_flag = True
                if "warnings" in feature:
                    for warning in feature["warnings"]:
                        if warning == "The feature coordinates order are suspect and the feature is not listed as being trans_splicing":
                            has_cds_warning = True
        self.assertTrue(has_cds_warning, "The position coordinates for CDS 'InvalidOrder' are out of order and this is not listed as a transpliced CDS.  It should have a warning.")     
        self.assertFalse(cds_transpliced_flag, "The trans_splicing flag for the cds InvalidOrder_CDS_1 was set, technically it appears it may be transpliced, but the file does not state it to be.")
        if "suspect" in genome:
            if genome["suspect"] == 1:
                genome_suspect = True
        if "warnings" in genome:
            for warning in genome["warnings"]:
                m = p.match(warning)
                if m:
                    genome_warning = True            
        self.assertTrue(genome_suspect, "This genome has invalid position order features in it. It should be deemed suspect.")
        self.assertTrue(genome_warning, "This Genome has feature(s) with invalid coordinates, and should have a genome level warning to reflect that.")
        self.assertTrue(found_gene, "The gene InvalidOrder was not found.")
        self.assertTrue(found_cds, "The CDS InvalidOrder_CDS_1 was not found.")   

    def test_for_zero_spanning_pos_strand_feature(self):
        genome = self.__class__.genome
        found_gene = False
        found_cds = False
        gene_transpliced_flag = False
        cds_transpliced_flag = False
        has_gene_warning = False
        has_cds_warning = False
        genome_warning = False
        for feature in genome["features"]:
            if feature['id'] == "RL4742A":
#                print "FEATURE::::" + str(feature)
                print "Found RL4742A"
                found_gene = True
                if "flags" in feature:
                    for flag in feature["flags"]:
                        if flag == "trans_splicing":
                            gene_transpliced_flag = True
                if "warnings" in feature:
                    for warning in feature["warnings"]:
                        if warning == "The feature coordinates order are suspect and the feature is not listed as being trans_splicing":
                            has_gene_warning = True
        self.assertFalse(has_gene_warning, "The position coordinates are out of order, but they start and stop at the start and end of a circular contig, therefore it is valid.")                 
        self.assertFalse(gene_transpliced_flag, "The trans_splicing flag was set, technically it appears it may be transpliced, but the file does not state it to be.")
        for feature in genome["cdss"]:
            if feature['id'] == "InvalidOrder_CDS_1":
#                print "FEATURE::::" + str(feature)
                print "Found RL4742A_CDS_1"
                found_cds = True
                if "flags" in feature:
                    for flag in feature["flags"]:
                        if flag == "trans_splicing":
                            cds_transpliced_flag = True
                if "warnings" in feature:
                    for warning in feature["warnings"]:
                        if warning == "The feature coordinates order are suspect and the feature is not listed as being trans_splicing":
                            has_cds_warning = True
        self.assertFalse(has_cds_warning, "The position coordinates for CDS 'RL4742A' are out of order, but they start and stop at the start and end of a circular contig, therefore it is valid.")      
        self.assertFalse(cds_transpliced_flag, "The trans_splicing flag for the cds RL4742A_CDS_1 was set, technically it appears it may be transpliced, but the file does not state it to be.")
        self.assertTrue(found_gene, "The gene InvalidOrder was not found.")
        self.assertTrue(found_cds, "The CDS InvalidOrder_CDS_1 was not found.")   

    def test_for_zero_spanning_neg_strand_feature(self):
        genome = self.__class__.genome
        found_gene = False
        found_cds = False
        gene_transpliced_flag = False
        cds_transpliced_flag = False
        has_gene_warning = False
        has_cds_warning = False
        genome_warning = False
        for feature in genome["features"]:
            if feature['id'] == "RL4742":
#                print "FEATURE::::" + str(feature)
                print "Found RL4742"
                found_gene = True
                if "flags" in feature:
                    for flag in feature["flags"]:
                        if flag == "trans_splicing":
                            gene_transpliced_flag = True
                if "warnings" in feature:
                    for warning in feature["warnings"]:
                        if warning == "The feature coordinates order are suspect and the feature is not listed as being trans_splicing":
                            has_gene_warning = True
        self.assertFalse(has_gene_warning, "The position coordinates gene 'RL4742' are out of order, but they start and stop at the start and end of a circular contig, therefore it is valid.")                 
        self.assertFalse(gene_transpliced_flag, "The trans_splicing flag for the gene RL4742 was set, technically it appears it may be transpliced, but the file does not state it to be.")
        for feature in genome["cdss"]:
            if feature['id'] == "RL4742_CDS_1":
#                print "FEATURE::::" + str(feature)
                print "Found RL4742_CDS_1"
                found_cds = True
                if "flags" in feature:
                    for flag in feature["flags"]:
                        if flag == "trans_splicing":
                            cds_transpliced_flag = True
                if "warnings" in feature:
                    for warning in feature["warnings"]:
                        if warning == "The feature coordinates order are suspect and the feature is not listed as being trans_splicing":
                            has_cds_warning = True
        self.assertFalse(has_cds_warning, "The position coordinates CDS 'RL4742' are out of order, but they start and stop at the start and end of a circular contig, therefore it is valid.")      
        self.assertFalse(cds_transpliced_flag, "The trans_splicing flag for the cds RL4742_CDS_1 was set, technically it appears it may be transpliced, but the file does not state it to be.")
        self.assertTrue(found_gene, "The gene InvalidOrder was not found.")
        self.assertTrue(found_cds, "The CDS InvalidOrder_CDS_1 was not found.")   

    def test_ensembl_ontology_terms(self):
        genome = self.__class__.genome
        #print "ONTOLOGIES PRESENT: " + str(genome["ontologies_present"])
        found_cds = False
        for feature in genome["cdss"]:
            if feature['id'] == "ArthCp001_CDS_1":
                found_cds = True
                self.assertTrue("ontology_terms" in feature,"There are 2 Ensembl style ontology terms that should be accounted for.")
                #print "ONTOLOGY TERMS: " + str(feature["ontology_terms"])
                self.assertTrue("GO" in feature["ontology_terms"],"There is 1 Ensembl style ontology GO term that should be accounted for.")
                self.assertTrue("GO:0009523" in feature["ontology_terms"]["GO"],"GO:0009523 should be in the feature's ontology terms map ")
                self.assertTrue("GO:0009523" in genome["ontologies_present"]["GO"],"GO:0009523 should be in the ontologies_present map")
                self.assertTrue("PO" in feature["ontology_terms"],"There is 1 Ensembl style ontology PO term that should be accounted for.")
                self.assertTrue("PO:0000005" in feature["ontology_terms"]["GO"],"PO:0000005 should be in the feature's ontology terms map ")
                self.assertTrue("PO:0000005" in genome["ontologies_present"]["GO"],"PO:0000005 should be in the ontologies_present map")
        self.assertTrue(found_cds, "The cds ArthCp001_CDS_1 was not found.")  
        self.assertTrue("ontology_events" in genome,"There should be an ontology event for the upload.")       

    def test_drastic_translation_length_difference(self):
        genome = self.__class__.genome
        found_cds = False
        has_CDS_warning = False
        has_translation_off_warning = False
        genome_suspect = False
        genome_translation_length_warning = False  
        genome_translation_warning = False        
        for feature in genome["cdss"]:
            if feature['id'] == "ArthCp004_CDS_1":
                found_cds = True
                if "warnings" in feature:
                    for warning in feature["warnings"]:
                        if warning == "This CDS has a length of " + str(feature('dna_sequence_length')) + ", the translation included is too large at a length of 44 amino acids.":
                            has_CDS_warning = True
                        if warning == "The CDS DNA sequence is does not translate to the supplied protein translation.":
                            has_translation_off_warning = True
            self.assertTrue(has_CDS_warning, "The AA translation length is significantly off in ArthCp004_CDS_1")
            self.assertTrue(has_CDS_warning, "The AA translation sequence is significantly off in ArthCp004_CDS_1") 
        self.assertTrue(found_cds, "The cds ArthCp004_CDS_1 was not found.")     
        if "warnings" in genome:
            p = re.compile("SUSPECT: This Genome has \d+ CDS features significantly off translation lengths.")
            for warning in genome["warnings"]:
                m = p.match(warning)
                if m:
                    genome_translation_length_warning = True
            p = re.compile("SUSPECT: This Genome has a high proportion (\d+ out of \d+) CDS features that do not translate the supplied translation.")
            for warning in genome["warnings"]:
                m = p.match(warning)
                if m:
                    genome_translation_warning = True
        if "suspect" in genome:
            if genome["suspect"] == 1:
                genome_suspect = True
        self.assertTrue(genome_translation_length_warning,"This does not have the translation length warning it should.")
        self.assertTrue(genome_translation_warning,"This does not have the CDS translation proportion warning it should.")
        self.assertTrue(genome_suspect,"This has significant translation issues both in length and AA sequence.")
    
    def test_translation_not_supplied(self):
        genome = self.__class__.genome
        found_cds = False
        has_cds_warning = False
        has_translation = False
        for feature in genome["cdss"]:
            if feature['id'] == "RL4742_CDS_1":
                found_cds = True
                if "warnings" in feature:
                    for warning in feature["warnings"]:
                        if warning == "This CDS did not have a supplied translation. The translation is derived directly from DNA sequence.":
                            has_cds_warning = True 
                if "protein_translation" in feature:
                    has_translation = True
        self.assertTrue(found_cds,"Did not find RL4742_CDS_1.")
        self.assertTrue(has_cds_warning,"Missing warning was derived from DNA Sequence.")
        self.assertTrue(has_translation,"The translation was derived and populated.")            
               
    def test_ensembl_coordinates(self):
        genome = self.__class__.genome
        refseq_gene = None
        ensembl_gene = None
        refseq_CDS = None
        ensembl_CDS = None
        for feature in genome["features"]:
            if feature['id'] == "ArthCp004":
                refseq_gene = feature
            if feature['id'] == "ArthCp004_ENSEMBL":
                ensembl_gene = feature
        self.assertTrue(refseq_gene and ensembl_gene, "One or both of the following genes were not found: ArthCp004, ArthCp004_ENSEMBL")
        self.assertTrue(refseq_gene["location"] == ensembl_gene["location"],"The Ensembl style coordinates did not result in the same gene location information")
        self.assertTrue(refseq_gene["dna_sequence"] == ensembl_gene["dna_sequence"],"The Ensembl style coordinates did not result in the same gene sequence information")
        for feature in genome["cdss"]:
            if feature['id'] == "ArthCp004_CDS_1":
                refseq_CDS = feature
            if feature['id'] == "ArthCp004_ENSEMBL_CDS_1":
                ensembl_CDS = feature
        self.assertTrue(refseq_CDS and ensembl_CDS,"One or both of the following CDSs were not found: ArthCp004_CDS_1, ArthCp004_ENSEMBL_CDS_1")
#        print "ENSEMBL CDS LOCATIONS : " + str(ensembl_CDS["location"])
#        print "REFSEQ CDS LOCATIONS : " + str(refseq_CDS["location"])
        self.assertTrue(refseq_CDS["location"] == ensembl_CDS["location"],"The Ensembl style coordinates did not result in the same CDS location information: ENSEMBL CDS LOCATIONS : " +
                                                                        str(ensembl_CDS["location"]) + " --- REFSEQ CDS LOCATIONS : " + str(refseq_CDS["location"])) 
        self.assertTrue(refseq_CDS["dna_sequence"] == ensembl_CDS["dna_sequence"], "The Ensembl style coordinates did not result in the same CDS sequence information")
        cds_translation = "MVKLRLKRCGRKQRAVYRILAIDVRYRREGRDLSKVGFYDPITNQTFLNLSAILDFLKKGAQPTRTAHDISKKAGIFTE"
        self.assertTrue(ensembl_CDS['protein_translation'] == cds_translation, "The AA sequence for the cds ArthCp004_ENSEMBL_CDS_1 was not as expected. It contained the following sequence : " + str(feature['protein_translation']))            

    def test_check_ribosomal_slippage(self):
        genome = self.__class__.genome
        found_cds = False
        ribosomal_slippage_flag = False
        has_warnings = False
        for feature in genome["cdss"]:
            if feature['id'] == "SPO1_87_CDS_1":
                found_CDS = True
                if "flags" in feature:
                    for flag in feature["flags"]:
                        if flag == "ribosomal_slippage":
                            ribosomal_slippage_flag = True
                    self.assertTrue(ribosomal_slippage_flag,"The CDS SPO1_87_CDS_1 was supposed to have a ribosomal slippage flag")
                self.assertFalse("warnings" in feature,"Since there is a ribosomal slippage, there should be any translation sequence off warnings.")

    def test_invalid_coordinates_off_of_contig(self):
        genome = self.__class__.genome
        found_cds = False
        found_gene = False
        found_gene_warning = False
        found_cds_warning = False
        for feature in genome["features"]:
            if feature['id'] == "ArthCp085":
                found_gene = True
        for feature in genome["cdss"]:
            if feature['id'] == "ArthCp085_CDS_1":
                found_CDS = True
        self.assertFalse(found_gene,"'ArthCp085' Should not have been added to the genes, since it has invalid coordinates over the end of the contig")
        self.assertFalse(found_cds,"'ArthCp085__CDS_1' Should not have been added to the CDS, since it has invalid coordinates over the end of the contig")
        if "warnings" in genome:
            for warning in genome["warnings"]:
                if warning == "gene with the locus_tag ArthCp085' has invalid coordinates off of the end of the contig. This gene was not included.":
                    found_gene_warning = True
                if warning == "CDS with the locus_tag ArthCp085' has invalid coordinates off of the end of the contig. This CDS was not included.":
                    found_cds_warning = True
            self.assertTrue(found_gene_warning,"SUSPECT: The warning for the invalid gene off the end of the contig was not found.")
            self.assertTrue(found_cds_warning,"SUSPECT: The warning for the invalid CDS off the end of the contig was not found.")


    def test_odd_coordinates_unknown_lower_bound(self):
        genome = self.__class__.genome
        found_cds = False
        found_gene = False
        found_cds_warning = False
        found_gene_warning = False
        for feature in genome["features"]:
            if feature['id'] == "ArthCp005":
                found_gene = True
                #print "Gene : " + str(feature)
                if "warnings" in feature:
                    for warning in feature["warnings"]:
                        if warning == "The coordinates supplied for this feature are non-exact. DNA or protein translations are approximate.":
                            found_gene_warning = True
        for feature in genome["cdss"]:
            if feature['id'] == "ArthCp005_CDS_1":
                found_cds = True
                #print "CDS : " + str(feature)
                if "warnings" in feature:
                    for warning in feature["warnings"]:
                        if warning == "The coordinates supplied for this feature are non-exact. DNA or protein translations are approximate.":
                            found_gene_warning = True
        self.assertTrue(found_gene,"'ArthCp005' Was not found in the genes")
        self.assertTrue(found_cds,"'ArthCp005_CDS_1' Was not found in the cdss")
        self.assertTrue(found_gene_warning,"Did not have gene unknown lower bound warning")
        self.assertTrue(found_cds_warning,"Did not have cds unknown lower bound warning")

    def test_odd_coordinates_unknown_upper_bound(self):
        genome = self.__class__.genome
        found_cds = False
        found_gene = False
        found_cds_warning = False
        found_gene_warning = False
        for feature in genome["features"]:
            if feature['id'] == "ArthCp006":
                found_gene = True
                #print "Gene : " + str(feature)
                if "warnings" in feature:
                    for warning in feature["warnings"]:
                        if warning == "This feature has coordinates that end somewhere after base 7693. Note the sequence is an approximation and ends at base 7693":
                            found_gene_warning = True
        for feature in genome["cdss"]:
            if feature['id'] == "ArthCp006_CDS_1":
                found_cds = True
                #print "CDS : " + str(feature)
                if "warnings" in feature:
                    for warning in feature["warnings"]:
                        if warning == "This feature has coordinates that end somewhere after base 7693. Note the sequence is an approximation and ends at base 7693":
                            found_gene_warning = True
        self.assertTrue(found_gene,"'ArthCp006' Was not found in the genes")
        self.assertTrue(found_cds,"'ArthCp006_CDS_1' Was not found in the cdss")
        self.assertTrue(found_gene_warning,"Did not have gene unknown upper bound warning")
        self.assertTrue(found_cds_warning,"Did not have cds unknown upper bound warning")

    def test_odd_coordinates_unknown_both_bounds(self):
        genome = self.__class__.genome
        found_cds = False
        found_gene = False
        found_cds_upper_warning = False
        found_cds_lower_warning = False
        found_gene_upper_warning = False
        found_gene_lower_warning = False
        for feature in genome["features"]:
            if feature['id'] == "ArthCp007":
                found_gene = True
                #print "Gene : " + str(feature)
                if "warnings" in feature:
                    for warning in feature["warnings"]:
                        if warning == "This feature has coordinates that end somewhere after base 11461. Note the sequence is an approximation and ends at base 11461":
                            found_gene_upper_warning = True
                if "warnings" in feature:
                    for warning in feature["warnings"]:
                        if warning == "This feature has coordinates that starts somewhere before base 9938. Note the sequence is an approximation and starts at base 9938":
                            found_gene_lower_warning = True                            
        for feature in genome["cdss"]:
            if feature['id'] == "ArthCp007_CDS_1":
                found_cds = True
                #print "CDS : " + str(feature)
                if "warnings" in feature:
                    for warning in feature["warnings"]:
                        if warning == "This feature has coordinates that end somewhere after base 11461. Note the sequence is an approximation and ends at base 11461":
                            found_cds_upper_warning = True
                if "warnings" in feature:
                    for warning in feature["warnings"]:
                        if warning == "This feature has coordinates that starts somewhere before base 9938. Note the sequence is an approximation and starts at base 9938":
                            found_cds_lower_warning = True 
        self.assertTrue(found_gene,"'ArthCp007' Was not found in the genes")
        self.assertTrue(found_cds,"'ArthCp007_CDS_1' Was not found in the cdss")
        self.assertTrue(found_gene_upper_warning,"Did not have gene unknown upper bound warning")
        self.assertTrue(found_cds_upper_warning,"Did not have cds unknown upper bound warning")
        self.assertTrue(found_gene_lower_warning,"Did not have gene unknown lower bound warning")
        self.assertTrue(found_cds_lower_warning,"Did not have cds unknown lower bound warning")

'''
    def test_reversed_position(self):
        genome = self.__class__.genome
        found_gene = False
        genome_suspect = False
        genome_warning = False
        p = re.compile("SUSPECT: This Genome has invalid \d+ features with coordinates where the first position value is greater than the second individual value.  These features will be ignored.")
        for feature in genome["features"]:
            if feature['id'] == "REVERSED":
#                print "FEATURE::::" + str(feature)
                print "Found REVERSED"
                found_gene = True
        if "suspect" in genome:
            if genome["suspect"] == 1:
                genome_suspect = True
        self.assertTrue(genome_suspect,"This genome has reversed position features in it. It should be deemed suspect.")
        if "warnings" in genome:
            for warning in genome["warnings"]:
##TODO check for correct warning message (NOT SURE ON NUMBER YET IN THIS FILE)
                m = p.match(warning)
                if m:
                    genome_warning = True                 
        self.assertTrue(genome_warning, "This Genome has feature(s) with reversed coordinates, and should have a genome level warning to reflect that.")          
        self.assertTrue(found_gene, "The gene REVERSED was found, it should have been excluded.")
'''


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

