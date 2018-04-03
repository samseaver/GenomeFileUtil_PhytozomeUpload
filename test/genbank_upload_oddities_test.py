import unittest
import time
import os
import re

try:
    from ConfigParser import ConfigParser  # py2
except:
    from configparser import ConfigParser  # py3

from Workspace.WorkspaceClient import Workspace as workspaceService
from GenomeFileUtil.GenomeFileUtilImpl import GenomeFileUtil
from GenomeFileUtil.GenomeFileUtilServer import MethodContext
from DataFileUtil.DataFileUtilClient import DataFileUtil
from GenomeFileUtil.core.GenomeUtils import warnings
import json

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
              'generate_ids_if_needed': 1,
              'source': "RefSeq Latest"
            })[0]
        data_file_cli = DataFileUtil(os.environ['SDK_CALLBACK_URL'],
                                token=cls.ctx['token'],
                                service_ver='dev')
        cls.genome = data_file_cli.get_objects({'object_refs': [result['genome_ref']]})['data'][0]['data']
        json.dump(cls.genome, open(cls.cfg['scratch']+"/test_genome.json", 'w'))

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    #Note tested duplicate locus tags and got a value error which is good.
    #TO induce this error again add in the following lines into the GenBank file.
    #note if these lines are added to the source file, it will fail uploading.
#    gene            1..6
#                     /gene="trnH_duplicate_locus_tag"
#                     /locus_tag="ArthCt088"
#                     /db_xref="GeneID:1466273"
#     CDS             1..6
#                     /gene="trnH_duplicate_locus_tag"
#                     /locus_tag="ArthCt088"
#                     /db_xref="GeneID:1466273"

    def test_refseq_latest_source_and_tiers(self):
        genome = self.__class__.genome
        has_genome_tiers = False
        has_external_db = False
        if "genome_tier" in genome:
            has_genome_tiers = True
            for tier in genome["genome_tier"]:           
                if tier == "ExternalDB" :
                    has_external_db = True
        self.assertTrue(genome.get("source") == "RefSeq", "Source is not RefSeq : " + str(genome.get("source")))
        self.assertTrue(has_genome_tiers, "Does not have Genome Tiers")
        self.assertTrue(len(genome["genome_tier"]) == 1, "Should only have 1 tier in it.")
        self.assertTrue(has_external_db, "Does not have ExternalDB Genome Tier")    

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

    def test_both_strand_trans_splicing(self):
        genome = self.__class__.genome
        found_gene = False
        found_cds = False
        gene_flag_found = False
        cds_flag_found = False
        for feature in genome["features"]:
            if feature['id'] == "ArthCp047":
                print "Found ArthCp047"
                found_gene = True
                if "flags" in feature:
                    for flag in feature["flags"]:
                        if flag == "trans_splicing":
                            gene_flag_found = True
                gene_sequence = "ATGCCAACCATTAAACAACTTATTAGAAATACAAGACAGCCAATCCGAAACGTCACGAAATCCCCAGCGCTTCGGGGATGCCCTCAGCGACGAGGAACATGTACTCGGGTGTATACTATCACCCCCAAAAAACCAAACTCTGCTTTACGTAAAGTTGCCAGAGTACGATTAACCTCGGGATTTGAAATCACTGCTTATATACCTGGTATTGGCCATAATTTACAAGAACATTCTGTAGTCTTAGTAAGAGGGGGAAGGGTTAAGGATTTACCCGGTGTGAGATATCACATTGTTCGAGGAACCCTAGATGCTGTCGGAGTAAAGGATCGTCAACAAGGGCGTTCTAGTGCGTTGTAGATTCTTATCCAAGACTTGTATCATTTGATGATGCCATGTGAATCGCTAGAAACATGTGAAGTGTATGGCTAACCCAATAACGAAAGTTTCGTAAGGGGACTGGAGCAGGCTACCACGAGACAAAAGATCTTCTTTCAAAAGAGATTCGATTCGGAACTCTTATATGTCCAAGGTTCAATATTGAAATAATTTCAGAGGTTTTCCCTGACTTTGTCCGTGTCAACAAACAATTCGAAATGCCTCGACTTTTTTAGAACAGGTCCGGGTCAAATAGCAATGATTCGAAGCACTTATTTTTACACTATTTCGGAAACCCAAGGACTCAATCGTATGGATATGTAAAATACAGGATTTCCAATCCTAGCAGGAAAAGGAGGGAAACGGATACTCAATTTAAAAGTGAGTAAACAGAATTCCATACTCGATTTCAGAGATACATATATAATTCTGTGGAAAGCCGTATTCGATGAAAGTCGTATGTACGGTTTGGAGGGAGATCTTTCATATCTTTCGAGATCCACCCTACAATATGGGGTCAAAAAGCCAAAATAA"
                self.assertTrue(feature['dna_sequence'] == gene_sequence, "The DNA sequence for the gene ArthCp047 was not as expected. It contained the following sequence : " + str(feature['dna_sequence']))
                if 'cdss' in feature:
                    self.assertTrue("ArthCp047_CDS_1" in feature['cdss'], "The child CDS of ArthCp047 was not found")
                else:
                    self.assertTrue('cdss' in feature, "There was no child CDS for ArthCp047")
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
        for feature in genome["features"]:
            if feature['id'] == "InvalidOrder":
                # print "FEATURE::::" + str(feature)
                print "Found InvalidOrder"
                found_gene = True
                if "flags" in feature:
                    for flag in feature["flags"]:
                        if flag == "trans_splicing":
                            gene_transpliced_flag = True
                if "warnings" in feature:
                    if warnings['not_trans_spliced'] in feature["warnings"]:
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
                    if warnings['not_trans_spliced'] in feature["warnings"]:
                        has_cds_warning = True
        self.assertTrue(has_cds_warning, "The position coordinates for CDS 'InvalidOrder' are out of order and this is not listed as a transpliced CDS.  It should have a warning.")
        self.assertFalse(cds_transpliced_flag, "The trans_splicing flag for the cds InvalidOrder_CDS_1 was set, technically it appears it may be transpliced, but the file does not state it to be.")
        if "suspect" in genome:
            if genome["suspect"] == 1:
                genome_suspect = True
        if "warnings" in genome:
            if warnings['genome_not_trans_spliced'].format(4) in genome["warnings"]:
                genome_warning = True
        self.assertTrue(genome_suspect, "This genome has invalid position order features in it. It should be deemed suspect.")
        self.assertTrue(genome_warning, "This Genome has feature(s) with invalid coordinates, and should have a genome level warning to reflect that." + str(genome["warnings"]))
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
                    if warnings['not_trans_spliced'] in feature["warnings"]:
                        has_gene_warning = True
        self.assertFalse(has_gene_warning, "The position coordinates are out of order, but they start and stop at the start and end of a circular contig, therefore it is valid.")
        self.assertFalse(gene_transpliced_flag, "The trans_splicing flag was set, technically it appears it may be transpliced, but the file does not state it to be.")
        for feature in genome["cdss"]:
            if feature['id'] == "RL4742A_CDS_1":
#                print "FEATURE::::" + str(feature)
                print "Found RL4742A_CDS_1"
                found_cds = True
                if "flags" in feature:
                    for flag in feature["flags"]:
                        if flag == "trans_splicing":
                            cds_transpliced_flag = True
                if "warnings" in feature:
                    if warnings['not_trans_spliced'] in feature["warnings"]:
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
                    if warnings['not_trans_spliced'] in feature["warnings"]:
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
                    if warnings['not_trans_spliced'] in feature["warnings"]:
                        has_cds_warning = True
        self.assertFalse(has_cds_warning, "The position coordinates CDS 'RL4742' are out of order, but they start and stop at the start and end of a circular contig, therefore it is valid.")
        self.assertFalse(cds_transpliced_flag, "The trans_splicing flag for the cds RL4742_CDS_1 was set, technically it appears it may be transpliced, but the file does not state it to be.")
        self.assertTrue(found_gene, "The gene RL4742 was not found.")
        self.assertTrue(found_cds, "The CDS RL4742_CDS_1 was not found.")

    def test_for_zero_spanning_two_exon_feature(self):
        genome = self.__class__.genome
        found_gene = False
        found_cds = False
        gene_transpliced_flag = False
        cds_transpliced_flag = False
        has_gene_warning = False
        has_cds_warning = False
        genome_warning = False
        for feature in genome["features"]:
            if feature['id'] == "Zero_Span_two_exon":
#                print "FEATURE::::" + str(feature)
                print "Found Zero_Span_two_exon"
                found_gene = True
                if "flags" in feature:
                    for flag in feature["flags"]:
                        if flag == "trans_splicing":
                            gene_transpliced_flag = True
                if "warnings" in feature:
                    if warnings['not_trans_spliced'] in feature["warnings"]:
                        has_gene_warning = True
        self.assertFalse(has_gene_warning, "The position coordinates gene 'Zero_Span_two_exon' are out of order, but they start and stop at the start and end of a circular contig, therefore it is valid.")
        self.assertFalse(gene_transpliced_flag, "The trans_splicing flag for the gene Zero_Span_two_exon was set, but it is not trans_spliced.")
        for feature in genome["cdss"]:
            if feature['id'] == "Zero_Span_two_exon_CDS_1":
#                print "FEATURE::::" + str(feature)
                print "Found Zero_Span_two_exon_CDS_1"
                found_cds = True
                if "flags" in feature:
                    for flag in feature["flags"]:
                        if flag == "trans_splicing":
                            cds_transpliced_flag = True
                if "warnings" in feature:
                    if warnings['not_trans_spliced'] in feature["warnings"]:
                        has_cds_warning = True
        self.assertFalse(has_cds_warning, "The position coordinates CDS 'Zero_Span_two_exon' are out of order, but they are the child of a gene that start and stop at the start and end of a circular contig, therefore it is valid.")
        self.assertFalse(cds_transpliced_flag, "The trans_splicing flag for the cds Zero_Span_two_exon_CDS_1 was set, technically it appears it may be transpliced, but the file does not state it to be.")
        self.assertTrue(found_gene, "The gene Zero_Span_two_exon was not found.")
        self.assertTrue(found_cds, "The CDS Zero_Span_two_exon_CDS_1 was not found.")

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
                self.assertTrue("PO:0000005" in feature["ontology_terms"]["PO"],"PO:0000005 should be in the feature's ontology terms map ")
                self.assertTrue("PO:0000005" in genome["ontologies_present"]["PO"],"PO:0000005 should be in the ontologies_present map")
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
                        if warning == warnings["inconsistent_CDS_length"].format(240,88):
                            has_CDS_warning = True
                        if warning == warnings["inconsistent_translation"]:
                            has_translation_off_warning = True
        self.assertTrue(found_cds, "The cds ArthCp004_CDS_1 was not found.")
        self.assertTrue(has_CDS_warning, "The AA translation length is significantly off in ArthCp004_CDS_1")
        self.assertTrue(has_translation_off_warning, "The AA translation sequence is significantly off in ArthCp004_CDS_1")
        if "warnings" in genome:
            for warning in genome["warnings"]:
                if warning == warnings["genome_inc_CDS_length"].format("ArthCp004_CDS_1",240,88):
                    genome_translation_length_warning = True
                if warning == warnings["genome_inc_translation"].format(2,25):
                    genome_translation_warning = True
        if "suspect" in genome:
            if genome["suspect"] == 1:
                genome_suspect = True
        self.assertTrue(genome_translation_length_warning,"This does not have the translation length warning it should.")
        self.assertTrue(genome_translation_warning,"This does not have the CDS translation proportion warning it should." + str(genome["warnings"]))
        self.assertTrue(genome_suspect,"This has significant translation issues both in length and AA sequence.")

    def test_translation_not_supplied(self):
        genome = self.__class__.genome
        found_cds = False
        has_cds_warning = False
        has_translation = False
        for feature in genome["cdss"]:
            if feature['id'] == "ArthCp015_CDS_1":
                found_cds = True
                if "warnings" in feature:
                    for warning in feature["warnings"]:
                        if warning == warnings["no_translation_supplied"]:
                            has_cds_warning = True
                if feature.get("protein_translation"):
                    has_translation = True
        self.assertTrue(found_cds,"Did not find ArthCp015_CDS_1.")
        self.assertTrue(has_cds_warning,"Missing warning was derived from DNA Sequence.")
        self.assertTrue(has_translation,"The translation was derived and populated.")            

    def test_translation_not_supplied_not_multiple_of_3(self):
        genome = self.__class__.genome
        found_cds = False
        has_cds_warning = False
        has_empty_translation = False
        for feature in genome["cdss"]:
            if feature['id'] == "ArthCp015_non3_CDS_1":
                found_cds = True
                if "warnings" in feature:
                    for warning in feature["warnings"]:
                        #Note the later portion of this error is thrown by biopython
                        #if they change their error message, this will need to be updated.
                        if warning == warnings["no_translation_supplied"] + "Sequence length 88 is not a multiple of three" :
                            has_cds_warning = True
                if "protein_translation" in feature:
                    if feature["protein_translation"] == "":
                        has_empty_translation = True
        self.assertTrue(found_cds,"Did not find ArthCp015_non3_CDS_1.")
        self.assertTrue(has_cds_warning,"Missing warning unable to do translation.")
        self.assertTrue(has_empty_translation,"No trnaslation supplied and not a multiple of 3, should be emoty.")

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
                #Has the folowing warning: This CDS did not have a supplied translation. The translation is derived directly from DNA sequence.
                #self.assertFalse("warnings" in feature, "Since there is a ribosomal slippage, there should be any translation sequence off warnings.")

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
                found_cds = True
        self.assertFalse(found_gene,"'ArthCp085' Should not have been added to the genes, since it has invalid coordinates over the end of the contig")
        self.assertFalse(found_cds,"'ArthCp085__CDS_1' Should not have been added to the CDS, since it has invalid coordinates over the end of the contig")
        if "warnings" in genome:
            for warning in genome["warnings"]:
                if warning == warnings["coordinates_off_end"].format('ArthCp085_gene'):
                    found_gene_warning = True
                if warning == warning == warnings["coordinates_off_end"].format('ArthCp085_CDS'):
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
                    if warnings["non_exact_coordinates"] in feature["warnings"]:
                        found_gene_warning = True
        for feature in genome["cdss"]:
            if feature['id'] == "ArthCp005_CDS_1":
                found_cds = True
                #print "CDS : " + str(feature)
                if "warnings" in feature:
                    if warnings["non_exact_coordinates"] in feature["warnings"]:
                        found_cds_warning = True
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
                    if warnings["non_exact_coordinates"] in feature["warnings"]:
                        found_gene_warning = True
        for feature in genome["cdss"]:
            if feature['id'] == "ArthCp006_CDS_1":
                found_cds = True
                #print "CDS : " + str(feature)
                if "warnings" in feature:
                    if warnings["non_exact_coordinates"] in feature["warnings"]:
                        found_cds_warning = True
        self.assertTrue(found_gene,"'ArthCp006' Was not found in the genes")
        self.assertTrue(found_cds,"'ArthCp006_CDS_1' Was not found in the cdss")
        self.assertTrue(found_gene_warning,"Did not have gene unknown upper bound warning")
        self.assertTrue(found_cds_warning,"Did not have cds unknown upper bound warning")

    def test_odd_coordinates_unknown_both_bounds(self):
        genome = self.__class__.genome
        found_cds = False
        found_gene = False
        found_cds_upper_warning = False
        found_gene_upper_warning = False
        for feature in genome["features"]:
            if feature['id'] == "ArthCp007":
                found_gene = True
                #print "Gene : " + str(feature)
                if "warnings" in feature:
                    if warnings["non_exact_coordinates"] in feature["warnings"]:
                        found_gene_upper_warning = True
        for feature in genome["cdss"]:
            if feature['id'] == "ArthCp007_CDS_1":
                found_cds = True
                #print "CDS : " + str(feature)
                if "warnings" in feature:
                    if warnings["non_exact_coordinates"] in feature["warnings"]:
                        found_cds_upper_warning = True
        self.assertTrue(found_gene,"'ArthCp007' Was not found in the genes")
        self.assertTrue(found_cds,"'ArthCp007_CDS_1' Was not found in the cdss")
        self.assertTrue(found_gene_upper_warning,"Did not have gene unknown upper bound warning")
        self.assertTrue(found_cds_upper_warning,"Did not have cds unknown upper bound warning")

'''
#BIOPYTHON CURRENTLY DOES NOT HANDLE THIS.
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


