import json
import os
import time
import unittest
from configparser import ConfigParser

from installed_clients.DataFileUtilClient import DataFileUtil
from GenomeFileUtil.GenomeFileUtilImpl import GenomeFileUtil
from GenomeFileUtil.GenomeFileUtilServer import MethodContext
from GenomeFileUtil.core.GenomeUtils import warnings
from installed_clients.WorkspaceClient import Workspace as workspaceService


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
        cls.genome = data_file_cli.get_objects(
            {'object_refs': [result['genome_ref']]}
        )['data'][0]['data']
        json.dump(cls.genome, open(cls.cfg['scratch']+"/test_genome.json", 'w'))
        cls.genome_features = {x['id']: x for x in cls.genome['features']}
        cls.genome_cdss = {x['id']: x for x in cls.genome['cdss']}

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
        self.assertIn('genome_tiers', genome)
        self.assertCountEqual(genome['genome_tiers'], ['ExternalDB'])
        self.assertTrue(genome.get("source") == "RefSeq", "Source is not RefSeq : " + str(genome.get("source")))
        self.assertTrue(len(genome["genome_tiers"]) == 1, "Should only have 1 tier in it.")

    def test_for_alias_colon(self):
        self.assertIn("ArthCp001_CDS_1", self.genome_cdss,"Did not find ArthCp001_CDS_1 in Genomes CDSs")
        cds = self.genome_cdss["ArthCp001_CDS_1"]
        colon_included = False
        for alias_tuple in cds["aliases"]:
            if alias_tuple[0] == "gene_synonym":
                if alias_tuple[1] == "TEST:COLON":
                    colon_included = True
        self.assertTrue(colon_included, "The synonym TEST:COLON was not found.")

    def test_for_db_xref_colon(self):
        self.assertIn("ArthCp001_CDS_1", self.genome_cdss,"Did not find ArthCp001_CDS_1 in Genomes CDSs")
        cds = self.genome_cdss["ArthCp001_CDS_1"]
        found_normal_colon = False
        found_1extra_colon = False
        found_2extra_colon = False
        for db_xref_tuple in cds["db_xrefs"]:
            if db_xref_tuple[0] == "MSI":
                if db_xref_tuple[1] == "123":
                    found_normal_colon = True
                elif db_xref_tuple[1] == "MSI:123":
                    found_1extra_colon = True
                elif db_xref_tuple[1] == "GB:MSI:123":
                    found_2extra_colon = True                                                    
        self.assertTrue(found_normal_colon, "The normal colon db_xref was not found.")
        self.assertTrue(found_1extra_colon, "The 1extra colon db_xref was not found.")
        self.assertTrue(found_2extra_colon, "The 2extra colon db_xref was not found.")

    def test_for_trans_splicing(self):
        genome = self.__class__.genome
        found_mRNA = False
        self.assertIn("ArthCp001", self.genome_features,"Did not find ArthCp001 in Genomes features")
        feature = self.genome_features["ArthCp001"]
        self.assertIn('trans_splicing', feature.get('flags', []),"The trans_splicing flag for the gene ArthCp001 was not found.")
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
        for feature in genome["mrnas"]:
            if feature['id'] == "ArthCp001_mRNA_1":
                found_mRNA = True
                if 'parent_gene' in feature:
                    self.assertTrue(feature['parent_gene'] == 'ArthCp001',"The parent gene for ArthCp001_CDS_1 was not as expected")
                else:
                    self.assertTrue('parent_gene' in feature, "The parent gene for ArthCp001_CDS_1 was not populated")
        self.assertTrue(found_mRNA, "The mRNA ArthCp001_mRNA_1 was not found.")
        self.assertIn("ArthCp001_CDS_1", self.genome_cdss,"Did not find ArthCp001_CDS_1 in Genomes CDSs")
        cds = self.genome_cdss["ArthCp001_CDS_1"]        
        self.assertIn('trans_splicing', cds.get('flags', []),"The trans_splicing flag for the CDS ArthCp001_CDS_1 was not found.")
        cds_sequence = "ATGCCAACCATTAAACAACTTATTAGAAATACAAGACAGCCAATCCGAAACGTCACGAAATCCCCAGCGCTTCGGGGATGCCCTCAGCGACGAGGAACATGTACTCGGGTGTATACTATCACCCCCAAAAAACCAAACTCTGCTTTACGTAAAGTTGCCAGAGTACGATTAACCTCGGGATTTGAAATCACTGCTTATATACCTGGTATTGGCCATAATTTACAAGAACATTCTGTAGTCTTAGTAAGAGGGGGAAGGGTTAAGGATTTACCCGGTGTGAGATATCACATTGTTCGAGGAACCCTAGATGCTGTCGGAGTAAAGGATCGTCAACAAGGGCGTTCTAAATATGGGGTCAAAAAGCCAAAATAA"
        self.assertTrue(cds['dna_sequence'] == cds_sequence, "The DNA sequence for the cds ArthCp001_CDS_1 was not as expected. It contained the following sequence : " + str(cds['dna_sequence']))
        cds_translation = "MPTIKQLIRNTRQPIRNVTKSPALRGCPQRRGTCTRVYTITPKKPNSALRKVARVRLTSGFEITAYIPGIGHNLQEHSVVLVRGGRVKDLPGVRYHIVRGTLDAVGVKDRQQGRSKYGVKKPK"
        self.assertTrue(cds['protein_translation'] == cds_translation, "The AA sequence for the cds ArthCp001_CDS_1 was not as expected. It contained the following sequence : " + str(cds['protein_translation']))
        if 'parent_gene' in cds:
            self.assertTrue(cds['parent_gene'] == 'ArthCp001',"The parent gene for ArthCp001_CDS_1 was not as expected")
        else:
            self.assertTrue('parent_gene' in cds, "The parent gene for ArthCp001_CDS_1 was not populated")

    def test_both_strand_trans_splicing(self):
        self.assertIn("ArthCp047", self.genome_features,"Did not find ArthCp047 in Genomes features")
        feature = self.genome_features["ArthCp047"]
        self.assertIn('trans_splicing', feature.get('flags', []),"The trans_splicing flag for the gene ArthCp047 was not found.")        
        gene_sequence = "ATGCCAACCATTAAACAACTTATTAGAAATACAAGACAGCCAATCCGAAACGTCACGAAATCCCCAGCGCTTCGGGGATGCCCTCAGCGACGAGGAACATGTACTCGGGTGTATACTATCACCCCCAAAAAACCAAACTCTGCTTTACGTAAAGTTGCCAGAGTACGATTAACCTCGGGATTTGAAATCACTGCTTATATACCTGGTATTGGCCATAATTTACAAGAACATTCTGTAGTCTTAGTAAGAGGGGGAAGGGTTAAGGATTTACCCGGTGTGAGATATCACATTGTTCGAGGAACCCTAGATGCTGTCGGAGTAAAGGATCGTCAACAAGGGCGTTCTAGTGCGTTGTAGATTCTTATCCAAGACTTGTATCATTTGATGATGCCATGTGAATCGCTAGAAACATGTGAAGTGTATGGCTAACCCAATAACGAAAGTTTCGTAAGGGGACTGGAGCAGGCTACCACGAGACAAAAGATCTTCTTTCAAAAGAGATTCGATTCGGAACTCTTATATGTCCAAGGTTCAATATTGAAATAATTTCAGAGGTTTTCCCTGACTTTGTCCGTGTCAACAAACAATTCGAAATGCCTCGACTTTTTTAGAACAGGTCCGGGTCAAATAGCAATGATTCGAAGCACTTATTTTTACACTATTTCGGAAACCCAAGGACTCAATCGTATGGATATGTAAAATACAGGATTTCCAATCCTAGCAGGAAAAGGAGGGAAACGGATACTCAATTTAAAAGTGAGTAAACAGAATTCCATACTCGATTTCAGAGATACATATATAATTCTGTGGAAAGCCGTATTCGATGAAAGTCGTATGTACGGTTTGGAGGGAGATCTTTCATATCTTTCGAGATCCACCCTACAATATGGGGTCAAAAAGCCAAAATAA"
        self.assertTrue(feature['dna_sequence'] == gene_sequence, "The DNA sequence for the gene ArthCp047 was not as expected. It contained the following sequence : " + str(feature['dna_sequence']))
        if 'cdss' in feature:
            self.assertTrue("ArthCp047_CDS_1" in feature['cdss'], "The child CDS of ArthCp047 was not found")
        else:
            self.assertTrue('cdss' in feature, "There was no child CDS for ArthCp047")
        self.assertIn("ArthCp047_CDS_1", self.genome_cdss,"Did not find ArthCp047_CDS_1 in Genomes CDSs")
        cds = self.genome_cdss["ArthCp047_CDS_1"]   
        self.assertIn('trans_splicing', cds.get('flags', []),"The trans_splicing flag for the gene ArthCp047_CDS_1 was not found.")      
        cds_sequence = "ATGCCAACCATTAAACAACTTATTAGAAATACAAGACAGCCAATCCGAAACGTCACGAAATCCCCAGCGCTTCGGGGATGCCCTCAGCGACGAGGAACATGTACTCGGGTGTATACTATCACCCCCAAAAAACCAAACTCTGCTTTACGTAAAGTTGCCAGAGTACGATTAACCTCGGGATTTGAAATCACTGCTTATATACCTGGTATTGGCCATAATTTACAAGAACATTCTGTAGTCTTAGTAAGAGGGGGAAGGGTTAAGGATTTACCCGGTGTGAGATATCACATTGTTCGAGGAACCCTAGATGCTGTCGGAGTAAAGGATCGTCAACAAGGGCGTTCTAAATATGGGGTCAAAAAGCCAAAATAA"
        self.assertTrue(cds['dna_sequence'] == cds_sequence, "The DNA sequence for the cds ArthCp047_CDS_1 was not as expected. It contained the following sequence : " + str(cds['dna_sequence']))
        cds_translation = "MPTIKQLIRNTRQPIRNVTKSPALRGCPQRRGTCTRVYTITPKKPNSALRKVARVRLTSGFEITAYIPGIGHNLQEHSVVLVRGGRVKDLPGVRYHIVRGTLDAVGVKDRQQGRSKYGVKKPK"
        self.assertTrue(cds['protein_translation'] == cds_translation, "The AA sequence for the cds ArthCp001_CDS_1 was not as expected. It contained the following sequence : " + str(cds['protein_translation']))
        if 'parent_gene' in cds:
            self.assertTrue(cds['parent_gene'] == 'ArthCp047',"The parent gene for ArthCp047_CDS_1 was not as expected")
        else:
            self.assertTrue('parent_gene' in cds, "The parent gene for ArthCp047_CDS_1 was not populated")

    def test_for_trans_splicing_multicontig(self):
        self.assertIn("MultiContigTransSpliced", self.genome_features,"Did not find MultiContigTransSpliced in Genomes features")
        feature = self.genome_features["MultiContigTransSpliced"]
        self.assertIn('trans_splicing', feature.get('flags', []),"The trans_splicing flag for the gene MultiContigTransSpliced was not found.")  
        gene_sequence = "ATGCCAACCATTAAACAACTTATTAGAAATACAAGACAGCCAATCCGAAACGTCACGAAATCCCCAGCGCTTCGGGGATGCCCTCAGCGACGAGGAACATGTACTCGGGTGTATACTATCACCCCCAAAAAACCAAACTCTGCTTTACGTAAAGTTGCCAGAGTACGATTAACCTCGGGATTTGAAATCACTGCTTATATACCTGGTATTGGCCATAATTTACAAGAACATTCTGTAGTCTTAGTAAGAGGGGGAAGGGTTAAGGATTTACCCGGTGTGAGATATCACATTGTTCGAGGAACCCTAGATGCTGTCGGAGTAAAGGATCGTCAACAAGGGCGTTCTAGTGCGTTGTAGATTCTTATCCAAGACTTGTATCATTTGATGATGCCATGTGAATCGCTAGAAACATGTGAAGTGTATGGCTAACCCAATAACGAAAGTTTCGTAAGGGGACTGGAGCAGGCTACCACGAGACAAAAGATCTTCTTTCAAAAGAGATTCGATTCGGAACTCTTATATGTCCAAGGTTCAATATTGAAATAATTTCAGAGGTTTTCCCTGACTTTGTCCGTGTCAACAAACAATTCGAAATGCCTCGACTTTTTTAGAACAGGTCCGGGTCAAATAGCAATGATTCGAAGCACTTATTTTTACACTATTTCGGAAACCCAAGGACTCAATCGTATGGATATGTAAAATACAGGATTTCCAATCCTAGCAGGAAAAGGAGGGAAACGGATACTCAATTTAAAAGTGAGTAAACAGAATTCCATACTCGATTTCAGAGATACATATATAATTCTGTGGAAAGCCGTATTCGATGAAAGTCGTATGTACGGTTTGGAGGGAGATCTTTCATATCTTTCGAGATCCACCCTACAATATGGGGTCAAAAAGCCAAAAATGTAG"
        self.assertTrue(feature['dna_sequence'] == gene_sequence, "The DNA sequence for the gene MultiContigTransSpliced was not as expected. It contained the following sequence : " + str(feature['dna_sequence']))
        if 'cdss' in feature:
            self.assertTrue("MultiContigTransSpliced_CDS_1" in feature['cdss'], "The child CDS of MultiContigTransSpliced was not found")
        else:
            self.assertTrue('cdss' in feature, "There was no child CDS for MultiContigTransSpliced")
        self.assertIn("MultiContigTransSpliced_CDS_1", self.genome_cdss,"Did not find MultiContigTransSpliced_CDS_1 in Genomes CDSs")
        cds = self.genome_cdss["MultiContigTransSpliced_CDS_1"]   
        self.assertIn('trans_splicing', cds.get('flags', []),"The trans_splicing flag for the gene MultiContigTransSpliced_CDS_1 was not found.")      
        cds_sequence = "ATGCCAACCATTAAACAACTTATTAGAAATACAAGACAGCCAATCCGAAACGTCACGAAATCCCCAGCGCTTCGGGGATGCCCTCAGCGACGAGGAACATGTACTCGGGTGTATACTATCACCCCCAAAAAACCAAACTCTGCTTTACGTAAAGTTGCCAGAGTACGATTAACCTCGGGATTTGAAATCACTGCTTATATACCTGGTATTGGCCATAATTTACAAGAACATTCTGTAGTCTTAGTAAGAGGGGGAAGGGTTAAGGATTTACCCGGTGTGAGATATCACATTGTTCGAGGAACCCTAGATGCTGTCGGAGTAAAGGATCGTCAACAAGGGCGTTCTAAATATGGGGTCAAAAAGCCAAAAATGTAG"
        self.assertTrue(cds['dna_sequence'] == cds_sequence, "The DNA sequence for the cds MultiContigTransSpliced_CDS_1 was not as expected. It contained the following sequence : " + str(cds['dna_sequence']))
        cds_translation = "MPTIKQLIRNTRQPIRNVTKSPALRGCPQRRGTCTRVYTITPKKPNSALRKVARVRLTSGFEITAYIPGIGHNLQEHSVVLVRGGRVKDLPGVRYHIVRGTLDAVGVKDRQQGRSKYGVKKPKM"
        self.assertTrue(cds['protein_translation'] == cds_translation, "The AA sequence for the cds MultiContigTransSpliced_CDS_1 was not as expected. It contained the following sequence : " + str(cds['protein_translation']))
        if 'parent_gene' in feature:
            self.assertTrue(cds['parent_gene'] == 'MultiContigTransSpliced',"The parent gene for MultiContigTransSpliced_CDS_1 was not as expected")
        else:
            self.assertTrue('parent_gene' in cds, "The parent gene for MultiContigTransSpliced_CDS_1 was not populated")

    def test_for_invalid_order(self):
        genome = self.genome
        self.assertIn('InvalidOrder', self.genome_features)
        gene = self.genome_features['InvalidOrder']
        self.assertNotIn('trans_splicing', gene.get('flags', []))
        self.assertIn(warnings['not_trans_spliced'], gene["warnings"])

        self.assertIn('InvalidOrder_CDS_1', self.genome_cdss)
        cds = self.genome_cdss['InvalidOrder_CDS_1']
        self.assertNotIn('trans_splicing', cds.get('flags', []))
        self.assertIn(warnings['not_trans_spliced'], cds["warnings"])

        self.assertTrue(genome.get("suspect") == 1,
                        "This genome has invalid position order features in it. It should be deemed suspect.")
        self.assertIn(warnings['genome_not_trans_spliced'].format(4), genome["warnings"])

    def test_for_zero_spanning_pos_strand_feature(self):
        self.assertIn('RL4742A', self.genome_features)
        gene = self.genome_features["RL4742A"]
        self.assertNotIn('trans_splicing', gene.get('flags', []), "The trans_splicing flag was set, technically it appears it may be transpliced, but the file does not state it to be.")
        self.assertNotIn(warnings['not_trans_spliced'], gene.get("warnings",[]), "The position coordinates are out of order, but they start and stop at the start and end of a circular contig, therefore it is valid.")

        self.assertIn('RL4742A_CDS_1', self.genome_cdss)
        cds = self.genome_cdss['RL4742A_CDS_1']
        self.assertNotIn('trans_splicing', cds.get('flags', []), "The trans_splicing flag was set, technically it appears it may be transpliced, but the file does not state it to be.")
        self.assertNotIn(warnings['not_trans_spliced'], cds.get("warnings",[]), "The trans_splicing flag for the cds RL4742A_CDS_1 was set, technically it appears it may be transpliced, but the file does not state it to be.")

    def test_for_zero_spanning_neg_strand_feature(self):
        self.assertIn('RL4742', self.genome_features)
        gene = self.genome_features["RL4742"]
        self.assertNotIn('trans_splicing', gene.get('flags', []), "The trans_splicing flag was set, technically it appears it may be transpliced, but the file does not state it to be.")
        self.assertNotIn(warnings['not_trans_spliced'], gene.get("warnings",[]), "The position coordinates are out of order, but they start and stop at the start and end of a circular contig, therefore it is valid.")

        self.assertIn('RL4742_CDS_1', self.genome_cdss)
        cds = self.genome_cdss['RL4742_CDS_1']
        self.assertNotIn('trans_splicing', cds.get('flags', []), "The trans_splicing flag was set, technically it appears it may be transpliced, but the file does not state it to be.")
        self.assertNotIn(warnings['not_trans_spliced'], cds.get("warnings",[]), "The trans_splicing flag for the cds RL4742_CDS_1 was set, technically it appears it may be transpliced, but the file does not state it to be.")

    def test_for_zero_spanning_two_exon_feature(self):
        self.assertIn('Zero_Span_two_exon', self.genome_features)
        gene = self.genome_features["Zero_Span_two_exon"]
        self.assertNotIn('trans_splicing', gene.get('flags', []), "The trans_splicing flag for the gene Zero_Span_two_exon was set, but it is not trans_spliced.")
        self.assertNotIn(warnings['not_trans_spliced'], gene.get("warnings",[]), "The position coordinates gene 'Zero_Span_two_exon' are out of order, but they start and stop at the start and end of a circular contig, therefore it is valid.")

        self.assertIn('Zero_Span_two_exon_CDS_1', self.genome_cdss)
        cds = self.genome_cdss['Zero_Span_two_exon_CDS_1']
        self.assertNotIn('trans_splicing', cds.get('flags', []), "The trans_splicing flag for the gene Zero_Span_two_exon_CDS_1 was set, but it is not trans_spliced.")
        self.assertNotIn(warnings['not_trans_spliced'], cds.get("warnings",[]), "The position coordinates gene 'Zero_Span_two_exon_CDS_1' are out of order, but they start and stop at the start and end of a circular contig, therefore it is valid.")

    def test_ensembl_ontology_terms(self):
        genome = self.__class__.genome
        self.assertIn('ArthCp001_CDS_1', self.genome_cdss)
        feature = self.genome_cdss['ArthCp001_CDS_1']
        self.assertTrue("ontology_terms" in feature,"There are 2 Ensembl style ontology terms that should be accounted for.")
        self.assertTrue("GO" in feature["ontology_terms"],"There is 1 Ensembl style ontology GO term that should be accounted for.")
        self.assertTrue("GO:0009523" in feature["ontology_terms"]["GO"],"GO:0009523 should be in the feature's ontology terms map ")
        self.assertTrue("GO:0009523" in genome["ontologies_present"]["GO"],"GO:0009523 should be in the ontologies_present map")
        self.assertTrue("PO" in feature["ontology_terms"],"There is 1 Ensembl style ontology PO term that should be accounted for.")
        self.assertTrue("PO:0000005" in feature["ontology_terms"]["PO"],"PO:0000005 should be in the feature's ontology terms map ")
        self.assertTrue("PO:0000005" in genome["ontologies_present"]["PO"],"PO:0000005 should be in the ontologies_present map")
        self.assertTrue("ontology_events" in genome,"There should be an ontology event for the upload.")

    def test_drastic_translation_length_difference(self):
        genome = self.__class__.genome
        self.assertIn('ArthCp004_CDS_1', self.genome_cdss)
        cds = self.genome_cdss['ArthCp004_CDS_1']
        self.assertIn(warnings["inconsistent_CDS_length"].format(240,88), cds.get("warnings",[]), "The AA translation length is significantly off in ArthCp004_CDS_1")
        self.assertIn(warnings["inconsistent_translation"], cds.get("warnings",[]), "The AA translation sequence is significantly off in ArthCp004_CDS_1")

        self.assertIn(warnings["genome_inc_CDS_length"].format("ArthCp004_CDS_1",240,88), genome.get("warnings",[]), "This does not have the translation length warning it should.")
        self.assertIn(warnings["genome_inc_translation"].format(2, 26), genome.get("warnings",[]), "This does not have the CDS translation proportion warning it should." + str(genome["warnings"]))
        self.assertTrue("suspect" in genome and genome["suspect"] == 1,"This has significant translation issues both in length and AA sequence.")

    def test_translation_not_supplied(self):
        self.assertIn('ArthCp015_CDS_1', self.genome_cdss)
        cds = self.genome_cdss['ArthCp015_CDS_1']
        self.assertIn(warnings["no_translation_supplied"], cds.get("warnings",[]), "Missing warning was derived from DNA Sequence.")
        self.assertTrue(cds.get("protein_translation"),"The translation was derived and populated.")  
        
    def test_translation_not_supplied_not_multiple_of_3(self):
        self.assertIn('ArthCp015_non3_CDS_1', self.genome_cdss)
        cds = self.genome_cdss['ArthCp015_non3_CDS_1']  
        self.assertIn(warnings["no_translation_supplied"] + "Sequence length 88 is not a multiple of three", cds.get("warnings",[]), "Missing warning unable to do translation.")  
        self.assertTrue(cds.get("protein_translation") == "","No translation supplied and not a multiple of 3, should be empty.")

    def test_ensembl_coordinates(self):
        self.assertIn('ArthCp004', self.genome_features)
        refseq_gene = self.genome_features["ArthCp004"]
        self.assertIn('ArthCp004_ENSEMBL', self.genome_features)
        ensembl_gene = self.genome_features["ArthCp004_ENSEMBL"]
        self.assertTrue(refseq_gene and ensembl_gene, "One or both of the following genes were not found: ArthCp004, ArthCp004_ENSEMBL")
        self.assertTrue(refseq_gene["location"] == ensembl_gene["location"],"The Ensembl style coordinates did not result in the same gene location information")
        self.assertTrue(refseq_gene["dna_sequence"] == ensembl_gene["dna_sequence"],"The Ensembl style coordinates did not result in the same gene sequence information")
 
        self.assertIn('ArthCp004_CDS_1', self.genome_cdss)
        refseq_CDS = self.genome_cdss['ArthCp004_CDS_1']  
        self.assertIn('ArthCp004_ENSEMBL_CDS_1', self.genome_cdss)
        ensembl_CDS = self.genome_cdss['ArthCp004_ENSEMBL_CDS_1'] 
        self.assertTrue(refseq_CDS and ensembl_CDS,"One or both of the following CDSs were not found: ArthCp004_CDS_1, ArthCp004_ENSEMBL_CDS_1")
        self.assertTrue(refseq_CDS["location"] == ensembl_CDS["location"],"The Ensembl style coordinates did not result in the same CDS location information: ENSEMBL CDS LOCATIONS : " +
                                                                        str(ensembl_CDS["location"]) + " --- REFSEQ CDS LOCATIONS : " + str(refseq_CDS["location"]))
        self.assertTrue(refseq_CDS["dna_sequence"] == ensembl_CDS["dna_sequence"], "The Ensembl style coordinates did not result in the same CDS sequence information")
        cds_translation = "MVKLRLKRCGRKQRAVYRILAIDVRYRREGRDLSKVGFYDPITNQTFLNLSAILDFLKKGAQPTRTAHDISKKAGIFTE"
        self.assertTrue(ensembl_CDS['protein_translation'] == cds_translation, "The AA sequence for the cds ArthCp004_ENSEMBL_CDS_1 was not as expected. It contained the following sequence : " + str(ensembl_CDS['protein_translation']))
 
    def test_check_ribosomal_slippage(self):
        self.assertIn('SPO1_87_CDS_1', self.genome_cdss)
        cds = self.genome_cdss['SPO1_87_CDS_1']  
        self.assertIn('ribosomal_slippage', cds.get('flags', []),"The CDS SPO1_87_CDS_1 was supposed to have a ribosomal slippage flag")      

    def test_invalid_coordinates_off_of_contig(self):
        genome = self.__class__.genome
        self.assertNotIn('ArthCp085', self.genome_features,"'ArthCp085' Should not have been added to the genes, since it has invalid coordinates over the end of the contig")
        self.assertNotIn('ArthCp085_CDS_1', self.genome_features,"'ArthCp085__CDS_1' Should not have been added to the CDS, since it has invalid coordinates over the end of the contig")
        found_gene_warning = False
        found_cds_warning = False
        if "warnings" in genome:
            for warning in genome["warnings"]:
                if warning == warnings["coordinates_off_end"].format('ArthCp085_gene'):
                    found_gene_warning = True
                if warning == warning == warnings["coordinates_off_end"].format('ArthCp085_CDS'):
                   found_cds_warning = True
            self.assertTrue(found_gene_warning,"SUSPECT: The warning for the invalid gene off the end of the contig was not found.")
            self.assertTrue(found_cds_warning,"SUSPECT: The warning for the invalid CDS off the end of the contig was not found.")

    def test_odd_coordinates_unknown_lower_bound(self):
        self.assertIn('ArthCp005', self.genome_features)
        self.assertIn(warnings["non_exact_coordinates"],
                      self.genome_features['ArthCp005']["warnings"])

        self.assertIn('ArthCp005_CDS_1', self.genome_cdss)
        self.assertIn(warnings["non_exact_coordinates"],
                      self.genome_cdss['ArthCp005_CDS_1']["warnings"])

    def test_odd_coordinates_unknown_upper_bound(self):
        self.assertIn('ArthCp006', self.genome_features)
        self.assertIn(warnings["non_exact_coordinates"],
                      self.genome_features['ArthCp006']["warnings"])

        self.assertIn('ArthCp006_CDS_1', self.genome_cdss)
        self.assertIn(warnings["non_exact_coordinates"],
                      self.genome_cdss['ArthCp006_CDS_1']["warnings"])

    def test_odd_coordinates_unknown_both_bounds(self):
        self.assertIn('ArthCp007', self.genome_features)
        self.assertIn(warnings["non_exact_coordinates"],
                      self.genome_features['ArthCp007']["warnings"])

        self.assertIn('ArthCp007_CDS_1', self.genome_cdss)
        self.assertIn(warnings["non_exact_coordinates"],
                      self.genome_cdss['ArthCp007_CDS_1']["warnings"])

    def test_long_misc_feature(self):
        genome = self.__class__.genome
        found_feature = False
        found_note = False
        has_sequence = False
        for feature in genome["non_coding_features"]:
            if feature["location"][0] == ['NC_000932.1', 10000, '+', 11001] and feature["type"] == "misc_feature":
                found_feature = True
                if "note" in feature and feature["note"] == "long_misc_feature":
                    found_note = True
                if "dna_sequence" in feature:
                    has_sequence = True
        self.assertTrue(found_feature,"Did not find long_misc_feature")
        self.assertTrue(found_note,"Did not find long_misc_feature note")
        self.assertFalse(has_sequence,"Too large should not have had sequence")

    def test_long_gene_feature(self):
        genome = self.__class__.genome
        found_feature = False
        found_note = False
        has_sequence = False
        for feature in genome["non_coding_features"]:
            if feature["location"][0] == ['NC_000932.1', 10000, '+', 11001] and feature["type"] == "gene":
                found_feature = True
                if "note" in feature and feature["note"] == "long_gene_feature":
                    found_note = True
                if "dna_sequence" in feature:
                    has_sequence = True
        self.assertTrue(found_feature,"Did not find long_misc_feature")
        self.assertTrue(found_note,"Did not find long_misc_feature note")
        self.assertTrue(has_sequence,"Is a gene should have sequence")

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


