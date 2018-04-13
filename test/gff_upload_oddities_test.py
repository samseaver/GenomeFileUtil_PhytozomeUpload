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
              'generate_ids_if_needed': 1,
              'strict': False
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

    def test_for_picking_up_other_fields(self):
        #tests for picking up ontologies. Other db_xrefs  (NOTE DB Xrefs different ways: in GenBank - /db_xref   in GFF I have seen: Dbxref=)
        #For picking up notes, function, product, locus_tag
        #THIS TEST NEEDS MORE WORK Check on ontologies and other stuff
        genome = self.__class__.genome
        found_rna6CDS_GO = False
        found_rna6CDS_PO = False
        found_rna5CDS_GO = False
        found_note = False
        found_function_product = False
        found_functional_descriptors = False
        found_locus_tag = False
        for feature in genome['cdss']:
            if feature["id"] == "rna6.CDS":
                #Do rna6.CDS first.
                #look for db_xrefs and ontologies.
                print "RNA6CDS: " + str(feature)
                #currently this is not picking up the ontologies in either form.
                if "ontology_terms" in feature:
                    print "GFF 6 ONTOLOGIES: " + str(feature["ontology_terms"])
                    if "GO" in feature["ontology_terms"]:
                        if feature["ontology_terms"]["GO"] == {"GO:0009523":[0]}:
                            found_rna6CDS_GO = True                            
                        print "GFF 6 G ONTOLOGY: " + str(feature["ontology_terms"]["GO"])
                    if "PO" in feature["ontology_terms"]:
                        if feature["ontology_terms"]["PO"] == {"PO:0000005":[0]}:
                            found_rna6CDS_PO = True  
                        print "GFF 6 P ONTOLOGY: " + str(feature["ontology_terms"]["PO"])                                        
            if feature["id"] == "rna5.CDS":
                #Do rna6.CDS first.
                #look for db_xrefs and ontologies.
                print "RNA5CDS: " + str(feature)
                if "ontology_terms" in feature:
                    print "GFF 5 ONTOLOGIES: " + str(feature["ontology_terms"])
                    if "GO" in feature["ontology_terms"]:
                        if feature["ontology_terms"]["GO"] == {"GO:0009523":[0]}:
                                found_rna5CDS_GO = True  
                        print "GFF 5 ONTOLOGY: " + str(feature["ontology_terms"]["GO"])
                if "note" in feature:
                    if feature["note"] == "Test":
                        found_note = True
                if "functional_descriptions" in feature:
                    if feature["functional_descriptions"][0] == "Test function":
                        found_functional_descriptors = True
                if "functions" in feature:
                    if feature["functions"][0] == "NAD-dependent sorbitol dehydrogenase":
                        found_function_product = True
                if "aliases" in feature:
                    for alias in feature["aliases"]:
                        if alias[0] == "locus_tag" and alias[1] == "Test_locus_tag":
                            found_locus_tag = True
        self.assertTrue(found_function_product)
        self.assertTrue(found_locus_tag)
        self.assertTrue(found_note)
        self.assertTrue(found_functional_descriptors)
        self.assertTrue(found_rna6CDS_GO)
        self.assertTrue(found_rna6CDS_PO)
        self.assertTrue(found_rna5CDS_GO)

    def test_off_contig(self):
        #test for feature off the end of the contig.
        #tested this it caused an error on upload, which is reasonable. So no object to actually test.
        #If change ID=gene117_off_contig; to have end coordinate of 422617 it fails.
        genome = self.__class__.genome
        found_gene = False
        for feature in genome["non_coding_features"]:
            if feature["id"] == 'gene117_off_contig':
                found_gene = True
        self.assertTrue(found_gene,"Did not find off contig gene")

    def test_2_contig_feature(self):
        #Test for a trans_spliced feature on 2 contigs
        genome = self.__class__.genome
        found_gene = False
        trans_splicing_flag = False
        for feature in genome["non_coding_features"]:
            if feature["id"] == 'gene_on_2contigs':
                found_gene = True
                #print "2Contig Feature: " + str(feature)
                self.assertTrue(len(feature["location"]) == 2, "@ contig gene does not have 2 locations.")
                self.assertTrue(feature["location"][0] == ['NC_010127.1', 13867, '+', 369], "First location incorrect.")
                self.assertTrue(feature["location"][1] == ['NC_010128.1', 6605, '+', 1959], "Second location incorrect.")
                if "flags" in feature:
                    if "trans_splicing" in feature["flags"]:
                        trans_splicing_flag = True      
        self.assertTrue(found_gene,"Did not find gene_on_2contigs gene")
        self.assertTrue(trans_splicing_flag, "gene trans_splicing flag not set.")

    def test_jgi_strand_plus(self):
        #Test JGI Plus strand, strand = "1"
        genome = self.__class__.genome
        found_gene = False
        for feature in genome["features"]:
            if feature["id"] == 'gene9':
                found_gene = True
                #print "JGI Plus Feature: " + str(feature)
                self.assertTrue(feature["location"][0] == ['NC_010127.1', 37271, '+', 960], "JGI Plus incorrect.")      
        self.assertTrue(found_gene,"Did not find JGI plus gene")

    def test_jgi_strand_minus(self):
        #Test JGI Minus strand, strand = "-1". Also check parent child for non-coding genes.
        #Also is testing "gene" with a child "transcript" and not "mRNA" is being put into "non_coding_features" list.
        #NOTE THIS IS FAILING NOW. THIS GENE CODES FOR A "transcript" feature type and not an "mRNA".
        #Thus it should be in non_coding_features. Currently being put into features.
        genome = self.__class__.genome
        found_gene = False
        found_transcript = False
        found_genes_transcript_child = False
#        for feature in genome["features"]:
#            if feature["id"] == 'gene8':
#                print "JGI Minus Feature in wrong place: " + str(feature)        
        for feature in genome["non_coding_features"]:
            if feature["id"] == 'gene8':
                found_gene = True
                #print "JGI Minus Feature: " + str(feature)
                self.assertTrue(feature["location"][0] == ['NC_010127.1', 37028, '-', 1238], "JGI Minus incorrect.") 
                if "children" in feature:
                    if "rna8" in feature["children"]:
                        found_genes_transcript_child = True
                self.assertTrue(found_genes_transcript_child,"The gene did not have the right child transcript.")
                self.assertTrue(feature.get("type") == "gene","The gene did not get the gene type")
            if feature["id"] == 'rna8':
                found_transcript = True
                self.assertTrue(feature.get("parent_gene") == "gene8","The transcript did not have the right gene parent")  
                self.assertTrue(feature.get("type") == "transcript","The transcript did not get the transcript type")                    
        self.assertTrue(found_gene,"Did not find JGI minus gene")
        self.assertTrue(found_transcript,"Did not find JGI minus transcript")

    def test_mRNA_child_fail_coordinate_validation(self):
        #Fail parent mRNA coordinate validation
        genome = self.__class__.genome
        found_mrna = False
        found_mrna_warning = False
        found_mrnas_cds = False
        found_CDS = False
        found_CDS_warning = False
        found_CDSs_mrna = False
        found_CDSs_gene = False
        found_CDS_gene_warning = False
        for feature in genome["mrnas"]:
            if feature["id"] == 'rna12':
                found_mrna = True
                #print "coordinate mRNA fail child Feature: " + str(feature)
                if "cds" in feature:
                    if feature["cds"] == "rna12.CDS":
                        #I think in terms of behavior the relationship should still exist even though the location validation fails.
                        #Since the relationship is defined in the file.  But should have a warning. 
                        found_mrnas_cds = True
                if "warnings" in feature:
                    if warnings["mRNA_fail_parent_coordinate_validation"].format('rna12.CDS') in feature["warnings"]:
                        found_mrna_warning = True
        for feature in genome["cdss"]:
            if feature["id"] == 'rna12.CDS':
                found_CDS = True
                #print "coordinate CDS fail parent Feature: " + str(feature)
                if "parent_mrna" in feature:
                    if feature["parent_mrna"] == "rna12":
                        #I think in terms of behavior the relationship should still exist even though the location validation fails.
                        #Since the relationship is defined in the file.  But should have a warning. 
                        found_CDSs_mrna = True
                if "parent_gene" in feature:
                    if feature["parent_gene"] == "gene12":
                        found_CDSs_gene = True
                if "warnings" in feature:
                    if warnings["CDS_fail_child_of_mRNA_coordinate_validation"].format('rna12') in feature["warnings"]:
                        found_CDS_warning = True
                    #CHECK IF PARENT GENE IS OK. NOTE IT PASSES
                    if warnings[ "CDS_fail_child_of_gene_coordinate_validation"].format('gene12') in feature["warnings"]:
                        found_CDS_gene_warning = True                                           
        self.assertTrue(found_mrna,"Did not find mRNA Feature")
        self.assertTrue(found_CDS,"Did not find CDS Feature")
        self.assertTrue(found_mrnas_cds,"Did not find mRNA's CDS Feature")
        self.assertTrue(found_CDSs_mrna,"Did not find CDS's mRNA Feature")
        self.assertTrue(found_CDSs_gene,"Did not find CDS's gene Feature")
        self.assertTrue(found_mrna_warning,"Did not find mRNA fail CDS coordinate validation warning")
        self.assertTrue(found_CDS_warning,"Did not find CDS fail mRNA coordinate validation warning")
        self.assertFalse(found_CDS_gene_warning,"Erroneous found warning for CDS failing coordinate validation to Gene")

    def test_genes_mRNA_child_fail_coordinate_validation(self):
        #Fail parent mRNA coordinate validation
        genome = self.__class__.genome
        found_gene = False
        found_genes_cds = False
        found_genes_mrna = False
        found_gene_warning = False
        found_mrna = False
        found_mrna_warning = False
        found_CDS = False
        found_CDS_warning = False
        for feature in genome["features"]:
            if feature["id"] == 'gene14':
                found_gene = True
                if "rna14.CDS" in feature["cdss"]:
                    found_genes_cds = True
                if "mrnas" in feature:
                    if "rna14" in feature["mrnas"]:
                        found_genes_mrna = True
                if "warnings" in feature:
                    if warnings["genes_mRNA_child_fails_location_validation"].format("rna14") in feature["warnings"]:
                        found_gene_warning = True
        for feature in genome["mrnas"]:
            if feature["id"] == 'rna14':
                found_mrna = True
                #print "coordinate mRNA fail Parent Gene: " + str(feature)
                self.assertTrue(feature.get("parent_gene") == 'gene14', "mRNA did not have the right parent gene")  
                self.assertTrue(feature.get("cds") == 'rna14.CDS', "mRNA did not have the right CDS")                 
                if "warnings" in feature:
                    if warnings["mRNAs_parent_gene_fails_location_validation"].format('gene14') in feature["warnings"]:
                        found_mrna_warning = True
        for feature in genome["cdss"]:
            if feature["id"] == 'rna14.CDS':
                found_CDS = True
                self.assertTrue(feature.get("parent_gene") == 'gene14', "CDS did not have the right parent gene")  
                self.assertTrue(feature.get("parent_mrna") == 'rna14', "CDS did not have the right parent mRNA")    
                #print "coordinate CDS fail parent Feature: " + str(feature)
                if "warnings" in feature:
                    if warnings["CDS_fail_child_coordinate_validation"].format('rna14') in feature["warnings"]:
                        found_CDS_warning = True
                #CHECK IF PARENT GENE IS OK.         
        self.assertTrue(found_gene,"Did not find Gene Feature")                    
        self.assertTrue(found_mrna,"Did not find mRNA Feature")
        self.assertTrue(found_CDS,"Did not find CDS Feature")
        self.assertTrue(found_genes_cds,"Did not find Gene's CDS") 
        self.assertTrue(found_genes_mrna,"Did not find Gene's mRNA")
        self.assertTrue(found_gene_warning,"Did not find gene fail mRNA coordinate validation warning")
        self.assertTrue(found_mrna_warning,"Did not find mRNA fail parnet coordinate validation warning")
        self.assertFalse(found_CDS_warning,"The CDS has bad warning. It is valid within the mRNA")

    def test_gene_child_CDS_fail_coordinate_validation(self):
        #CDS and gene relationship checking
        genome = self.__class__.genome
        found_gene = False
        found_gene_warning = False
        found_genes_cds = False
        found_CDS = False
        found_CDS_warning = False
        for feature in genome["features"]:
            if feature["id"] == 'gene21':
                found_gene = True
                #print "coordinate mRNA fail child Feature: " + str(feature)
                if "cdss" in feature:
                    if "gene21.CDS" in feature["cdss"]:
                        #I think in terms of behavior the relationship should still exist even though the location validation fails.
                        #Since the relationship is defined in the file.  But should have a warning. 
                        found_genes_cds = True
                if "warnings" in feature:
                    if warnings["genes_CDS_child_fails_location_validation"].format('gene21.CDS') in feature["warnings"]:
                        found_gene_warning = True
        for feature in genome["cdss"]:
            if feature["id"] == 'gene21.CDS':
                found_CDS = True
                #print "coordinate CDS fail parent Feature: " + str(feature)
                self.assertFalse(feature.get("parent_mrna",False),"This CDS is not supposed to have a parent mRNA")
                self.assertTrue(feature.get("parent_gene") == "gene21","This CDS has no parent Gene")
                if "warnings" in feature:
                    #CHECK IF PARENT GENE IS OK. NOTE IT PASSES
                    if warnings[ "CDS_fail_child_of_gene_coordinate_validation"].format('gene21') in feature["warnings"]:
                        found_CDS_warning = True                                           
        self.assertTrue(found_gene,"Did not find gene Feature")
        self.assertTrue(found_CDS,"Did not find CDS Feature")
        self.assertTrue(found_genes_cds,"Did not find Gene's CDS Feature")
        self.assertTrue(found_gene_warning,"Did not find Gene fail CDS coordinate validation warning")
        self.assertTrue(found_CDS_warning,"Did not find CDS fail mRNA coordinate validation warning")

    def test_sequences(self):
        #CHECK SEQUENCE AND PROTEIN TRANSLATION.
        genome = self.__class__.genome    
        found_cds = False
        for feature in genome['cdss']:
            if feature['id'] == "rna19.CDS":
                self.assertTrue(feature["protein_translation"] == "MMRRTLLSFWGLWCVIALVQSAFARPQAGTIQELLEAPASQYVE" +
                     "NIFVDGVVALSTAAGSVQAYEGTRNPLLNSVGSAALYFTLKKCAILEPHVHTNTPEFYFVISGTGTFSLWSSNGSVVHLKTPITNGSFM" +
                     "IIPAGWPHMITGPETSSTPLVLLANYLSGLPQVYFLASRASVFEQTSPSVMASVFNVSTAEYSTFFGAESGVGIVFNSSCLAS")
                self.assertTrue(feature["dna_sequence"] == "ATGATGCGCAGAACTTTGCTTTCGTTTTGGGGTTTGTGGTGCGTCATAGCACTCGTGCAGT" +
                    "CTGCGTTCGCGCGTCCACAAGCTGGTACGATCCAAGAGCTGCTGGAAGCTCCGGCTTCTCAGTACGTAGAGAATATCTTTGTCGACGGTGTGGTTGCACT" +
                    "CTCAACTGCGGCAGGGTCCGTACAGGCATATGAAGGGACACGGAATCCGTTGTTGAACTCTGTCGGCTCGGCTGCTTTGTATTTCACGTTGAAGAAGTGC" +
                    "GCGATTCTGGAGCCGCACGTTCACACCAATACACCCGAGTTCTATTTCGTCATTTCCGGTACGGGCACCTTCAGCCTTTGGTCCTCGAACGGCTCTGTTG" +
                    "TCCACCTGAAGACGCCCATAACCAACGGGAGCTTCATGATAATCCCGGCCGGGTGGCCGCACATGATCACCGGGCCGGAGACATCATCGACACCGTTGGT" +
                    "GCTCCTCGCCAACTATCTGAGCGGTCTTCCACAGGTCTATTTCCTTGCGAGCCGGGCAAGCGTCTTTGAACAAACCTCACCATCTGTGATGGCAAGCGTG" +
                    "TTCAATGTCTCGACAGCGGAGTACTCGACGTTCTTCGGCGCTGAATCCGGCGTCGGGATCGTCTTTAACTCGTCGTGTCTGGCGAGCTAG")
                found_cds = True
        self.assertTrue(found_cds,"Did not find CDS")
                
    def test_mRNA_coordinates_picking_up_exons_plus_strand(self):
        #accuracy of mRNA coodinates since mRNA is not taken directly from the GFF.
        genome = self.__class__.genome
        for feature in genome['mrnas']:
            if feature["id"] == "rna27a":
                self.assertTrue(feature["location"] == [['NC_010127.1', 101800, '+', 21], ['NC_010127.1', 101846, '+', 1655]])
                self.assertTrue(feature.get("cds") == "rna27a.CDS")
            if feature["id"] == "rna27b":
                self.assertTrue(feature["location"] == [['NC_010127.1', 101800, '+', 201], 
                                                        ['NC_010127.1', 102007, '+', 1423],
                                                        ['NC_010127.1', 103450, '+', 51]])
                self.assertTrue(feature.get("cds") == "rna27b.CDS")            
        for feature in genome['cdss']:
            if feature["id"] == "rna27a.CDS":
                self.assertTrue(feature["location"] == [['NC_010127.1', 101846, '+', 1584]],str(feature["location"]))
                self.assertTrue(feature.get("parent_mrna") == "rna27a")
            if feature["id"] == "rna27b.CDS":
                self.assertTrue(feature["location"] == [['NC_010127.1', 101846, '+', 155], ['NC_010127.1', 102007, '+', 1423]],
                                str(feature["location"]))
                self.assertTrue(feature.get("parent_mrna") == "rna27b")   

    def test_mRNA_coordinates_picking_up_exons_minus_strand(self):
        #accuracy of mRNA coodinates since mRNA is not taken directly from the GFF.
        genome = self.__class__.genome
        for feature in genome['mrnas']:
            if feature["id"] == "rna27a":
                self.assertTrue(feature["location"] == [['NC_010127.1', 101800, '+', 21], ['NC_010127.1', 101846, '+', 1655]])
                self.assertTrue(feature.get("cds") == "rna27a.CDS")
            if feature["id"] == "rna27b":
                self.assertTrue(feature["location"] == [['NC_010127.1', 101800, '+', 201], 
                                                        ['NC_010127.1', 102007, '+', 1423],
                                                        ['NC_010127.1', 103450, '+', 51]])
                self.assertTrue(feature.get("cds") == "rna27b.CDS")            
        for feature in genome['cdss']:
            if feature["id"] == "rna27a.CDS":
                self.assertTrue(feature["location"] == [['NC_010127.1', 101846, '+', 1584]],str(feature["location"]))
                self.assertTrue(feature.get("parent_mrna") == "rna27a")
            if feature["id"] == "rna27b.CDS":
                self.assertTrue(feature["location"] == [['NC_010127.1', 101846, '+', 155], ['NC_010127.1', 102007, '+', 1423]],
                                str(feature["location"]))
                self.assertTrue(feature.get("parent_mrna") == "rna27b")   

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
                        gene_trans_splicing_flag = True           
        for feature in genome['mrnas']:
            if feature['id'] == 'rna4ts':
                #print "RNA4ts Feature: " + str(feature)
                self.assertTrue(len(feature["location"]) == 3, "mRNA is not 3 locations as it should be.")
                self.assertTrue(feature['location'][0] == ['NC_010127.1', 69724, '-', 114])
                self.assertTrue(feature['location'][1] == ['NC_010127.1', 139856, '+', 232])   
                self.assertTrue(feature['location'][2] == ['NC_010127.1', 140625, '+', 26]) 
                if "flags" in feature:
                    if "trans_splicing" in feature["flags"]:
                        mRNA_trans_splicing_flag = True   
        for feature in genome['cdss']:
            if feature['id'] == 'rna4ts.CDS':
                #print "RNA4ts.CDS Feature: " + str(feature)
                self.assertTrue(len(feature["location"]) == 3, "CDS is not 3 locations as it should be.")
                if "flags" in feature:
                    if "trans_splicing" in feature["flags"]:
                        CDS_trans_splicing_flag = True
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
                        #print "FEATURE WITH A WARNING: " + str(feature)
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

    def test_entire_contigs_feature_types_not_included(self):
        #Looks at region, chromosome, scaffold and crazy full gene.
        genome = self.__class__.genome
        found_region = False
        found_chromosome = False
        found_scaffold = False
        found_crazy_gene = False
        found_crazy_gene_warning = False
        found_something_for_whole_contig = False
        for feature in genome["non_coding_features"]:
            if feature["id"] == "id0":
                found_region = True
            if feature["id"] == "id0_chromosome":
                found_chromosome = True
            if feature["id"] == "id0_scaffold":
                found_scaffold = True
            if feature['location'][0] == ['NC_010127.1', 1, '-', 422616]:
                found_something_for_whole_contig = True
        for feature in genome["non_coding_features"]:
            if feature["id"] == "id0_gene":
                found_crazy_gene = True
                if "warnings" in feature:
                    for warning in feature["warnings"]:
                        if warning == warnings["contig_length_feature"]:
                            found_crazy_gene_warning = True
        self.assertFalse(found_region,"This region should not have been included")
        self.assertFalse(found_chromosome,"This chromosome should not have been included")
        self.assertFalse(found_scaffold,"This scaffold should not have been included")
        self.assertFalse(found_something_for_whole_contig,"A feature for the whole contig should not have been included")
        self.assertTrue(found_crazy_gene)
        self.assertTrue(found_crazy_gene_warning)

    def test_odd_strands(self):
        #Testing cases where "." or "?" is used for the strand column.
        genome = self.__class__.genome
        found_dot_gene = False
        found_question_mRNA = False
        found_dot_made_plus = False
        found_question_made_plus = False
#not throwing a warning. So checking that it is defaulting it to "+" strand. Perhaps a warning in the future.
#        found_dot_warning = False
#        found_question_warning = False
        for feature in genome["features"]:
            if feature["id"] == "gene1":
                found_dot_gene = True
                if feature["location"][0][2] == "+":
                    found_dot_made_plus = True
#                if "warnings" in feature:
#                    for warning in feature["warnings"]:
#                        if warning == warnings["gff_odd_strand_type"].format("."):
#                            found_dot_warning = True
        for feature in genome["mrnas"]:
            if feature["id"] == "rna1":
                found_question_mRNA = True
                if feature["location"][0][2] == "+":
                    found_question_made_plus = True
#                if "warnings" in feature:
#                    for warning in feature["warnings"]:
#                        if warning == warnings["gff_odd_strand_type"].format("?"):
#                            found_question_warning = True 
        self.assertTrue(found_dot_gene) 
        self.assertTrue(found_question_mRNA) 
        self.assertTrue(found_dot_made_plus) 
        self.assertTrue(found_question_made_plus) 
#        self.assertTrue(found_dot_warning) 
#        self.assertTrue(found_question_warning) 
                            



