import unittest
import time
import os
import json

try:
    from ConfigParser import ConfigParser  # py2
except:
    from configparser import ConfigParser  # py3

from Workspace.WorkspaceClient import Workspace as workspaceService
from GenomeFileUtil.GenomeFileUtilImpl import GenomeFileUtil
from GenomeFileUtil.GenomeFileUtilServer import MethodContext
from DataFileUtil.DataFileUtilClient import DataFileUtil
from GenomeFileUtil.core.GenbankToGenome import warnings
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
        gbk_path = "data/Arabidopsis_gbff/A_thaliana_Ensembl_TAIR10_38_chr4_minus_xref.gbff"
        ws_obj_name = 'Yeast_chromosome1'
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
              'source': "Ensembl"
            })[0]
#        print("HERE IS THE RESULT:")
        data_file_cli = DataFileUtil(os.environ['SDK_CALLBACK_URL'], 
                                token=cls.ctx['token'],
                                service_ver='dev')
        cls.genome = data_file_cli.get_objects({'object_refs': [result['genome_ref']]})['data'][0]['data']
        json.dump(cls.genome,
                  open(cls.cfg['scratch'] + "/relationship_test_genome.json",
                       'w'))
#        print("GENE 1: ")
#        pprint(cls.genome['features'][0])
#        pprint(result)



    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def test_easy_1exon_relationship(self):
        #1 exon, 1 splice variant
        genome = self.__class__.genome
        found_gene = False
        found_CDS = False
        found_mRNA = False
        found_child_mRNA = False
        found_child_CDS = False
        found_sibling_mRNA = False
        found_sibling_CDS = False
        found_CDS_parent_gene = False
        found_mRNA_parent_gene = False
        for feature in genome["features"]:
            if feature['id'] == "AT4G00005":
                # print "Found AT4G00005  "  + str(feature)
                found_gene = True
                if "cdss" in feature:
                    for cds in feature["cdss"]:
                        if cds == "AT4G00005_CDS_1":
                            found_child_CDS = True
                if "mrnas" in feature:
                    for cds in feature["mrnas"]:
                        if cds == "AT4G00005_mRNA_1":
                            found_child_mRNA = True
        for feature in genome["mrnas"]:
            if feature['id'] == "AT4G00005_mRNA_1":
                found_mRNA = True
                if "cds" in feature:
                    if feature["cds"] == "AT4G00005_CDS_1":
                        found_sibling_CDS = True
                if "parent_gene" in feature:
                    if feature["parent_gene"] == "AT4G00005":
                        found_mRNA_parent_gene = True
        for feature in genome["cdss"]:
            if feature['id'] == "AT4G00005_CDS_1":
                found_CDS = True
                # print "AT4G00005_CDS_1 : " + str(feature)
                if "parent_mrna" in feature:
                    if feature["parent_mrna"] == "AT4G00005_mRNA_1":
                        found_sibling_mRNA = True
                if "parent_gene" in feature:
                    if feature["parent_gene"] == "AT4G00005":
                        found_CDS_parent_gene = True
        self.assertTrue(found_gene, "The gene AT4G00005 was not found.")
        self.assertTrue(found_mRNA, "The mRNA AT4G00005_mRNA_1 was not found.")
        self.assertTrue(found_CDS, "The CDS AT4G00005_CDS_1 was not found.")
        self.assertTrue(found_child_mRNA, "Did not find child_mRNA.")
        self.assertTrue(found_child_CDS, "Did not find child_CDS.")
        self.assertTrue(found_sibling_mRNA, "Did not find sibing_mRNA.")
        self.assertTrue(found_sibling_CDS, "Did not find sibling_CDS.")
        self.assertTrue(found_CDS_parent_gene, "Did not find CDS_Parent gene.")
        self.assertTrue(found_mRNA_parent_gene, "Did not find mRNA parent gene.")

    def test_easy_2variants_relationships(self):
        #Negative strand, 2 variants
        #Also played with inter mixing ensembl and refseq style coordinates.
        genome = self.__class__.genome
        found_gene = False
        found_CDS1 = False
        found_mRNA1 = False
        found_CDS2 = False
        found_mRNA2 = False
        found_child_mRNA1 = False
        found_child_CDS1 = False
        found_child_mRNA2 = False
        found_child_CDS2 = False
        found_sibling_mRNA1 = False
        found_sibling_CDS1 = False
        found_sibling_mRNA2 = False
        found_sibling_CDS2 = False
        found_CDS1_parent_gene = False
        found_mRNA1_parent_gene = False
        found_CDS2_parent_gene = False
        found_mRNA2_parent_gene = False
        for feature in genome["features"]:
            if feature['id'] == "AT4G00070":
                #print "Found AT4G00070  " + str(feature)
                found_gene = True
                if "cdss" in feature:
                    for cds in feature["cdss"]:
                        if cds == "AT4G00070_CDS_1":
                            found_child_CDS1 = True
                        if cds == "AT4G00070_CDS_2":
                            found_child_CDS2 = True
                if "mrnas" in feature:
                    for cds in feature["mrnas"]:
                        if cds == "AT4G00070_mRNA_1":
                            found_child_mRNA1 = True
                        if cds == "AT4G00070_mRNA_2":
                            found_child_mRNA2 = True
        for feature in genome["mrnas"]:
            if feature['id'] == "AT4G00070_mRNA_1":
                found_mRNA1 = True
                if "cds" in feature:
                    if feature["cds"] == "AT4G00070_CDS_1":
                        found_sibling_CDS1 = True
                if "parent_gene" in feature:
                    if feature["parent_gene"] == "AT4G00070":
                        found_mRNA1_parent_gene = True
            if feature['id'] == "AT4G00070_mRNA_2":
                found_mRNA2 = True
                if "cds" in feature:
                    if feature["cds"] == "AT4G00070_CDS_2":
                        found_sibling_CDS2 = True
                if "parent_gene" in feature:
                    if feature["parent_gene"] == "AT4G00070":
                        found_mRNA2_parent_gene = True
        for feature in genome["cdss"]:
            if feature['id'] == "AT4G00070_CDS_1":
                found_CDS1 = True
                if "parent_mrna" in feature:
                    if feature["parent_mrna"] == "AT4G00070_mRNA_1":
                        found_sibling_mRNA1 = True
                if "parent_gene" in feature:
                    if feature["parent_gene"] == "AT4G00070":
                        found_CDS1_parent_gene = True
            if feature['id'] == "AT4G00070_CDS_2":
                found_CDS2 = True
                if "parent_mrna" in feature:
                    if feature["parent_mrna"] == "AT4G00070_mRNA_2":
                        found_sibling_mRNA2 = True
                if "parent_gene" in feature:
                    if feature["parent_gene"] == "AT4G00070":
                        found_CDS2_parent_gene = True
        self.assertTrue(found_gene, "The gene AT4G00070 was not found.")
        self.assertTrue(found_mRNA1, "The mRNA AT4G00070_mRNA_1 was not found.")
        self.assertTrue(found_mRNA2, "The mRNA AT4G00070_mRNA_2 was not found.")
        self.assertTrue(found_CDS1, "The CDS AT4G00070_CDS_1 was not found.")
        self.assertTrue(found_CDS2, "The CDS AT4G00070_CDS_2 was not found.")
        self.assertTrue(found_child_mRNA1, "Did not find child_mRNA1.")
        self.assertTrue(found_child_mRNA2, "Did not find child_mRNA2.")
        self.assertTrue(found_child_CDS1, "Did not find child_CDS1.")
        self.assertTrue(found_child_CDS2, "Did not find child_CDS2.")
        self.assertTrue(found_sibling_mRNA1, "Did not find sibing_mRNA1.")
        self.assertTrue(found_sibling_mRNA2, "Did not find sibing_mRNA2.")
        self.assertTrue(found_sibling_CDS1, "Did not find sibling_CDS1.")
        self.assertTrue(found_sibling_CDS2, "Did not find sibling_CDS2.")
        self.assertTrue(found_CDS1_parent_gene, "Did not find CDS1_Parent gene.")
        self.assertTrue(found_CDS2_parent_gene, "Did not find CDS2_Parent gene.")
        self.assertTrue(found_mRNA1_parent_gene, "Did not find mRNA1 parent gene.")
        self.assertTrue(found_mRNA2_parent_gene, "Did not find mRNA2 parent gene.")


    def test_easy_3variants_relationships(self):
        #Positive strand, 3 variants
        genome = self.__class__.genome
        found_gene = False
        found_CDS1 = False
        found_mRNA1 = False
        found_CDS2 = False
        found_mRNA2 = False
        found_CDS3 = False
        found_mRNA3 = False
        found_child_mRNA1 = False
        found_child_CDS1 = False
        found_child_mRNA2 = False
        found_child_CDS2 = False
        found_child_mRNA3 = False
        found_child_CDS3 = False
        found_sibling_mRNA1 = False
        found_sibling_CDS1 = False
        found_sibling_mRNA2 = False
        found_sibling_CDS2 = False
        found_sibling_mRNA3 = False
        found_sibling_CDS3 = False
        found_CDS1_parent_gene = False
        found_mRNA1_parent_gene = False
        found_CDS2_parent_gene = False
        found_mRNA2_parent_gene = False
        found_CDS3_parent_gene = False
        found_mRNA3_parent_gene = False
        for feature in genome["features"]:
            if feature['id'] == "AT4G00050":
                #print "Found AT4G00050  " + str(feature)
                found_gene = True
                if "cdss" in feature:
                    for cds in feature["cdss"]:
                        if cds == "AT4G00050_CDS_1":
                            found_child_CDS1 = True
                        if cds == "AT4G00050_CDS_2":
                            found_child_CDS2 = True
                        if cds == "AT4G00050_CDS_3":
                            found_child_CDS3 = True
                if "mrnas" in feature:
                    for cds in feature["mrnas"]:
                        if cds == "AT4G00050_mRNA_1":
                            found_child_mRNA1 = True
                        if cds == "AT4G00050_mRNA_2":
                            found_child_mRNA2 = True
                        if cds == "AT4G00050_mRNA_3":
                            found_child_mRNA3 = True
        for feature in genome["mrnas"]:
            if feature['id'] == "AT4G00050_mRNA_1":
                found_mRNA1 = True
                if "cds" in feature:
                    if feature["cds"] == "AT4G00050_CDS_1":
                        found_sibling_CDS1 = True
                if "parent_gene" in feature:
                    if feature["parent_gene"] == "AT4G00050":
                        found_mRNA1_parent_gene = True
            if feature['id'] == "AT4G00050_mRNA_2":
                found_mRNA2 = True
                if "cds" in feature:
                    if feature["cds"] == "AT4G00050_CDS_2":
                        found_sibling_CDS2 = True
                if "parent_gene" in feature:
                    if feature["parent_gene"] == "AT4G00050":
                        found_mRNA2_parent_gene = True
            if feature['id'] == "AT4G00050_mRNA_3":
                found_mRNA3 = True
                if "cds" in feature:
                    if feature["cds"] == "AT4G00050_CDS_3":
                        found_sibling_CDS3 = True
                if "parent_gene" in feature:
                    if feature["parent_gene"] == "AT4G00050":
                        found_mRNA3_parent_gene = True
        for feature in genome["cdss"]:
            if feature['id'] == "AT4G00050_CDS_1":
                found_CDS1 = True
                if "parent_mrna" in feature:
                    if feature["parent_mrna"] == "AT4G00050_mRNA_1":
                        found_sibling_mRNA1 = True
                if "parent_gene" in feature:
                    if feature["parent_gene"] == "AT4G00050":
                        found_CDS1_parent_gene = True
            if feature['id'] == "AT4G00050_CDS_2":
                found_CDS2 = True
                if "parent_mrna" in feature:
                    if feature["parent_mrna"] == "AT4G00050_mRNA_2":
                        found_sibling_mRNA2 = True
                if "parent_gene" in feature:
                    if feature["parent_gene"] == "AT4G00050":
                        found_CDS2_parent_gene = True
            if feature['id'] == "AT4G00050_CDS_3":
                found_CDS3 = True
                if "parent_mrna" in feature:
                    if feature["parent_mrna"] == "AT4G00050_mRNA_3":
                        found_sibling_mRNA3 = True
                if "parent_gene" in feature:
                    if feature["parent_gene"] == "AT4G00050":
                        found_CDS3_parent_gene = True
        self.assertTrue(found_gene, "The gene AT4G00050 was not found.")
        self.assertTrue(found_mRNA1, "The mRNA AT4G00050_mRNA_1 was not found.")
        self.assertTrue(found_mRNA2, "The mRNA AT4G00050_mRNA_2 was not found.")
        self.assertTrue(found_mRNA3, "The mRNA AT4G00050_mRNA_3 was not found.")
        self.assertTrue(found_CDS1, "The CDS AT4G00050_CDS_1 was not found.")
        self.assertTrue(found_CDS2, "The CDS AT4G00050_CDS_2 was not found.")
        self.assertTrue(found_CDS3, "The CDS AT4G00050_CDS_3 was not found.")
        self.assertTrue(found_child_mRNA1, "Did not find child_mRNA1.")
        self.assertTrue(found_child_mRNA2, "Did not find child_mRNA2.")
        self.assertTrue(found_child_mRNA3, "Did not find child_mRNA3.")
        self.assertTrue(found_child_CDS1, "Did not find child_CDS1.")
        self.assertTrue(found_child_CDS2, "Did not find child_CDS2.")
        self.assertTrue(found_child_CDS3, "Did not find child_CDS3.")
        self.assertTrue(found_sibling_mRNA1, "Did not find sibing_mRNA1.")
        self.assertTrue(found_sibling_mRNA2, "Did not find sibing_mRNA2.")
        self.assertTrue(found_sibling_mRNA3, "Did not find sibing_mRNA3.")
        self.assertTrue(found_sibling_CDS1, "Did not find sibling_CDS1.")
        self.assertTrue(found_sibling_CDS2, "Did not find sibling_CDS2.")
        self.assertTrue(found_sibling_CDS3, "Did not find sibling_CDS3.")
        self.assertTrue(found_CDS1_parent_gene, "Did not find CDS1_Parent gene.")
        self.assertTrue(found_CDS2_parent_gene, "Did not find CDS2_Parent gene.")
        self.assertTrue(found_CDS3_parent_gene, "Did not find CDS3_Parent gene.")
        self.assertTrue(found_mRNA1_parent_gene, "Did not find mRNA1 parent gene.")
        self.assertTrue(found_mRNA2_parent_gene, "Did not find mRNA2 parent gene.")
        self.assertTrue(found_mRNA3_parent_gene, "Did not find mRNA3 parent gene.")

    def test_noncoding_gene_relationships(self):
        #Gene and misc_RNA
        genome = self.__class__.genome
        found_gene = False
        found_non_coding_gene = False
        found_misc_RNA = False
        found_child_RNA = False
        found_parent_gene = False
        for feature in genome["features"]:
            if feature['id'] == "AT4G00085":
                #print "Found AT4G00085  " + str(feature)
                found_gene = True
        for feature in genome["non_coding_features"]:
            if feature['id'] == "AT4G00085":
                #print "Found AT4G00085  " + str(feature)
                found_non_coding_gene = True        
                if "children" in feature:
                    for child in feature["children"]:
                        if child == "AT4G00085_misc_RNA_1":
                            found_child_RNA = True
            if feature['id'] == "AT4G00085_misc_RNA_1":
                found_misc_RNA= True        
                if "parent_gene" in feature:
                    if feature["parent_gene"] == "AT4G00085":
                        found_parent_gene = True
        self.assertFalse(found_gene, "The gene AT4G00085 is non-coding and should not be in the features list.")
        self.assertTrue(found_non_coding_gene, "The gene AT4G00085 was not found.")
        self.assertTrue(found_misc_RNA, "The misc_RNA AT4G00085_misc_RNA_1 was not found.")
        self.assertTrue(found_parent_gene, "Did not find AT4G00085_misc_RNA_1 parent gene.")
        self.assertTrue(found_child_RNA, "Did not find AT4G00085_misc_RNA_1 as a child of AT4G00085.")

    def test_easy_2variants_same_CDS(self):
        #Negative strand, 2 variants
        #2 different mRNA but bot result in same CDS
        genome = self.__class__.genome
        found_gene = False
        found_CDS1 = False
        found_mRNA1 = False
        found_CDS2 = False
        found_mRNA2 = False
        found_child_mRNA1 = False
        found_child_CDS1 = False
        found_child_mRNA2 = False
        found_child_CDS2 = False
        found_sibling_mRNA1 = False
        found_sibling_CDS1 = False
        found_sibling_mRNA2 = False
        found_sibling_CDS2 = False
        found_CDS1_parent_gene = False
        found_mRNA1_parent_gene = False
        found_CDS2_parent_gene = False
        found_mRNA2_parent_gene = False
        for feature in genome["features"]:
            if feature['id'] == "AT4G00165":
                found_gene = True
                if "cdss" in feature:
                    for cds in feature["cdss"]:
                        if cds == "AT4G00165_CDS_1":
                            found_child_CDS1 = True
                        if cds == "AT4G00165_CDS_2":
                            found_child_CDS2 = True
                if "mrnas" in feature:
                    for cds in feature["mrnas"]:
                        if cds == "AT4G00165_mRNA_1":
                            found_child_mRNA1 = True
                        if cds == "AT4G00165_mRNA_2":
                            found_child_mRNA2 = True
        for feature in genome["mrnas"]:
            if feature['id'] == "AT4G00165_mRNA_1":
                found_mRNA1 = True
                if "cds" in feature:
                    if feature["cds"] == "AT4G00165_CDS_1":
                        found_sibling_CDS1 = True
                if "parent_gene" in feature:
                    if feature["parent_gene"] == "AT4G00165":
                        found_mRNA1_parent_gene = True
            if feature['id'] == "AT4G00165_mRNA_2":
                found_mRNA2 = True
                if "cds" in feature:
                    if feature["cds"] == "AT4G00165_CDS_2":
                        found_sibling_CDS2 = True
                if "parent_gene" in feature:
                    if feature["parent_gene"] == "AT4G00165":
                        found_mRNA2_parent_gene = True
        for feature in genome["cdss"]:
            if feature['id'] == "AT4G00165_CDS_1":
                found_CDS1 = True
                if "parent_mrna" in feature:
                    if feature["parent_mrna"] == "AT4G00165_mRNA_1":
                        found_sibling_mRNA1 = True
                if "parent_gene" in feature:
                    if feature["parent_gene"] == "AT4G00165":
                        found_CDS1_parent_gene = True
            if feature['id'] == "AT4G00165_CDS_2":
                found_CDS2 = True
                if "parent_mrna" in feature:
                    if feature["parent_mrna"] == "AT4G00165_mRNA_2":
                        found_sibling_mRNA2 = True
                if "parent_gene" in feature:
                    if feature["parent_gene"] == "AT4G00165":
                        found_CDS2_parent_gene = True
        self.assertTrue(found_gene, "The gene AT4G00165 was not found.")
        self.assertTrue(found_mRNA1, "The mRNA AT4G00165_mRNA_1 was not found.")
        self.assertTrue(found_mRNA2, "The mRNA AT4G00165_mRNA_2 was not found.")
        self.assertTrue(found_CDS1, "The CDS AT4G00165_CDS_1 was not found.")
        self.assertTrue(found_CDS2, "The CDS AT4G00165_CDS_2 was not found.")
        self.assertTrue(found_child_mRNA1, "Did not find child_mRNA1.")
        self.assertTrue(found_child_mRNA2, "Did not find child_mRNA2.")
        self.assertTrue(found_child_CDS1, "Did not find child_CDS1.")
        self.assertTrue(found_child_CDS2, "Did not find child_CDS2.")
        self.assertTrue(found_sibling_mRNA1, "Did not find sibing_mRNA1.")
        self.assertTrue(found_sibling_mRNA2, "Did not find sibing_mRNA2.")
        self.assertTrue(found_sibling_CDS1, "Did not find sibling_CDS1.")
        self.assertTrue(found_sibling_CDS2, "Did not find sibling_CDS2.")
        self.assertTrue(found_CDS1_parent_gene, "Did not find CDS1_Parent gene.")
        self.assertTrue(found_CDS2_parent_gene, "Did not find CDS2_Parent gene.")
        self.assertTrue(found_mRNA1_parent_gene, "Did not find mRNA1 parent gene.")
        self.assertTrue(found_mRNA2_parent_gene, "Did not find mRNA2 parent gene.")

    def test_CDS_not_inside_gene(self):
        genome = self.__class__.genome
        found_gene_in_features = False
        found_gene_in_non_coding = False
        found_mRNA = False
        found_CDS = False
        genome_suspect = False
        found_warning = False
        for feature in genome["features"]:
            if feature['id'] == "AT4G12617":
                found_gene_in_features = True
        for feature in genome["non_coding_features"]:
            if feature['id'] == "AT4G12617":
                found_gene_in_non_coding = True
        for feature in genome["mrnas"]:
            if feature['id'] == "AT4G12617_mRNA_1":
                found_mRNA = True
        for feature in genome["cdss"]:
            if feature['id'] == "AT4G12617_CDS_1":
                found_CDS = True
        self.assertFalse(found_gene_in_features, "The gene AT4G12617 was found in features, it should have been excluded.")
        self.assertFalse(found_gene_in_non_coding, "The gene AT4G12617 was found in non_coding, it should have been excluded.")
        self.assertFalse(found_mRNA, "The mRNA AT4G12617_mRNA_1 was found, it should have been excluded.")
        self.assertFalse(found_CDS, "The CDS AT4G12617_CDS_1 was found, it should have been excluded.")
        if "suspect" in genome:
            if genome["suspect"] == 1:
                genome_suspect = True
        self.assertTrue(genome_suspect, "The Genome should have been made suspect because the CDS is not within the gene.")
        if "warnings" in genome:
            for warning in genome["warnings"]:
                if warning == warnings['gene_excluded'].format("AT4G12617"):
                    found_warning = True
        self.assertTrue(found_warning, "The Genome level warning for the failed coordinates matching.")

    def test_mRNA_not_inside_gene(self):
        genome = self.__class__.genome
        found_gene_in_features = False
        found_gene_in_non_coding = False
        found_mRNA = False
        found_CDS = False
        genome_suspect = False
        found_gene_warning = False
        found_warning = False
        for feature in genome["features"]:
            if feature['id'] == "AT4G12620":
                found_gene_in_features = True
                if "warnings" in feature:
                    for warning in feature["warnings"]:
                        if warning == "The child mRNA failed location validation. That mRNA has been excluded.":
                            found_gene_warning = True
        for feature in genome["non_coding_features"]:
            if feature['id'] == "AT4G12620":
                found_gene_in_non_coding = True
        for feature in genome["mrnas"]:
            if feature['id'] == "AT4G12620_mRNA_1":
                found_mRNA = True
        for feature in genome["cdss"]:
            if feature['id'] == "AT4G12620_CDS_1":
                found_CDS = True
        self.assertTrue(found_gene_in_features, "The gene AT4G12620 was not found in features")
        self.assertFalse(found_gene_in_non_coding, "The gene AT4G12620 was found in non_coding.")
        self.assertFalse(found_mRNA, "The mRNA AT4G12620_mRNA_1 was found, it should have been excluded.")
        self.assertTrue(found_CDS, "The CDS AT4G12620_CDS_1 was not found.")
        self.assertTrue(found_gene_warning,"Gene warning for failing mRNA not found.")
        if "suspect" in genome:
            if genome["suspect"] == 1:
                genome_suspect = True
        self.assertTrue(genome_suspect, "The Genome should be suspect because the mRNA is not within the gene.")
        if "warnings" in genome:
            for warning in genome["warnings"]:
                if warning == "SUSPECT mRNA AT4G12620 was excluded because the associated mRNA failed coordinates validation":
                    found_warning = True
        self.assertTrue(found_warning, "The Genome level warning for the failed coordinates matching.")

    def test_CDS_and_mRNA_not_inside_gene(self):
        genome = self.__class__.genome
        found_gene_in_features = False
        found_gene_in_non_coding = False
        found_mRNA = False
        found_CDS = False
        genome_suspect = False
        found_warning = False
        for feature in genome["features"]:
            if feature['id'] == "AT4G12485":
                found_gene_in_features = True
        for feature in genome["non_coding_features"]:
            if feature['id'] == "AT4G12485":
                found_gene_in_non_coding = True
        for feature in genome["mrnas"]:
            if feature['id'] == "AT4G12485_mRNA_1":
                found_mRNA = True
        for feature in genome["cdss"]:
            if feature['id'] == "AT4G12485_CDS_1":
                found_CDS = True
        self.assertFalse(found_gene_in_features, "The gene AT4G12485 was found in features, it should have been excluded.")
        self.assertFalse(found_gene_in_non_coding, "The gene AT4G12485 was found in non_coding, it should have been excluded.")
        self.assertFalse(found_mRNA, "The mRNA AT4G12485_mRNA_1 was found, it should have been excluded.")
        self.assertFalse(found_CDS, "The CDS AT4G12485_CDS_1 was found, it should have been excluded.")
        if "suspect" in genome:
            if genome["suspect"] == 1:
                genome_suspect = True
        self.assertTrue(genome_suspect, "The Genome should be suspect because the CDS is not within the gene.")
        if "warnings" in genome:
            for warning in genome["warnings"]:
                if warning == warnings['gene_excluded'].format("AT4G12485"):
                    found_warning = True
        self.assertTrue(found_warning, "The Genome level warning for the failed coordinates matching.")

    def test_CDS_not_inside_mRNA(self):
        #Thoughts on behavior welcome. 
        genome = self.__class__.genome
        found_gene_in_features = False
        found_gene_in_non_coding = False
        found_mRNA = False
        found_CDS = False
        CDS_has_mRNA = False
        mRNA_has_CDS = False
        gene_has_CDS = False
        gene_has_mRNA = False
        mRNA_has_gene = False
        CDS_has_gene = False
        mRNA_warning = False
        CDS_warning = False
        for feature in genome["features"]:
            if feature['id'] == "AT4G12490":
                found_gene_in_features = True
                if "cdss" in feature:
                    if feature["cdss"][0] == "AT4G12490_CDS_1":
                        gene_has_CDS = True
                if "mrnas" in feature:
                    if feature["mrnas"][0] == "AT4G12490_mRNA_1":
                        gene_has_mRNA = True                    
        for feature in genome["non_coding_features"]:
            if feature['id'] == "AT4G12490":
                found_gene_in_non_coding = True
        for feature in genome["mrnas"]:
            if feature['id'] == "AT4G12490_mRNA_1":
                found_mRNA = True
                if "cds" in feature:
                    mRNA_has_CDS = True
                if "parent_gene" in feature:
                    mRNA_has_gene = True
                if "warnings" in feature:
                    for warning in feature["warnings"]:
                        if warning == "Potential child CDS relationship failed due to location validation.":
                            mRNA_warning = True
        for feature in genome["cdss"]:
            if feature['id'] == "AT4G12490_CDS_1":
                found_CDS = True
                if "mrna" in feature:
                    CDS_has_mRNA = True
                if "parent_gene" in feature:
                    CDS_has_gene = True 
                if "warnings" in feature:
                    for warning in feature["warnings"]:
                        if warning == warnings['cds_mrna_cds'].format('AT4G12490_mRNA_1'):
                            CDS_warning = True               
        self.assertTrue(found_gene_in_features, "The gene AT4G12490 was not found in features")
        self.assertFalse(found_gene_in_non_coding, "The gene AT4G12490 was found in non_coding, it should be in features.")
        self.assertTrue(found_mRNA, "The mRNA AT4G12490_mRNA_1 was not found.")
        self.assertTrue(found_CDS, "The CDS AT4G12490_CDS_1 was not found.")
        self.assertTrue(gene_has_CDS, "The gene's CDS relationship was not found.")
        self.assertTrue(gene_has_mRNA, "The gene's mRNA relationship was not found.")
        self.assertFalse(mRNA_has_CDS, "The mRNA's CDS relationship should not exist.")
        self.assertTrue(mRNA_has_gene, "The mRNA's gene relationship was not found.")
        self.assertTrue(mRNA_warning, "The mRNA's warning was not found.")
        self.assertFalse(CDS_has_mRNA, "The CDS's mRNA relationship should not exist.")
        self.assertTrue(CDS_has_gene, "The CDS's gene relationship was not found.")
        self.assertTrue(CDS_warning, "The CDS's warning was not found.")        

    def test_2variants_1CDS_not_inside_gene(self):
        #Gene with 2 variants. 1 Variant passes location checking, the other one fails.
        genome = self.__class__.genome
        found_gene = False
        found_mRNA1 = False
        found_mRNA2 = False
        found_CDS1 = False
        found_CDS2 = False
        gene_has_CDS1 = False
        gene_has_CDS2 = False
        gene_has_mRNA1 = False
        gene_has_mRNA2 = False
        found_gene_mRNA_warning = False
        found_gene_CDS_warning = False
        found_mRNA1_parent = False
        found_CDS1_parent = False
        found_mRNA1_CDS1 = False
        found_CDS1_mRNA1 = False
        found_genome_warning = False
        for feature in genome["features"]:
            if feature['id'] == "AT4G12580":
                found_gene = True
                if "cdss" in feature:
                    if feature["cdss"][0] == "AT4G12580_CDS_1":
                        gene_has_CDS1 = True
                    if feature["cdss"][0] == "AT4G12580_CDS_2":
                        gene_has_CDS2 = True
                if "warnings" in feature:
                    for warning in feature["warnings"]:
                        if warning == "The child mRNA failed location validation. That mRNA has been excluded.":
                            found_gene_mRNA_warning = True
                        if warning == "The child CDS failed location validation. That CDS has been excluded.":
                                found_gene_CDS_warning = True
                if "mrnas" in feature:
                    if feature["mrnas"][0] == "AT4G12580_mRNA_1":
                        gene_has_mRNA1 = True   
                    if feature["mrnas"][0] == "AT4G12580_mRNA_2":
                        gene_has_mRNA2 = True            
        for feature in genome["mrnas"]:
            if feature['id'] == "AT4G12580_mRNA_1":
                found_mRNA1 = True
                if feature["parent_gene"] == "AT4G12580":
                    found_mRNA1_parent = True
                if "cds" in feature:
                    if feature["cds"] == "AT4G12580_CDS_1":
                        found_mRNA1_CDS1 = True
            if feature['id'] == "AT4G12580_mRNA_2":
                found_mRNA2 = True
        for feature in genome["cdss"]:
            if feature['id'] == "AT4G12580_CDS_1":
                found_CDS1 = True
                if feature["parent_gene"] == "AT4G12580":
                    found_CDS1_parent = True
                if "parent_mrna" in feature:
                    if feature["parent_mrna"] == "AT4G12580_mRNA_1":
                        found_CDS1_mRNA1 = True
            if feature['id'] == "AT4G12580_CDS_2":
                found_CDS2 = True
        self.assertTrue(found_gene, "The gene AT4G12580 was not found in features.")
        self.assertTrue(found_mRNA1, "The mRNA AT4G12580_mRNA_1 was not found.")
        self.assertFalse(found_mRNA2, "The mRNA AT4G12580_mRNA_2 was found, it should have been excluded.")
        self.assertTrue(found_CDS1, "The CDS AT4G12580_CDS_1 was not found.")
        self.assertFalse(found_CDS2, "The CDS AT4G12580_CDS_2 was found, it should have been excluded.")
        self.assertTrue(gene_has_CDS1,"The gene did not have the good CDS1")
        self.assertFalse(gene_has_CDS2, "The gene had CDS AT4G12580_CDS_2 was found, it should have been excluded.")
        self.assertTrue(gene_has_mRNA1,"The gene did not have the good mRNA1.")
        self.assertFalse(gene_has_mRNA2, "The gene had mRNA AT4G12580_mRNA_2 was found, it should have been excluded.")
        self.assertTrue(found_gene_mRNA_warning,"The gene did not have the mRNA related warning.")
        self.assertTrue(found_gene_CDS_warning,"The gene did not have the CDS related warning.")
        self.assertTrue(found_mRNA1_parent,"The mRNA did not have the parent gene.")
        self.assertTrue(found_CDS1_parent,"The CDS did not have the parent gene.")
        self.assertTrue(found_mRNA1_CDS1,"The mRNA did not have the corresponding CDS.")
        self.assertTrue(found_CDS1_mRNA1,"The CDS did not have the correspondig mRNA")
        self.assertTrue(found_genome_warning,"The proper genome level wearning was not found.")
        if "warnings" in genome:
            for warning in genome["warnings"]:
                if warning == "SUSPECT gene AT4G12580 had some of its child features (CDS and/or mRNAs) excuded because of failed coordinates validation":
                    found_genome_warning = True
        self.assertTrue(found_genome_warning, "The Genome level warning for the failed coordinates matching.")

    def test_CDS_not_sharing_mRNA_internal_boundaries(self):
        #CDS not sharing internal boundaries with mRNA.
        genome = self.__class__.genome
        found_gene = False
        found_mRNA = False
        found_CDS = False
        gene_has_CDS = False
        gene_has_mRNA = False
        found_mRNA_warning = False
        found_CDS_warning = False
        found_mRNA_parent = False
        found_CDS_parent = False
        found_mRNA_CDS = False
        found_CDS_mRNA = False
        for feature in genome["features"]:
            if feature['id'] == "AT4G12560":
                found_gene = True
                if "cdss" in feature:
                    if feature["cdss"][0] == "AT4G12560_CDS_1":
                        gene_has_CDS = True
                if "mrnas" in feature:
                    if feature["mrnas"][0] == "AT4G12560_mRNA_1":
                        gene_has_mRNA = True
        for feature in genome["mrnas"]:
            if feature['id'] == "AT4G12560_mRNA_1":
                found_mRNA = True
                if feature["parent_gene"] == "AT4G12560":
                    found_mRNA_parent = True
            if "cds" in feature:
                if feature["cds"] == "AT4G12560_CDS_1":
                    found_mRNA_CDS = True
            if "warnings" in feature:
                if "Potential child CDS relationship failed due to location validation." in feature["warnings"]:
                    found_mRNA_warning = True
        for feature in genome["cdss"]:
            if feature['id'] == "AT4G12560_CDS_1":
                found_mRNA = True
                if feature["parent_gene"] == "AT4G12560":
                    found_CDS_parent = True
            if "parent_mrna" in feature:
                if feature["parent_mrna"] == "AT4G12560_mRNA_1":
                    found_CDS_mRNA = True
            if "warnings" in feature:
                if warnings['cds_mrna_cds'].format('AT4G12560_mRNA_1') in feature["warnings"]:
                    found_CDS_warning = True
        self.assertTrue(found_gene, "The gene AT4G12560 was not found in features.")
        self.assertTrue(found_mRNA, "The mRNA AT4G12560_mRNA_1 was not found.")
        self.assertTrue(found_CDS, "The CDS AT4G12560_CDS_1 was not found.")
        self.assertTrue(gene_has_CDS,"The gene did not have the good CDS1")
        self.assertTrue(gene_has_mRNA,"The gene did not have the good mRNA1.")
        self.assertTrue(found_mRNA_warning,"No mRNA related warning.")
        self.assertTrue(found_CDS_warning,"No CDS related warning.")
        self.assertTrue(found_mRNA_parent,"The mRNA did not have the parent gene.")
        self.assertTrue(found_CDS_parent,"The CDS did not have the parent gene.")
        self.assertFalse(found_mRNA_CDS,"The mRNA should not have had a corresponding CDS.")
        self.assertFalse(found_CDS_mRNA,"The CDS should not have had a corresponding mRNA.") 
            
    def test_CDS_minus1_internal_exon_mRNA(self):
        #CDS has 1 less internal exon than the parent mRNA. 
        #Thus the CDS sequence will not be in the mRNA sequence.
        genome = self.__class__.genome
        found_gene = False
        found_mRNA = False
        found_CDS = False
        gene_has_CDS = False
        gene_has_mRNA = False
        found_mRNA_warning = False
        found_CDS_warning = False
        found_mRNA_parent = False
        found_CDS_parent = False
        found_mRNA_CDS = False
        found_CDS_mRNA = False
        for feature in genome["features"]:
            if feature['id'] == "AT4G12600":
                found_gene = True
                if "cdss" in feature:
                    if feature["cdss"][0] == "AT4G12600_CDS_1":
                        gene_has_CDS = True
                if "mrnas" in feature:
                    if feature["mrnas"][0] == "AT4G12600_mRNA_1":
                        gene_has_mRNA = True
        for feature in genome["mrnas"]:
            if feature['id'] == "AT4G12560_mRNA_1":
                found_mRNA = True
                if feature["parent_gene"] == "AT4G12600":
                    found_mRNA_parent = True
            if "cds" in feature:
                if feature["cds"] == "AT4G12600_CDS_1":
                    found_mRNA_CDS = True
            if "warnings" in feature:
                if "Potential child CDS relationship failed due to location validation." in feature["warnings"]:
                    found_mRNA_warning = True
        for feature in genome["cdss"]:
            if feature['id'] == "AT4G12600_CDS_1":
                found_mRNA = True
                if feature["parent_gene"] == "AT4G12600":
                    found_CDS_parent = True
            if "parent_mrna" in feature:
                if feature["parent_mrna"] == "AT4G12600_mRNA_1":
                    found_CDS_mRNA = True
            if "warnings" in feature:
                if warnings['cds_mrna_cds'].format('AT4G12600_mRNA_1') in feature["warnings"]:
                    found_CDS_warning = True
        self.assertTrue(found_gene, "The gene AT4G12600 was not found in features.")
        self.assertTrue(found_mRNA, "The mRNA AT4G12600_mRNA_1 was not found.")
        self.assertTrue(found_CDS, "The CDS AT4G12600_CDS_1 was not found.")
        self.assertTrue(gene_has_CDS,"The gene did not have the good CDS1")
        self.assertTrue(gene_has_mRNA,"The gene did not have the good mRNA1.")
        self.assertTrue(found_mRNA_warning,"No mRNA related warning.")
        self.assertTrue(found_CDS_warning,"No CDS related warning.")
        self.assertTrue(found_mRNA_parent,"The mRNA did not have the parent gene.")
        self.assertTrue(found_CDS_parent,"The CDS did not have the parent gene.")
        self.assertFalse(found_mRNA_CDS,"The mRNA should not have had a corresponding CDS.")
        self.assertFalse(found_CDS_mRNA,"The CDS should not have had a corresponding mRNA.") 
                  

'''
TO DO
Same CDS valid to multiple mRNA and CDS (happy still, valid) - DONE
CDS not a child of gene by location
mRNA not a child of gene by location
CDS and mRNA not a child of gene (but CDS child of mRNA)
CDS not a child of mRNA
Gene has 2 variants one with good CDS, one with invalid CDS
internal boundaries of CDS with mRNA
CDS with mRNA parent with 1 more exon 
CDS seq not in mRNA sequence.

transpliced - regular, multiple strands, multiple contigs
'''