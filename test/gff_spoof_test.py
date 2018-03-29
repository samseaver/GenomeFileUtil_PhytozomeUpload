import unittest
import time
import os

try:
    from ConfigParser import ConfigParser  # py2
except:
    from configparser import ConfigParser  # py3

from Workspace.WorkspaceClient import Workspace as workspaceService
from GenomeFileUtil.GenomeFileUtilImpl import GenomeFileUtil
from GenomeFileUtil.GenomeFileUtilServer import MethodContext
from DataFileUtil.DataFileUtilClient import DataFileUtil
from GenomeFileUtil.core.GenomeUtils import warnings


class GenomeFileUtilTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        token = os.environ.get('KB_AUTH_TOKEN', None)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx['token'] = token
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
        

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        if hasattr(self.__class__, 'wsName'):
            return self.__class__.wsName
        suffix = int(time.time() * 1000)
        wsName = "test_GenomeFileUtil_" + str(suffix)
        self.getWsClient().create_workspace({'workspace': wsName})
        self.__class__.wsName = wsName
        return wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    def test_spoof_on(self):
        gff_path = "data/fasta_gff/JGI/Bacterial_Data/91705.assembled.gff"
        fna_path = "data/fasta_gff/JGI/Bacterial_Data/91705.assembled.fna"
        ws_obj_name = 'JGI_Spoof'
        result = self.getImpl().fasta_gff_to_genome(
            self.getContext(),
            {
              'gff_file': {
                  'path': gff_path},
              'fasta_file': {
                  'path': fna_path},
              'workspace_name': self.getWsName(),
              'genome_name': ws_obj_name,
              'generate_missing_genes' : 1,
              'generate_ids_if_needed': 1
            })[0]
        data_file_cli = DataFileUtil(os.environ['SDK_CALLBACK_URL'])
        genome = data_file_cli.get_objects({'object_refs': [result['genome_ref']]})['data'][0]['data']
        found_spoofed_gene = False
        found_spoofed_gene_warning = False
        found_cds = False
        found_gene_cds = False
        found_cds_gene = False
        found_genome_warning = False
        suspect_genome = False

        has_genome_tiers = False
        has_user = False
        has_external_db = False
        if "genome_tier" in genome:
            has_genome_tiers = True
            for tier in genome["genome_tier"]:           
                if tier == "user":
                    has_user = True
        self.assertTrue(genome.get("source") == "User", "Source is not User : " + str(genome.get("source")))
        self.assertTrue(has_genome_tiers, "Does not have Genome Tiers")
        self.assertTrue(has_user, "Does not have User Genome Tier")   

        if "features" in genome:
            for feature in genome["features"]:
                print "Spoof GFF Feature ID : " + str(feature["id"])
                if feature['id'] == "Ga0123678_11.3":
                    found_spoofed_gene = True
                    if warnings['spoofed_gene'] in feature.get("warnings", []):
                        found_spoofed_gene_warning = True
                    if "cdss" in feature:
                        if feature["cdss"][0] == "Ga0123678_11.3_CDS_1":
                            found_gene_cds = True
                    temp_location = feature['location'][0]
                    self.assertEqual(temp_location[1], 1020)
                    self.assertEqual(temp_location[3], 627)  # first chunk
        if "cdss" in genome:
            for feature in genome["cdss"]:
                if feature['id'] == "Ga0123678_11.3_CDS_1":
                    found_cds = True
                    if feature["parent_gene"] == "Ga0123678_11.3":
                        found_cds_gene = True
        if "suspect" in genome:
            if int(genome["suspect"]) == 1:
                suspect_genome = True
        if warnings['spoofed_genome'].format(1) in genome.get("warnings", []):
            found_genome_warning = True
        self.assertTrue(found_cds,"The CDS was not found.")
        self.assertTrue(found_spoofed_gene,"The gene did not get spoofed.")
        self.assertTrue(found_spoofed_gene_warning,"The gene warning was not present.")
        self.assertTrue(found_cds_gene,"The Gene for the CDS was not found.")                                    
        self.assertTrue(found_gene_cds,"The CDS for the Gene was not found.") 
        self.assertTrue(found_genome_warning,"The genome warning was not present.")  
        self.assertTrue(suspect_genome,"The genome was not labeled as being suspect.")  

    def test_spoof_off(self):
        gff_path = "data/fasta_gff/JGI/Bacterial_Data/91705.assembled.gff"
        fna_path = "data/fasta_gff/JGI/Bacterial_Data/91705.assembled.fna"
        data_file_cli = DataFileUtil(os.environ['SDK_CALLBACK_URL'])
        ws_obj_name = 'JGI_Spoof_FAIL'
        with self.assertRaisesRegexp(
                        ValueError, warnings['no_spoof']):
            self.getImpl().fasta_gff_to_genome(
                self.getContext(),
                {
                'gff_file': {
                    'path': gff_path},
                'fasta_file': {
                    'path': fna_path},
                'workspace_name': self.getWsName(),
                'genome_name': ws_obj_name,
                'generate_ids_if_needed': 1
            })


