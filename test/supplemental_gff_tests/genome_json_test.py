# -*- coding: utf-8 -*-
import os  # noqa: F401
import re
import shutil
import time
import unittest
from configparser import ConfigParser
from os import environ

from installed_clients.DataFileUtilClient import DataFileUtil
from GenomeFileUtil.GenomeFileUtilImpl import GenomeFileUtil
from GenomeFileUtil.GenomeFileUtilImpl import SDKConfig
from GenomeFileUtil.GenomeFileUtilServer import MethodContext
from GenomeFileUtil.authclient import KBaseAuth as _KBaseAuth
from GenomeFileUtil.core.FastaGFFToGenome import FastaGFFToGenome
from installed_clients.WorkspaceClient import Workspace as workspaceService


class FastaGFFToGenomeUploadTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        print('setting up class')
        cls.token = environ.get('KB_AUTH_TOKEN')
        config_file = environ['KB_DEPLOYMENT_CONFIG']
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('GenomeFileUtil'):
            cls.cfg[nameval[0]] = nameval[1]
        authServiceUrl = cls.cfg.get('auth-service-url',
                                     "https://kbase.us/services/authorization/Sessions/Login")
        auth_client = _KBaseAuth(authServiceUrl)
        cls.user_id = auth_client.get_user(cls.token)
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': cls.token,
                        'user_id': cls.user_id,
                        'provenance': [
                            {'service': 'GenomeFileUtil',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = workspaceService(cls.wsURL, token=cls.token)
        cls.serviceImpl = GenomeFileUtil(cls.cfg)
        cls.dfu = DataFileUtil(os.environ['SDK_CALLBACK_URL'], token=cls.token)
        cls.scratch = cls.cfg['scratch']
        cls.shockURL = cls.cfg['shock-url']
        cls.gfu_cfg = SDKConfig(cls.cfg)
        cls.prepare_data()

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
        ret = self.getWsClient().create_workspace({'workspace': wsName})  # noqa
        self.__class__.wsName = wsName
        return wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    @classmethod
    def prepare_data(cls):
        cls.importer = FastaGFFToGenome(cls.gfu_cfg)

        cls.gff_filename = 'Test_v1.0.gene.gff3.gz'
        cls.gff_path = os.path.join(cls.scratch, cls.gff_filename)
        shutil.copy(os.path.join("data", "fasta_gff", "JGI", "Plant_Data", cls.gff_filename), cls.gff_path)

        cls.fa_filename = 'Test_v1.0.fa.gz'
        cls.fa_path = os.path.join(cls.scratch, cls.fa_filename)
        shutil.copy(os.path.join("data", "fasta_gff", "JGI", "Plant_Data", cls.fa_filename), cls.fa_path)

        cls.fungal_gff_filename = 'Neucr2.filtered_proteins.BroadModels.gff3.gz'
        cls.fungal_gff_path = os.path.join(cls.scratch, cls.fungal_gff_filename)
        shutil.copy(os.path.join("data", "fasta_gff", "JGI", "Fungal_Data", cls.fungal_gff_filename),
                    cls.fungal_gff_path)

        cls.fungal_fa_filename = 'Neucr2_AssemblyScaffolds.fasta.gz'
        cls.fungal_fa_path = os.path.join(cls.scratch, cls.fungal_fa_filename)
        shutil.copy(os.path.join("data", "fasta_gff", "JGI", "Fungal_Data", cls.fungal_fa_filename),
                    cls.fungal_fa_path)

        cls.jgi_bacterial_gff_filename = '2547132501.gff.gz'
        cls.jgi_bacterial_gff_path = os.path.join(cls.scratch, cls.jgi_bacterial_gff_filename)
        shutil.copy(os.path.join("data", "fasta_gff", "JGI", "Bacterial_Data", cls.jgi_bacterial_gff_filename),
                    cls.jgi_bacterial_gff_path)

        cls.jgi_bacterial_fa_filename = '2547132501.fna.gz'
        cls.jgi_bacterial_fa_path = os.path.join(cls.scratch, cls.jgi_bacterial_fa_filename)
        shutil.copy(os.path.join("data", "fasta_gff", "JGI", "Bacterial_Data", cls.jgi_bacterial_fa_filename),
                    cls.jgi_bacterial_fa_path)

        cls.jgi_bacterial_gff2_filename = '91705.assembled.gff'
        cls.jgi_bacterial_gff2_path = os.path.join(cls.scratch, cls.jgi_bacterial_gff2_filename)
        shutil.copy(os.path.join("data", "fasta_gff", "JGI", "Bacterial_Data", cls.jgi_bacterial_gff2_filename),
                    cls.jgi_bacterial_gff2_path)

        cls.jgi_bacterial_fa2_filename = '91705.assembled.fna'
        cls.jgi_bacterial_fa2_path = os.path.join(cls.scratch, cls.jgi_bacterial_fa2_filename)
        shutil.copy(os.path.join("data", "fasta_gff", "JGI", "Bacterial_Data", cls.jgi_bacterial_fa2_filename),
                    cls.jgi_bacterial_fa2_path)

        cls.patric_bacterial_gff_filename = '1240778.3.PATRIC.gff.gz'
        cls.patric_bacterial_gff_path = os.path.join(cls.scratch, cls.patric_bacterial_gff_filename)
        shutil.copy(os.path.join("data", "fasta_gff", "PATRIC", "Ecoli_O104", cls.patric_bacterial_gff_filename),
                    cls.patric_bacterial_gff_path)

        cls.patric_bacterial_fa_filename = '1240778.3.fna.gz'
        cls.patric_bacterial_fa_path = os.path.join(cls.scratch, cls.patric_bacterial_fa_filename)
        shutil.copy(os.path.join("data", "fasta_gff", "PATRIC", "Ecoli_O104", cls.patric_bacterial_fa_filename),
                    cls.patric_bacterial_fa_path)

        cls.refseq_bacterial_gff_filename = 'NC_021490.gff.gz'
        cls.refseq_bacterial_gff_path = os.path.join(cls.scratch, cls.refseq_bacterial_gff_filename)
        shutil.copy(os.path.join("data", "fasta_gff", "RefSeq", "Bacterial_Data", cls.refseq_bacterial_gff_filename),
                    cls.refseq_bacterial_gff_path)

        cls.refseq_bacterial_fa_filename = 'NC_021490.fasta.gz'
        cls.refseq_bacterial_fa_path = os.path.join(cls.scratch, cls.refseq_bacterial_fa_filename)
        shutil.copy(os.path.join("data", "fasta_gff", "RefSeq", "Bacterial_Data", cls.refseq_bacterial_fa_filename),
                    cls.refseq_bacterial_fa_path)

    def check_minimal_items_exist(self, result):
        self.assertTrue('genome_info' in result)
        self.assertTrue('genome_ref' in result)
        genome_info = result['genome_info']
        self.assertEqual(genome_info[10]['Domain'], 'Unknown')
        self.assertEqual(genome_info[10]['Genetic code'], '11')
        self.assertEqual(genome_info[10]['Name'], 'unknown_taxon')
        self.assertEqual(genome_info[10]['Source'], 'Genbank')
        self.assertTrue('GC content' in genome_info[10])
        self.assertTrue(re.match(r"^\d+?\.\d+?$", genome_info[10]['GC content']) is not None)
        self.assertTrue('Number of Protein Encoding Genes' in genome_info[10])
        self.assertTrue(genome_info[10]['Number of Protein Encoding Genes'].isdigit())
        self.assertTrue('Size' in genome_info[10])
        self.assertTrue(genome_info[10]['Size'].isdigit())
        self.assertEqual(genome_info[10]['Taxonomy'], 'Unconfirmed Organism: unknown_taxon')

    def print_genome_warnings(self, result):
        data_file_cli = DataFileUtil(os.environ['SDK_CALLBACK_URL'],
                                     service_ver='dev')
        genome = data_file_cli.get_objects({'object_refs': [result['genome_ref']]})['data'][0]['data']
        if 'warnings' in genome:
            print("Genome warnings:" + str(genome['warnings']))

    def check_CDS_warnings(self, result, test_name):
        data_file_cli = DataFileUtil(os.environ['SDK_CALLBACK_URL'],
                                     service_ver='dev')
        genome = data_file_cli.get_objects({'object_refs': [result['genome_ref']]})['data'][0]['data']
        print("IN TEST NAME : " + str(test_name))
        cds_warning_count = 0
        cds_with_warning_count = 0
        if 'cdss' in genome:
            total_cds_count = len(genome['cdss'])
            for feature in genome['cdss']:
                if 'warnings' in feature:
                    if test_name == "test_jgi_bacterial_fasta_gff2_to_genome":
                        print(str(feature['id']) + " warnings:" + str(feature['warnings']))
                        print("Location: " + str(feature['location']))
                        print("Translation: " + feature['protein_translation'])
                        print("DNA Sequence: " + feature["dna_sequence"])
                    cds_with_warning_count = cds_with_warning_count + 1
                    cds_warning_count = cds_warning_count + len(feature['warnings'])

            print("Total CDS: " + str(total_cds_count))
            print("CDS Warning Count: " + str(cds_warning_count))
            print("CDSs with a warning Count: " + str(cds_with_warning_count))
            print("Percent CDS with warning: " + str((cds_with_warning_count/float(total_cds_count)) * 100))

    def test_simple_fasta_gff_to_genome_w_null_params(self):
        input_params = {
            "fasta_file": {'path': self.fa_path},
            "gff_file": {'path': self.gff_path},
            "workspace_name": self.getWsName(),
            "genome_name": 'MyGenome',
            "scientific_name": None,
            "taxon_reference": None,
            "genetic_code": None,
            "source": None,
            "taxon_wsname": None,
            "release": None,
            "type": None
        }
        result = self.getImpl().fasta_gff_to_genome(self.getContext(), input_params)[0]
        self.assertTrue('genome_info' in result)
        self.assertTrue('genome_ref' in result)
        genome_info = result['genome_info']
        self.assertEqual(genome_info[10]['Domain'], 'Unknown')
        self.assertEqual(genome_info[10]['Genetic code'], '11')
        self.assertEqual(genome_info[10]['Name'], 'unknown_taxon')
        self.assertEqual(genome_info[10]['Source'], 'User')
        self.assertTrue('GC content' in genome_info[10])
        self.assertTrue(re.match(r"^\d+?\.\d+?$", genome_info[10]['GC content']) is not None)
        self.assertTrue('Number of Protein Encoding Genes' in genome_info[10])
        self.assertTrue(genome_info[10]['Number of Protein Encoding Genes'].isdigit())
        self.assertTrue('Size' in genome_info[10])
        self.assertTrue(genome_info[10]['Size'].isdigit())
        self.assertEqual(genome_info[10]['Taxonomy'], 'Unconfirmed Organism')

    def test_fasta_gff_to_genome_json(self):
        input_params = {
            'fasta_file': {'path': self.fa_path},
            'gff_file': {'path': self.gff_path},
            'genome_name': 'Plant',
            'workspace_name': self.getWsName(),
            'source': 'Genbank',
            'type': 'Reference',
            'taxon_id': 3694,
            'scientific_name': 'Populus trichocarpa'
        }
        genome_json = self.getImpl().fasta_gff_to_genome_json(self.getContext(), input_params)[0][0]
        assert 'features' in genome_json
        assert 'feature_counts' in genome_json
        assert 'genome_tiers' in genome_json
        self.assertEqual(genome_json['domain'], 'Eukaryota')
        self.assertEqual(genome_json['genetic_code'], 1)
        self.assertEqual(genome_json['scientific_name'], 'Populus trichocarpa')
        self.assertEqual(genome_json['source'], 'Genbank')
        self.assertTrue('gc_content' in genome_json)
        self.assertTrue('dna_size' in genome_json)
        self.assertEqual(genome_json['taxonomy'],
                         'cellular organisms; Eukaryota; Viridiplantae; Streptophyta; ' +
                         'Streptophytina; Embryophyta; Tracheophyta; Euphyllophyta; ' +
                         'Spermatophyta; Magnoliopsida; Mesangiospermae; eudicotyledons; ' +
                         'Gunneridae; Pentapetalae; rosids; fabids; Malpighiales; Salicaceae; ' +
                         'Saliceae; Populus')
