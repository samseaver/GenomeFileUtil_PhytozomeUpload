# -*- coding: utf-8 -*-
import unittest
import os  # noqa: F401
import json  # noqa: F401
import time
import shutil
import re

from os import environ
try:
    from ConfigParser import ConfigParser  # py2
except:
    from configparser import ConfigParser  # py3

from pprint import pprint  # noqa: F401

from biokbase.workspace.client import Workspace as workspaceService
from GenomeFileUtil.GenomeFileUtilImpl import GenomeFileUtil
from GenomeFileUtil.GenomeFileUtilImpl import SDKConfig
from GenomeFileUtil.GenomeFileUtilServer import MethodContext
from GenomeFileUtil.authclient import KBaseAuth as _KBaseAuth
from GenomeFileUtil.core.FastaGFFToGenome import FastaGFFToGenome
from DataFileUtil.DataFileUtilClient import DataFileUtil


class FastaGFFToGenomeUploadTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        print('setting up class')
        cls.token = environ.get('KB_AUTH_TOKEN', None)
        config_file = environ.get('KB_DEPLOYMENT_CONFIG', None)
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
        shutil.copy(os.path.join("data", "fasta_gff", "Plant_Data", cls.gff_filename), cls.gff_path)

        cls.fa_filename = 'Test_v1.0.fa.gz'
        cls.fa_path = os.path.join(cls.scratch, cls.fa_filename)
        shutil.copy(os.path.join("data", "fasta_gff", "Plant_Data", cls.fa_filename), cls.fa_path)

        cls.fungal_gff_filename = 'Neucr2.filtered_proteins.BroadModels.gff3.gz'
        cls.fungal_gff_path = os.path.join(cls.scratch, cls.fungal_gff_filename)
        shutil.copy(os.path.join("data", "fasta_gff", "Fungal_Data", cls.fungal_gff_filename),
                    cls.fungal_gff_path)

        cls.fungal_fa_filename = 'Neucr2_AssemblyScaffolds.fasta.gz'
        cls.fungal_fa_path = os.path.join(cls.scratch, cls.fungal_fa_filename)
        shutil.copy(os.path.join("data", "fasta_gff", "Fungal_Data", cls.fungal_fa_filename),
                    cls.fungal_fa_path)

    def check_minimal_items_exist(self, result):

        self.assertTrue('genome_info' in result)
        self.assertTrue('genome_ref' in result)
        self.assertTrue('report_name' in result)
        self.assertTrue('report_ref' in result)

        genome_info = result['genome_info']
        self.assertEquals(genome_info[10]['Domain'], 'Eukaryota')
        self.assertEquals(genome_info[10]['Genetic code'], '1')
        self.assertEquals(genome_info[10]['Name'], 'unknown_taxon')
        self.assertEquals(genome_info[10]['Source'], 'Genbank')
        self.assertTrue('GC content' in genome_info[10])
        self.assertTrue(re.match("^\d+?\.\d+?$", genome_info[10]['GC content']) is not None)
        self.assertTrue('Number features' in genome_info[10])
        self.assertTrue(genome_info[10]['Number features'].isdigit())
        self.assertTrue('Size' in genome_info[10])
        self.assertTrue(genome_info[10]['Size'].isdigit())
        self.assertEquals(genome_info[10]['Taxonomy'], 'Unknown')

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
        self.assertTrue('report_name' in result)
        self.assertTrue('report_ref' in result)

        genome_info = result['genome_info']
        self.assertEquals(genome_info[10]['Domain'], 'Eukaryota')
        self.assertEquals(genome_info[10]['Genetic code'], '1')
        self.assertEquals(genome_info[10]['Name'], 'unknown_taxon')
        self.assertEquals(genome_info[10]['Source'], 'User')
        self.assertTrue('GC content' in genome_info[10])
        self.assertTrue(re.match("^\d+?\.\d+?$", genome_info[10]['GC content']) is not None)
        self.assertTrue('Number features' in genome_info[10])
        self.assertTrue(genome_info[10]['Number features'].isdigit())
        self.assertTrue('Size' in genome_info[10])
        self.assertTrue(genome_info[10]['Size'].isdigit())
        self.assertEquals(genome_info[10]['Taxonomy'], 'Unknown')

    def test_simple_fasta_gff_to_genome(self):
        input_params = {
            'fasta_file': {'path': self.fa_path},
            'gff_file': {'path': self.gff_path},
            'genome_name': 'MyGenome',
            'workspace_name': self.getWsName(),
            'source': 'Genbank',
            'type': 'Reference',
            'scientific_name': 'Populus trichocarpa'
        }

        result = self.getImpl().fasta_gff_to_genome(self.getContext(), input_params)[0]

        self.assertTrue('genome_info' in result)
        self.assertTrue('genome_ref' in result)
        self.assertTrue('report_name' in result)
        self.assertTrue('report_ref' in result)

        genome_info = result['genome_info']
        self.assertEquals(genome_info[10]['Number features'], '1028')
        self.assertEquals(genome_info[10]['Domain'], 'Eukaryota')
        self.assertEquals(genome_info[10]['Genetic code'], '1')
        self.assertEquals(genome_info[10]['Name'], 'Populus trichocarpa')
        self.assertEquals(genome_info[10]['Source'], 'Genbank')
        self.assertTrue('GC content' in genome_info[10])
        self.assertTrue(re.match("^\d+?\.\d+?$", genome_info[10]['GC content']) is not None)
        self.assertTrue('Number features' in genome_info[10])
        self.assertTrue(genome_info[10]['Number features'].isdigit())
        self.assertTrue('Size' in genome_info[10])
        self.assertTrue(genome_info[10]['Size'].isdigit())
        self.assertEquals(genome_info[10]['Taxonomy'],
                          'cellular organisms; Eukaryota; Viridiplantae; Streptophyta; ' +
                          'Streptophytina; Embryophyta; Tracheophyta; Euphyllophyta; ' +
                          'Spermatophyta; Magnoliophyta; Mesangiospermae; eudicotyledons; ' +
                          'Gunneridae; Pentapetalae; rosids; fabids; Malpighiales; Salicaceae; ' +
                          'Saliceae; Populus')

    def test_taxon_reference_fasta_gff_to_genome(self):
        taxon_wsname = 'ReferenceTaxons'
        taxon_object_name = "unknown_taxon"
        taxon_info = self.dfu.get_objects({'object_refs': [taxon_wsname+"/"+taxon_object_name],
                                           'ignore_errors': 0})['data'][0]
        taxon_reference = "{}/{}/{}".format(taxon_info['info'][6],
                                            taxon_info['info'][0],
                                            taxon_info['info'][4])

        input_params = {
            'fasta_file': {'path': self.fa_path},
            'gff_file': {'path': self.gff_path},
            'genome_name': 'MyGenome',
            'workspace_name': self.getWsName(),
            'source': 'Genbank',
            'taxon_reference': taxon_reference,
            'type': 'Reference'
        }

        result = self.getImpl().fasta_gff_to_genome(self.getContext(), input_params)[0]

        self.check_minimal_items_exist(result)

    def test_shock_fasta_gff_to_genome(self):
        gff_shock_id = self.dfu.file_to_shock({'file_path': self.gff_path})['shock_id']
        fa_shock_id = self.dfu.file_to_shock({'file_path': self.fa_path})['shock_id']

        input_params = {
            'workspace_name': self.getWsName(),
            'genome_name': 'MyGenome',
            'fasta_file': {'shock_id': fa_shock_id},
            'gff_file': {'shock_id': gff_shock_id},
            'source': 'Genbank',
            'type': 'Reference'
        }

        result = self.getImpl().fasta_gff_to_genome(self.getContext(), input_params)[0]

        self.check_minimal_items_exist(result)

    def test_fungal_fasta_gff_to_genome(self):
        input_params = {
            'workspace_name': self.getWsName(),
            'genome_name': 'MyGenome',
            'fasta_file': {'path': self.fungal_fa_path},
            'gff_file': {'path': self.fungal_gff_path},
            'source': 'Genbank',
            'type': 'Reference'
        }

        result = self.getImpl().fasta_gff_to_genome(self.getContext(), input_params)[0]

        self.check_minimal_items_exist(result)

    def test_bad_fasta_gff_to_genome_params(self):
        invalidate_input_params = {
          'missing_workspace_name': 'workspace_name',
          'genome_name': 'genome_name',
          'fasta_file': {'path': 'fasta_file'},
          'gff_file': {'path': 'gff_file'}
        }
        with self.assertRaisesRegexp(
                    ValueError,
                    '"workspace_name" parameter is required, but missing'):
            self.getImpl().fasta_gff_to_genome(self.getContext(), invalidate_input_params)

        invalidate_input_params = {
          'workspace_name': 'workspace_name',
          'missing_genome_name': 'genome_name',
          'fasta_file': {'path': 'fasta_file'},
          'gff_file': {'path': 'gff_file'}
        }
        with self.assertRaisesRegexp(
                    ValueError,
                    '"genome_name" parameter is required, but missing'):
            self.getImpl().fasta_gff_to_genome(self.getContext(), invalidate_input_params)

        invalidate_input_params = {
          'workspace_name': 'workspace_name',
          'genome_name': 'genome_name',
          'missing_fasta_file': {'path': 'fasta_file'},
          'gff_file': {'path': 'gff_file'}
        }
        with self.assertRaisesRegexp(
                ValueError,
                '"fasta_file" parameter is required, but missing'):
            self.getImpl().fasta_gff_to_genome(self.getContext(), invalidate_input_params)

        invalidate_input_params = {
          'workspace_name': 'workspace_name',
          'genome_name': 'genome_name',
          'fasta_file': {'path': 'fasta_file'},
          'missing_gff_file': {'path': 'gff_file'}
        }
        with self.assertRaisesRegexp(
                ValueError,
                '"gff_file" parameter is required, but missing'):
            self.getImpl().fasta_gff_to_genome(self.getContext(), invalidate_input_params)

        invalidate_input_params = {
          'workspace_name': 'workspace_name',
          'genome_name': 'genome_name',
          'fasta_file': 'not_a_dict',
          'gff_file': {'path': 'gff_file'}
        }
        with self.assertRaisesRegexp(
                ValueError,
                'Required "fasta_file" field must be a map/dict'):
            self.getImpl().fasta_gff_to_genome(self.getContext(), invalidate_input_params)

        invalidate_input_params = {
          'workspace_name': 'workspace_name',
          'genome_name': 'genome_name',
          'fasta_file': {'path': 'fasta_file'},
          'gff_file': 'not_a_dict'
        }
        with self.assertRaisesRegexp(
                ValueError,
                'Required "gff_file" field must be a map/dict'):
            self.getImpl().fasta_gff_to_genome(self.getContext(), invalidate_input_params)

        invalidate_input_params = {
          'workspace_name': 'workspace_name',
          'genome_name': 'genome_name',
          'fasta_file': {'path': 'fasta_file'},
          'gff_file': {'ftp_url': 'gff_file'}
        }
        with self.assertRaisesRegexp(
                ValueError,
                'FTP link is currently not supported for FastaGFFToGenome'):
            self.getImpl().fasta_gff_to_genome(self.getContext(), invalidate_input_params)

        invalidate_input_params = {
          'workspace_name': 'workspace_name',
          'genome_name': 'genome_name',
          'fasta_file': {'missing_path': 'fasta_file'},
          'gff_file': {'path': 'gff_file'}
        }
        with self.assertRaisesRegexp(
                ValueError,
                'Required "fasta_file" field must include one source: path | shock_id'):
            self.getImpl().fasta_gff_to_genome(self.getContext(), invalidate_input_params)

        invalidate_input_params = {
          'workspace_name': 'workspace_name',
          'genome_name': 'genome_name',
          'fasta_file': {'path': 'fasta_file', 'shock_id': 'shock_id'},
          'gff_file': {'path': 'gff_file'}
        }
        with self.assertRaisesRegexp(
                ValueError,
                'Required "fasta_file" field has too many sources specified: '):
            self.getImpl().fasta_gff_to_genome(self.getContext(), invalidate_input_params)

        invalidate_input_params = {
          'workspace_name': 'workspace_name',
          'genome_name': 'genome_name',
          'fasta_file': {'path': 'fasta_file'},
          'gff_file': {'path': 'gff_file'},
          'type': 'invalid_type'
        }
        error_msg = 'Entered value for type is not one of the valid entries of '
        error_msg += '\["Reference", "User upload", "Representative"\]'
        with self.assertRaisesRegexp(
                ValueError, error_msg):
            self.getImpl().fasta_gff_to_genome(self.getContext(), invalidate_input_params)

    def test_FastaGFFToGenome_stage_input(self):
        # test absolute file path
        input_params = {
            'workspace_name': self.getWsName(),
            'genome_name': 'MyGenome',
            'fasta_file': {'path': self.fa_path},
            'gff_file': {'path': self.gff_path},
        }

        input_directory = os.path.join(self.scratch, 'test_FastaGFFToGenome_stage_input')
        os.makedirs(input_directory)

        file_paths = self.importer._stage_input(input_params, input_directory)

        self.assertTrue(self.gff_filename.rpartition('.')[0] in os.listdir(input_directory))
        self.assertTrue(self.fa_filename.rpartition('.')[0] in os.listdir(input_directory))
        self.assertTrue('gff_file' in file_paths)
        self.assertTrue('fasta_file' in file_paths)
        self.assertEquals(file_paths.get('gff_file'),
                          os.path.join(input_directory, self.gff_filename).rpartition('.')[0])
        self.assertEquals(file_paths.get('fasta_file'),
                          os.path.join(input_directory, self.fa_filename).rpartition('.')[0])

        # test shock id
        gff_shock_id = self.dfu.file_to_shock({'file_path': self.gff_path})['shock_id']
        fa_shock_id = self.dfu.file_to_shock({'file_path': self.fa_path})['shock_id']

        input_params = {
            'workspace_name': self.getWsName(),
            'genome_name': 'MyGenome',
            'fasta_file': {'shock_id': fa_shock_id},
            'gff_file': {'shock_id': gff_shock_id},
        }
        shutil.rmtree(input_directory)
        os.makedirs(input_directory)

        file_paths = self.importer._stage_input(input_params, input_directory)

        self.assertTrue(self.gff_filename.rpartition('.')[0] in os.listdir(input_directory))
        self.assertTrue(self.fa_filename.rpartition('.')[0] in os.listdir(input_directory))
        self.assertTrue('gff_file' in file_paths)
        self.assertTrue('fasta_file' in file_paths)
        self.assertEquals(file_paths.get('gff_file'),
                          os.path.join(input_directory, self.gff_filename).rpartition('.')[0])
        self.assertEquals(file_paths.get('fasta_file'),
                          os.path.join(input_directory, self.fa_filename).rpartition('.')[0])

    def test_FastaGFFToGenome_set_parsed_params(self):
        input_params = {
            'workspace_name': self.getWsName(),
            'genome_name': 'MyGenome',
            'fasta_file': {'path': self.fa_path},
            'gff_file': {'path': self.gff_path},
            'source': 'Genbank',
            'type': 'Reference'
        }

        parsed_params = self.importer._set_parsed_params(input_params)

        expect_param_keys = ['source', 'taxon_wsname', 'taxon_reference', 'release',
                             'type', 'metadata', 'workspace_name',
                             'genome_name', 'scientific_name', 'gff_file', 'fasta_file']
        self.assertItemsEqual(parsed_params.keys(), expect_param_keys)
        self.assertEquals(parsed_params['genome_name'], 'MyGenome')
        self.assertEquals(parsed_params['source'], 'Genbank')
        self.assertEquals(parsed_params['type'], 'Reference')

    def test_FastaGFFToGenome_retrieve_taxon(self):
        # test default none taxon_reference and none scientific_name
        taxon_reference = None
        taxon_wsname = 'ReferenceTaxons'
        scientific_name = 'unknown_taxon'

        expect_taxonomy, expect_taxon_reference = self.importer._retrieve_taxon(taxon_reference,
                                                                                taxon_wsname,
                                                                                scientific_name)
        self.assertEquals(expect_taxonomy, 'Unknown')
        self.assertTrue(expect_taxon_reference)

        # test given scientific_name
        scientific_name = "Populus trichocarpa"

        expect_taxonomy, expect_taxon_reference = self.importer._retrieve_taxon(taxon_reference,
                                                                                taxon_wsname,
                                                                                scientific_name)
        self.assertEquals(expect_taxonomy,
                          'cellular organisms; Eukaryota; Viridiplantae; Streptophyta; ' +
                          'Streptophytina; Embryophyta; Tracheophyta; Euphyllophyta; ' +
                          'Spermatophyta; Magnoliophyta; Mesangiospermae; eudicotyledons; ' +
                          'Gunneridae; Pentapetalae; rosids; fabids; Malpighiales; Salicaceae; ' +
                          'Saliceae; Populus')
        self.assertTrue(expect_taxon_reference)

        # test taxon_reference is not None
        taxon_object_name = "unknown_taxon"
        taxon_info = self.dfu.get_objects({'object_refs': [taxon_wsname+"/"+taxon_object_name],
                                           'ignore_errors': 0})['data'][0]
        taxon_reference = "{}/{}/{}".format(taxon_info['info'][6],
                                            taxon_info['info'][0],
                                            taxon_info['info'][4])

        expect_taxonomy, expect_taxon_reference = self.importer._retrieve_taxon(taxon_reference,
                                                                                taxon_wsname,
                                                                                scientific_name)
        self.assertEquals(expect_taxonomy, 'Unknown')
        self.assertEquals(expect_taxon_reference, taxon_reference)
        self.assertTrue(expect_taxon_reference)

    def test_FastaGFFToGenome_retrieve_fasta_file(self):

        input_fasta_file = self.dfu.unpack_file({'file_path': self.fa_path})['file_path']
        core_genome_name = 'MyGenome'
        scientific_name = 'unknown_taxon'
        source = 'Genbank'

        assembly = self.importer._retrieve_fasta_file(input_fasta_file, core_genome_name,
                                                      scientific_name, source)

        expect_assembly_keys = ['contigs', 'external_source', 'name',
                                'external_source_origination_date', 'assembly_id', 'notes',
                                'external_source_id', 'num_contigs', 'base_counts', 'gc_content',
                                'type', 'dna_size', 'md5']
        except_base_counts = {'A': 3697359, 'C': 1865502, 'T': 3674284, 'G': 1860999, 'N': 101776}
        expect_md5 = '638475f2ccf3cc8f3cd14e4881f0a979'
        expect_gc_content = 0.33

        self.assertItemsEqual(assembly.keys(), expect_assembly_keys)
        self.assertEquals(assembly['external_source'], source)
        self.assertEquals(assembly['name'], scientific_name)
        self.assertEquals(assembly['external_source_origination_date'],
                          str(os.stat(input_fasta_file).st_ctime))
        self.assertEquals(assembly['assembly_id'], core_genome_name+"_assembly")
        self.assertEquals(assembly['notes'],
                          'Note MD5s are generated from uppercasing the sequences')
        self.assertEquals(assembly['external_source_id'], os.path.basename(input_fasta_file))
        self.assertIsInstance(assembly['num_contigs'], int)
        self.assertIsInstance(assembly['base_counts'], dict)
        self.assertDictContainsSubset(assembly['base_counts'], except_base_counts)
        self.assertDictContainsSubset(except_base_counts, assembly['base_counts'])
        self.assertIsInstance(assembly['gc_content'], float)
        self.assertEquals(assembly['gc_content'], expect_gc_content)
        self.assertEquals(assembly['type'], 'Unknown')
        self.assertIsInstance(assembly['dna_size'], int)
        self.assertEquals(assembly['md5'], expect_md5)
        self.assertIsInstance(assembly['contigs'], dict)

    def test_FastaGFFToGenome_retrieve_gff_file(self):

        input_gff_file = self.dfu.unpack_file({'file_path': self.gff_path})['file_path']

        feature_list = self.importer._retrieve_gff_file(input_gff_file)

        self.assertIsInstance(feature_list, dict)

    def test_FastaGFFToGenome_retrieve_feature_identifiers(self):
        feature_list = {'Chr01':
                        [
                            {'end': 8201443,
                             'Name': 'Potri.001G102800',
                             'start': 8200895,
                             'score': '.',
                             'phase': '.',
                             'contig': 'Chr01',
                             'type': 'gene',
                             'ID': 'Potri.001G102800.v3.0',
                             'strand': '-'},
                            {'end': 8201443,
                             'Name': 'Potri.001G102800.1',
                             'Parent': 'Potri.001G102800.v3.0',
                             'pacid': '27047128',
                             'start': 8200895,
                             'score': '.',
                             'longest': '1',
                             'phase': '.',
                             'contig': 'Chr01',
                             'type': 'mRNA',
                             'ID': 'Potri.001G102800.1.v3.0',
                             'strand': '-'},
                            {'end': 8201443,
                             'Parent': 'Potri.001G102800.1.v3.0',
                             'pacid': '27047128',
                             'start': 8200895,
                             'score': '.',
                             'phase': '0',
                             'contig': 'Chr01',
                             'type': 'CDS',
                             'ID': 'Potri.001G102800.1.v3.0.CDS.1',
                             'strand': '-'}
                        ]}

        (features_identifiers_dict,
         features_identifiers_list,
         features_identifiers_count) = self.importer._retrieve_feature_identifiers(feature_list)

        expect_features_identifiers_dict = {'Potri.001G102800':
                                            {'Potri.001G102800.1':
                                             {'Potri.001G102800.1.CDS.1': 1}}}

        expect_features_identifiers_count = {'Potri.001G102800.1.CDS.1': 2,
                                             'Potri.001G102800.1': 1,
                                             'Potri.001G102800': 0, }

        self.assertEquals(features_identifiers_dict, expect_features_identifiers_dict)
        self.assertEquals(features_identifiers_count, expect_features_identifiers_count)
        self.assertEquals(features_identifiers_list, feature_list['Chr01'])

    def test_FastaGFFToGenome_update_feature_identifiers(self):
        features_identifiers_dict = {'Potri.001G102800':
                                     {'Potri.001G102800.1':
                                      {'Potri.001G102800.1.CDS.1': 1}}}
        features_identifiers_count = {'Potri.001G102800.1.CDS.1': 2,
                                      'Potri.001G102800.1': 1,
                                      'Potri.001G102800': 0, }
        features_identifiers_list = [{'end': 8201443,
                                      'Name': 'Potri.001G102800',
                                      'start': 8200895,
                                      'score': '.',
                                      'phase': '.',
                                      'contig': 'Chr01',
                                      'type': 'gene',
                                      'ID': 'Potri.001G102800.v3.0',
                                      'strand': '-'},
                                     {'end': 8201443,
                                      'Name': 'Potri.001G102800.1',
                                      'Parent': 'Potri.001G102800.v3.0',
                                      'pacid': '27047128',
                                      'start': 8200895,
                                      'score': '.',
                                      'longest': '1',
                                      'phase': '.',
                                      'contig': 'Chr01',
                                      'type': 'mRNA',
                                      'ID': 'Potri.001G102800.1.v3.0',
                                      'strand': '-'},
                                     {'end': 8201443,
                                      'Parent': 'Potri.001G102800.1.v3.0',
                                      'pacid': '27047128',
                                      'start': 8200895,
                                      'score': '.',
                                      'phase': '0',
                                      'contig': 'Chr01',
                                      'type': 'CDS',
                                      'ID': 'Potri.001G102800.1.v3.0.CDS.1',
                                      'strand': '-'}]

        (updated_features_identifiers_dict,
         updated_features_list,
         updated_features_identifiers_count) = self.importer._update_feature_identifiers(
                                                                features_identifiers_dict,
                                                                features_identifiers_list,
                                                                features_identifiers_count)

        expect_updated_features_identifiers_dict = {'Potri.001G102800.v3.0':
                                                    {'Potri.001G102800.1.v3.0':
                                                     {'Potri.001G102800.1.v3.0.CDS.1': 1}}}

        expect_updated_features_identifiers_count = {'Potri.001G102800.1.v3.0.CDS.1': 2,
                                                     'Potri.001G102800.1.v3.0': 1,
                                                     'Potri.001G102800.v3.0': 0}

        self.assertEquals(updated_features_identifiers_dict,
                          expect_updated_features_identifiers_dict)
        self.assertEquals(updated_features_identifiers_count,
                          expect_updated_features_identifiers_count)
        self.assertEquals(updated_features_list, features_identifiers_list)
