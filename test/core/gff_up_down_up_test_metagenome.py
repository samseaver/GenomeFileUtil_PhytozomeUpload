import os
import time
import json
import unittest
from configparser import ConfigParser

from installed_clients.DataFileUtilClient import DataFileUtil
from GenomeFileUtil.GenomeFileUtilImpl import GenomeFileUtil
from GenomeFileUtil.GenomeFileUtilServer import MethodContext
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
        # get metagenome data.
        cls.gff_path = "data/metagenomes/ebi/59111.assembled.gff"
        cls.fasta_path = "data/metagenomes/ebi/59111.assembled.fna"
        if not os.path.isfile(cls.fasta_path) or not os.path.isfile(cls.gff_path):
            raise InputError(f'Files {cls.gff_path} and/or {cls.fasta_path} not in test directory ')

        suffix = int(time.time() * 1000)
        cls.wsName = "test_GenomeFileUtil_" + str(suffix)
        cls.wsClient.create_workspace({'workspace': cls.wsName})

        print('Uploading GFF file')
        result = cls.serviceImpl.fasta_gff_to_metagenome(cls.ctx, {
            'workspace_name': cls.wsName,
            'genome_name': 'MyGenome',
            'fasta_file': {'path': cls.fasta_path},
            'gff_file': {'path': cls.gff_path},
            'source': 'GFF',
            'taxon_id': '3702',
            'type': 'Reference',
            'genome_type': 'Metagenome',
            'is_metagenome': True,
            'generate_missing_genes': True
        })[0]
        data_file_cli = DataFileUtil(os.environ['SDK_CALLBACK_URL'])
        cls.metagenome_ref = result['metagenome_ref']
        cls.genome_orig = data_file_cli.get_objects(
            {'object_refs': [result['metagenome_ref']]})['data'][0]['data']

        print('testing GFF download by building the file')
        down_result = cls.serviceImpl.metagenome_to_gff(
            cls.ctx, {'metagenome_ref': result['metagenome_ref']})[0]

        print('Reuploading GFF file')
        new_result = cls.serviceImpl.fasta_gff_to_metagenome(cls.ctx, {
            'workspace_name': cls.wsName,
            'genome_name': 'MyGenome',
            'fasta_file': {'path': cls.fasta_path},
            'gff_file': {'path': down_result['file_path']},
            'source': 'GFF',
            'type': 'Reference',
            'genome_type': 'Metagenome',
            'is_metagenome': True,
            'generate_missing_genes': True,
            'taxon_id': '3702',
        })[0]
        cls.genome_new = data_file_cli.get_objects({'object_refs': [new_result['metagenome_ref']]})['data'][0]['data']

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def test_gff_and_metagenome_to_metagenome(self):
        dfu = DataFileUtil(os.environ['SDK_CALLBACK_URL'])
        result = self.serviceImpl.ws_obj_gff_to_metagenome(self.ctx, {
            'workspace_name': self.wsName,
            'genome_name': 'MyGenome',
            'gff_file': {'path': self.gff_path},
            'ws_ref': self.metagenome_ref,
            'source': 'GFF',
            'type': 'Reference',
            'genome_type': 'Metagenome',
            'is_metagenome': True,
            'generate_missing_genes': True,
            'taxon_id': '3702',
        })[0]
        metagenome = dfu.get_objects({'object_refs': [result['metagenome_ref']]})['data'][0]['data']
        # make sure its same as original
        self._compare_features(self.genome_orig, metagenome)

    def _compare_features(self, metagenome_orig, metagenome_new):
        scratch_dir = self.cfg['scratch']

        dfu = DataFileUtil(os.environ['SDK_CALLBACK_URL'])
        orig_file_name = dfu.shock_to_file({'file_path': scratch_dir,
                                            'handle_id': metagenome_orig['features_handle_ref'],
                                            'unpack': 'unpack'
                                            })['file_path']
        new_file_name = dfu.shock_to_file({'file_path': scratch_dir,
                                           'handle_id': metagenome_new['features_handle_ref'],
                                           'unpack': 'unpack'
                                           })['file_path']

        # open json files
        with open(orig_file_name) as fid:
            metagenome_orig_data = json.load(fid)
        with open(new_file_name) as fid:
            metagenome_new_data = json.load(fid)

        print('Testing length or original vs new genome')
        self.assertTrue(len(metagenome_orig_data) == len(metagenome_new_data),
                        "list is not of equal length in Original and New Genomes.")
        print("\n\n" + " TOTAL NUMBER:" + str(len(metagenome_orig_data)))

        orig_dict = dict([(x['id'], x) for x in metagenome_orig_data])
        new_dict = dict([(x['id'], x) for x in metagenome_new_data])

        first_pass_matches = 0
        first_pass_non_match = 0
        second_pass_matches = 0

        print('Testing keys in metagenomes....')
        for key in orig_dict:
            orig_feature = orig_dict[key]
            new_feature = new_dict[key]
            if "aliases" in orig_feature:
                orig_feature['aliases'] = sorted(orig_feature.get('aliases', []))
                new_feature['aliases'] = sorted(new_feature.get('aliases', []))
            if "db_xrefs" in orig_feature:
                orig_feature['db_xrefs'] = sorted(orig_feature.get('db_xrefs', []))
                new_feature['db_xrefs'] = sorted(new_feature.get('db_xrefs', []))
            if "functions" in orig_feature:
                orig_feature["functions"] = sorted(orig_feature.get('functions', []))
                new_feature["functions"] = sorted(new_feature.get('functions', []))
            if orig_feature == new_feature:
                first_pass_matches += 1
            else:
                first_pass_non_match += 1
                orig_feature.pop("note", None)
                new_feature.pop("note", None)
                orig_feature.pop('inference_data', None)
                new_feature.pop('inference_data', None)
                if "warnings" in orig_feature and "warnings" not in new_feature:
                    del(orig_feature["warnings"])
                if orig_feature == new_feature:
                    second_pass_matches += 1
                else:
                    self.maxDiff = None
                    self.assertEqual(orig_feature, new_feature)
        self.assertEqual(
            len(orig_dict),
            (first_pass_matches + second_pass_matches),
            (f"There were {first_pass_matches} first pass matches "
             f"and {second_pass_matches} second pass matches out of "
             f"{len(orig_dict)} items in features")
        )

    def test_feature_list_comparison(self):
        metagenome_orig = self.genome_orig
        metagenome_new = self.genome_new
        self._compare_features(metagenome_orig, metagenome_new)

    def test_metagenome_gff_export(self):
        # fetch the test files and set things up
        res = self.serviceImpl.export_metagenome_as_gff(
            self.ctx, {'input_ref': self.metagenome_ref})[0]
        self.assertTrue('shock_id' in res)
