import os
import time
import json
import gzip
import shutil
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
        gff_path   = "data/metagenomes/ebi/59111.assembled.gff"
        fasta_path = "data/metagenomes/ebi/59111.assembled.fna"
        ws_obj_name = 'metagenome_test_objects'
        suffix = int(time.time() * 1000)
        cls.wsName = "test_GenomeFileUtil_" + str(suffix)
        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})

        print('Uploading GFF file')
        result = cls.serviceImpl.fasta_gff_to_genome(
            cls.ctx,
            {
                'workspace_name': cls.wsName,
                'genome_name': 'MyGenome',
                'fasta_file': {'path': fasta_path},
                'gff_file': {'path': gff_path},
                'source': 'GFF',
                'type': 'Reference',
                'genome_type': 'Metagenome',
                'is_metagenome': True,
                'generate_missing_genes': True
            })[0]
        data_file_cli = DataFileUtil(os.environ['SDK_CALLBACK_URL'])
        cls.genome_orig = data_file_cli.get_objects(
            {'object_refs': [result['genome_ref']]})['data'][0]['data']

        print('testing GFF download by building the file')
        down_result = cls.serviceImpl.metagenome_to_gff(
            cls.ctx, {'genome_ref': result['genome_ref']})[0]

        print('Reuploading GFF file')
        new_result = cls.serviceImpl.fasta_gff_to_genome(
            cls.ctx,
            {
                'workspace_name': cls.wsName,
                'genome_name': 'MyGenome',
                'fasta_file': {'path': fasta_path},
                'gff_file': {'path': down_result['file_path']},
                'source': 'GFF',
                'type': 'Reference',
                'genome_type': 'Metagenome',
                'is_metagenome': True,
                'generate_missing_genes': True
            })[0]
        cls.genome_new = data_file_cli.get_objects({'object_refs': [new_result['genome_ref']]})['data'][0]['data']

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def test_feature_list_comparison(self):
        metagenome_orig = self.genome_orig
        metagenome_new = self.genome_new

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
                note_orig = orig_feature.pop("note",None)
                note_new = new_feature.pop("note",None)
                inference_orig = orig_feature.pop('inference_data',None)
                inference_new = new_feature.pop('inference_data',None)
                if "warnings" in orig_feature and "warnings" not in new_feature:
                    del(orig_feature["warnings"])
                if orig_feature == new_feature:
                    second_pass_matches += 1
                else:
                    self.maxDiff = None
                    self.assertEqual(orig_feature,new_feature)
        self.assertEqual(len(orig_dict),(first_pass_matches + second_pass_matches),
                        "There were %d first pass matches and %d second pass matches out of %d items in features" %
                        (first_pass_matches, second_pass_matches, len(orig_dict)))
