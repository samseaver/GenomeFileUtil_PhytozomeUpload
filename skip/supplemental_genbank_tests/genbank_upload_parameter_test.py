import unittest
import os
import time
import shutil
from Bio import SeqIO

from os import environ
try:
    from ConfigParser import ConfigParser  # py2
except:
    from configparser import ConfigParser  # py3


from biokbase.workspace.client import Workspace as workspaceService
from GenomeFileUtil.GenomeFileUtilImpl import GenomeFileUtil
from GenomeFileUtil.GenomeFileUtilServer import MethodContext


class MinimalGenbankUploadTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        print('setting up class')
        token = environ.get('KB_AUTH_TOKEN', None)
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
        config_file = environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('GenomeFileUtil'):
            cls.cfg[nameval[0]] = nameval[1]
        cls.wsURL = cls.cfg['workspace-url']

        cls.ws = workspaceService(cls.wsURL, token=token)
        cls.impl = GenomeFileUtil(cls.cfg)

        cls.MINIMAL_TEST_FILE = os.path.join( cls.cfg['scratch'], 'minimal.gbff')
        shutil.copy('data/minimal.gbff', cls.MINIMAL_TEST_FILE )

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.ws.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def getWsClient(self):
        return self.ws

    def getWsName(self):
        if hasattr(self.__class__, 'wsName'):
            return self.__class__.wsName
        suffix = int(time.time() * 1000)
        wsName = "test_GenomeFileUtil_" + str(suffix)
        ret = self.getWsClient().create_workspace({'workspace': wsName})
        self.__class__.wsName = wsName
        return wsName

    def getImpl(self):
        return self.__class__.impl

    def getContext(self):
        return self.__class__.ctx

    def test_upload(self):
        return
        # fetch the test files and set things up
        genomeFileUtil = self.getImpl()
        gbk_path = self.MINIMAL_TEST_FILE

        # ok, first test with minimal options
        result = genomeFileUtil.genbank_to_genome(self.getContext(),
                                    {
                                        'file':{'path': gbk_path},
                                        'workspace_name': self.getWsName(),
                                        'genome_name': 'something',
                                    })[0]
        self.check_minimal_items_exist(result)

        # test without using ontologies, and with setting a taxon_reference directly
        result = genomeFileUtil.genbank_to_genome(self.getContext(),
                                    {
                                        'file':{'path': gbk_path},
                                        'workspace_name': self.getWsName(),
                                        'genome_name': 'something',
                                        'exclude_ontologies':1,
                                        'taxon_reference':'ReferenceTaxons/4932_taxon'
                                    })[0]
        self.check_minimal_items_exist(result)

        # test setting additional metadata
        result = genomeFileUtil.genbank_to_genome(self.getContext(),
                                    {
                                        'file':{'path': gbk_path},
                                        'workspace_name': self.getWsName(),
                                        'genome_name': 'something',
                                        'exclude_ontologies':1,
                                        'taxon_reference':'ReferenceTaxons/4932_taxon',
                                        'metadata': { 'mydata' : 'yay', 'otherdata':'ok' }
                                    })[0]
        self.check_minimal_items_exist(result)
        metadata_saved = result['genome_info'][10]
        self.assertTrue('mydata' in metadata_saved)
        self.assertTrue('otherdata' in metadata_saved)
        self.assertEquals(metadata_saved['mydata'], 'yay')

    def check_minimal_items_exist(self, result):

        self.assertTrue('genome_info' in result)
        self.assertTrue('genome_ref' in result)
        self.assertTrue('report_name' in result)
        self.assertTrue('report_ref' in result)

        genome_info = result['genome_info']
        self.assertEquals(genome_info[10]['Number contigs'],'1')
        self.assertEquals(genome_info[10]['Number of Protein Encoding Genes'],'2')
        self.assertEquals(genome_info[10]['Domain'],'Eukaryota')
        self.assertEquals(genome_info[10]['Genetic code'],'1')
        self.assertEquals(genome_info[10]['Name'],'Saccharomyces cerevisiae')
        self.assertEquals(genome_info[10]['Source'], 'Genbank')
        self.assertEquals(genome_info[10]['GC content'], '0.37967')
        self.assertEquals(genome_info[10]['Size'], '5028')
        self.assertEquals(genome_info[10]['Taxonomy'],
            'cellular organisms; Eukaryota; Opisthokonta; Fungi; Dikarya; Ascomycota; '+
            'saccharomyceta; Saccharomycotina; Saccharomycetes; Saccharomycetales; '+
            'Saccharomycetaceae; Saccharomyces')

    def test_supply_assembly(self):
        genomeFileUtil = self.getImpl()
        """Warning: This test will fail if not run against CI"""
        gbk_path = self.MINIMAL_TEST_FILE
        with self.assertRaisesRegexp(ValueError, "not a valid format."):
            result = genomeFileUtil.genbank_to_genome(
                self.getContext(), {
                                      'file': {'path': gbk_path},
                                      'workspace_name': self.getWsName(),
                                      'genome_name': 'something',
                                      'use_existing_assembly': "1",
                                  })[0]
        with self.assertRaisesRegexp(ValueError, "not a reference to an assembly"):
            result = genomeFileUtil.genbank_to_genome(
                self.getContext(), {
                    'file': {'path': gbk_path},
                    'workspace_name': self.getWsName(),
                    'genome_name': 'something',
                    'use_existing_assembly': "6976/923/6",
                })[0]
        with self.assertRaisesRegexp(ValueError, "contains contigs which are not present"):
            result = genomeFileUtil.genbank_to_genome(
                self.getContext(), {
                    'file': {'path': gbk_path},
                    'workspace_name': self.getWsName(),
                    'genome_name': 'something',
                    'use_existing_assembly': "31767/5/1",
                })[0]

    def test_translation(self):
        import string
        record = next(SeqIO.parse(open(self.MINIMAL_TEST_FILE), 'genbank'))
        f_seq = str(record.seq)
        r_seq = f_seq.translate(string.maketrans("CTAG", "GATC"))

        def _location(feat):
            strand_trans = ("", "+", "-")
            loc = []
            for part in feat.location.parts:
                if part.strand >= 0:
                    begin = int(part.start) + 1
                else:
                    begin = int(part.end)
                loc.append((
                        record.id,
                        begin,
                        strand_trans[part.strand],
                        len(part)))
            return loc

        def get_seq(feat):
            seq = []
            strand = 1
            for part in feat.location.parts:
                strand = part.strand
                if strand >= 0:
                    seq.append(f_seq[part.start:part.end])
                else:
                    seq.insert(0, r_seq[part.start:part.end])
            if strand >= 0:
                return "".join(seq)
            else:
                return "".join(seq)[::-1]

        for feat in record.features:
            print feat.id
            seq1 = feat.extract(record)
            seq2 = get_seq(feat)
            self.assertEqual(str(seq1.seq), seq2)

