import unittest
import os
import json
import time
import shutil
import filecmp

from os import environ
try:
    from ConfigParser import ConfigParser  # py2
except:
    from configparser import ConfigParser  # py3

from Workspace.WorkspaceClient import Workspace as workspaceService
from GenomeFileUtil.GenomeFileUtilImpl import GenomeFileUtil
from GenomeFileUtil.GenomeFileUtilServer import MethodContext
from AssemblyUtil.AssemblyUtilClient import AssemblyUtil
from GenomeFileUtil.core.GenomeInterface import GenomeInterface
from GenomeFileUtil.GenomeFileUtilImpl import SDKConfig


class GenomeFileUtilTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
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
        cls.serviceImpl = GenomeFileUtil(cls.cfg)
        gi_config = SDKConfig(cls.cfg)
        cls.genome_interface = GenomeInterface(gi_config)
        # create one WS for all tests
        suffix = int(time.time() * 1000)
        wsName = "test_GenomeAnnotationAPI_" + str(suffix)
        ret = cls.ws.create_workspace({'workspace': wsName})
        cls.wsName = wsName

        # save new genome
        assembly_file_path = os.path.join(cls.cfg['scratch'],
                                          'Rhodo_velvet_assembly.fa')
        shutil.copy('data/Rhodo_velvet_assembly.fa', assembly_file_path)
        au = AssemblyUtil(os.environ['SDK_CALLBACK_URL'])
        cls.assembly_ref = au.save_assembly_from_fasta({
            'workspace_name': cls.wsName,
            'assembly_name': 'ecoli.assembly',
            'file': {'path': assembly_file_path}
        })

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.ws.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def getWsClient(self):
        return self.__class__.ws

    def getWsName(self):
        if hasattr(self.__class__, 'wsName'):
            return self.__class__.wsName
        suffix = int(time.time() * 1000)
        wsName = "test_GenomeFileUtil_" + str(suffix)
        ret = self.getWsClient().create_workspace({'workspace': wsName})
        self.__class__.wsName = wsName
        return wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    def test_upgrade_genome(self):
        old_rhodobacter = json.load(open('data/Rhodo_SPAdes_RAST.json'))
        old_rhodobacter['assembly_ref'] = self.assembly_ref
        new_rhodobacter = self.genome_interface._update_genome(old_rhodobacter)
        json.dump(new_rhodobacter, open(self.cfg['scratch']+'/new_genome', 'w'))
        save_info = {
            'workspace': self.getWsName(),
            'objects': [{
                'type': 'NewTempGenomes.Genome',
                'data': new_rhodobacter,
                'name': 'rhodobacter'
            }]
        }
        result = self.ws.save_objects(save_info)

