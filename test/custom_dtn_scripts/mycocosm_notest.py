# -*- coding: utf-8 -*-
import unittest
import os  # noqa: F401
import json  # noqa: F401
import time
import shutil
import re
import sys

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
        cls.wsName = "Mycocosm_Genomes"
        cls.prepare_data()

#    @classmethod
#    def tearDownClass(cls):
#        if hasattr(cls, 'wsName'):
#            cls.wsClient.delete_workspace({'workspace': cls.wsName})
#            print('Test workspace was deleted')

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
#        cls.wsClient.create_workspace({'workspace': cls.wsName})

        cls.dtn_phytozome_root = "/kb/module/data/PhytozomeV11_[PhytozomeV11]/"
        cls.dtn_custom_root = "/kb/module/data/custom/"
        cls.dtn_mycocosm_root = "/kb/module/data/Fungi_[fungi]/"

    def test_phytozome_to_genome(self):
        
        #Read Species Names
        Species_Names_File = "/kb/module/data/Mycocosm_Names.txt"
        Species_Names_Dict={}
        input_file_handle = open(Species_Names_File,'rb')
        current_line = input_file_handle.readline()
        while ( current_line != '' ):
            current_line=current_line.strip()
            array = current_line.split('\t')
            Species_Names_Dict[array[0]]={'name':array[1],'folder':array[2]}
            current_line = input_file_handle.readline()

        for Species in sorted(Species_Names_Dict):
            Folder = Species_Names_Dict[Species]['folder']
            Species_Dir = os.path.join(self.dtn_mycocosm_root,Folder,"Files")
            Annotation_Dir = os.path.join(Species_Dir,"Annotation","Filtered_Models___best__","Genes") 
            Assembly_Dir = os.path.join(Species_Dir,"Assembly","Assembled_scaffolds__unmasked_")

            if(os.path.isdir(Annotation_Dir)==False or os.path.isdir(Assembly_Dir)==False):
                print "Skipping without annotation or assembly: "+Species_Dir
                continue

            Annotation_List = os.listdir(Annotation_Dir)
            Annotation_Files = []
            for File in Annotation_List:
                if('gff' not in File):
                    continue

                if('genes' in File):
                    Annotation_Files.append(File)

            if(len(Annotation_Files)==0):
                for File in Annotation_List:
                    if('gff' not in File):
                        continue

                    if('Gene' in File and 'protein' not in File):
                        Annotation_Files.append(File)

            if(len(Annotation_Files)>1):
                print "Too many annotation files: "+Species,Annotation_Files
                continue

            Annotation_File=Annotation_Files[0]

            Assembly_List = os.listdir(Assembly_Dir)
            Assembly_Files = []
            for File in Assembly_List:
                if('fasta' not in File or 'Mito' in File):
                    continue
                
                Assembly_Files.append(File)

            if(len(Assembly_Files)>1):
                print "Too many assembly files: "+Species,Assembly_Files
                continue

            Assembly_File=Assembly_Files[0]

            Genome_Name = Species

            #Special case for Zmays
            Genome_Name = Genome_Name.replace("+","")

            Fa_Path = os.path.join(Assembly_Dir,Assembly_File)
            Gff_Path = os.path.join(Annotation_Dir,Annotation_File)

            input_params = {
                'fasta_file': {'path': Fa_Path},
                'gff_file': {'path': Gff_Path},
                'genome_name': Genome_Name,
                'workspace_name': self.getWsName(),
                'source': 'JGI',
                'type': 'Reference',
                'scientific_name': Species_Names_Dict[Species]['name']
            }

            print input_params
            result = self.getImpl().fasta_gff_to_genome(self.getContext(), input_params)[0]
            genome_info = result['genome_info']
            print Species,genome_info[10]
            break
