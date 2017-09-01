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
        cls.wsName = "Phytozome_Genomes"
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

    def test_phytozome_to_genome(self):
        
        #Read Species Names
        Species_Names_File = "/kb/module/data/Phytozome_Names.txt"
        Species_Names_Dict={}
        input_file_handle = open(Species_Names_File,'rb')
        current_line = input_file_handle.readline()
        while ( current_line != '' ):
            current_line=current_line.strip()
            array = current_line.split('\t')
            Species_Names_Dict[array[0]]=array[1]
            current_line = input_file_handle.readline()

        Species_List = os.listdir(self.dtn_phytozome_root)

        for Species in Species_List:
            if(Species != "Cpapaya"):
                continue
            Species_Dir = os.path.join(self.dtn_phytozome_root,Species)
            Annotation_Dir = os.path.join(Species_Dir,"annotation")    
            Assembly_Dir = os.path.join(Species_Dir,"assembly")    

            if(os.path.isdir(Annotation_Dir)==False or os.path.isdir(Assembly_Dir)==False):
#                print "Skipping: "+Species
                continue

            Files_List  = os.listdir(Species_Dir)
            GeneModel_Version = ""
            for File in Files_List:
                if('.readme.txt' in File):
                    #Special Exception for Zmays
                    if(Species == "Zmays" and "MaizeSequence" in File):
                        continue

                    File=File.replace('.readme.txt','')
                    array = File.split('_')
                    GeneModel_Version=array[-1]
                    break

            Annotation_List = os.listdir(Annotation_Dir)
            Annotation_File = ""
            for File in Annotation_List:
                if(GeneModel_Version in File and '.gene.gff3.gz' in File):
                    Annotation_File = File
                    break

            Assembly_List = os.listdir(Assembly_Dir)
            Assembly_File = ""
            for File in Assembly_List:
                if('.fa.gz' in File):
                    if('masked' in File):
                        continue
                    Assembly_File = File
                    break

            Genome_Name = Species+"_PhytozomeV11_"+GeneModel_Version
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
                'scientific_name': Species_Names_Dict[Species]
            }

            result = self.getImpl().fasta_gff_to_genome(self.getContext(), input_params)[0]
