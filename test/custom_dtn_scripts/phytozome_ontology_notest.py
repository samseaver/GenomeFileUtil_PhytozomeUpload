# -*- coding: utf-8 -*-
import unittest
import os  # noqa: F401
import json  # noqa: F401
import time
import shutil
import re
import datetime
import simplejson

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
#        cls.wsName = "Phytozome_Genomes"
        cls.wsName = "sunita:narrative_1500591494735"
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

    def compile_ontology(self, Annotation_File):
        annotations = dict()
        #hardcoded header for now
        annotation_header=[]

        identifier_column=1
        ontology_column=9

        input_file_handle = open(Annotation_File,'rb')
        current_line = input_file_handle.readline()
        while ( current_line != '' ):
            current_line=current_line.strip()

            if(current_line.startswith("#pacId")):
                #Store header
                annotation_header=current_line.split('\t')
            else:
                annotation_items=current_line.split('\t')

                #Skip empty lines
                if(len(annotation_items) <= 1 or len(annotation_items)<=ontology_column):
                    current_line = input_file_handle.readline()
                    continue

                #Skip empty ontology
                if(annotation_items[ontology_column]==""):
                    current_line = input_file_handle.readline()
                    continue

                annotation_dict=dict()
                for entry in annotation_items[ontology_column].split(","):
                    if(entry == ''):
                        continue
                    entry=entry.replace("GO:GO:","GO:")
                    annotation_dict[entry]=1
                annotations[annotation_items[identifier_column]]=annotation_dict

            current_line = input_file_handle.readline()
        return annotations

    def compile_functions(self, files):
        functions = dict()
        if 'input_function_file' in files and files['input_function_file'] is not None:
            eprint("Compiling functions")
            input_file_handle = gzip.open(files['input_function_file'],'rb')
            current_line = input_file_handle.readline()
            while ( current_line != '' ):
                current_line=current_line.strip()

                #possible problems with unicode characters, re-encode
                #replace mis-encoded character
                #0x96
                try:
                    current_line=unicode(current_line, errors='strict')
                except UnicodeDecodeError:
                    old_line=current_line
                    current_line=""
                    for char in old_line:
                        if(char == "\x96"):
                            char="-"
                        current_line+=char

                function_items=current_line.split('\t')

                if(len(function_items) <= 1 or function_items[1] == ""):
                    current_line = input_file_handle.readline()
                    continue

                functions[function_items[0]]=function_items[1]
                current_line = input_file_handle.readline()

        return functions

    def compile_aliases(self, files):
        aliases = dict()
        if 'input_alias_file' in files and files['input_alias_file'] is not None:
            eprint("Compiling aliases")
            input_file_handle = gzip.open(files['input_alias_file'],'rb')
            current_line = input_file_handle.readline()
            while ( current_line != '' ):
                current_line=current_line.strip()
                aliases_items=current_line.split('\t')

                if(len(aliases_items) <= 1 or aliases_items[1] == ""):
                    current_line = input_file_handle.readline()
                    continue

                identifier = aliases_items.pop(0)
                if(identifier not in aliases):
                    aliases[identifier]=[]
                for alias in aliases_items:
                    aliases[identifier].append(alias)

                current_line = input_file_handle.readline()

        return aliases

    def compile_proteins(self, files):
        proteins = dict()
        if 'input_protein_file' in files and files['input_protein_file'] is not None:
            eprint("Compiling proteins")
            input_file_handle = gzip.open(files['input_protein_file'],'rb')
            faiter = (x[1] for x in itertools.groupby(input_file_handle, lambda line: line[0] == ">"))
            for header in faiter:
                # drop the ">"
                header = header.next()[1:].strip()
                # join all sequence lines to one.
                seq = "".join(s.strip() for s in faiter.next())

                try:
                    fasta_header,fasta_description = header.split(' ',1)
                except:
                    fasta_header = header
                    fasta_description = None

                protein_alias=fasta_header
                description_dict=dict()
                if(fasta_description is not None):
                    description_dict=dict(x.split('=') for x in fasta_description.split(' '))

                if('transcript' in description_dict.keys()):
                    fasta_header=description_dict['transcript']

                seq=seq.upper()
                proteins[fasta_header]={'protein_id':protein_alias,'sequence':seq}

        return proteins

    def test_ontology_to_genome(self):

        Custom_Annotation=""
        Custom_Annotation=self.dtn_custom_root+"GO_terms_Ptr_V2.2.txt"

#        use cls.dfu

        Species_List = os.listdir(self.dtn_phytozome_root)

        for Species in Species_List:
            if(Species != "Ptrichocarpa"):
                continue

            Species_Dir = os.path.join(self.dtn_phytozome_root,Species)

            Files_List  = os.listdir(Species_Dir)
            GeneModel_Version = ""
            for File in Files_List:
                if('.readme.txt' in File):
                    #Special Exception for Zmays
                    if(Species == "Zmays" and "MaizeSequence" in File):
                        continue

                    File=File.replace('.readme.txt','')
                    array = File.split('_')
                    print array
                    GeneModel_Version=array[-1]
                    break

            Genome_Name = Species+"_PhytozomeV11_"+GeneModel_Version

            print Genome_Name

            if(Genome_Name != "Ptrichocarpa_PhytozomeV11_v2.2"):
                continue

            Annotation_File = ""
            if(Custom_Annotation != ""):

                Annotation_Dir = os.path.join(Species_Dir,"annotation")    

                #Skip species altogether if no annotation dir present
                if(os.path.isdir(Annotation_Dir)==False):
                    continue

                Annotation_List = os.listdir(Annotation_Dir)
                for File in Annotation_List:
                    if(GeneModel_Version in File and '.annotation_info.txt' in File):
                        Annotation_File = File
                        break

                Annotation_File = os.path.join(Annotation_Dir,Annotation_File)
            else:
                Annotation_File = Custom_Annotation

            #Parse Annotation File
            Ontology = self.compile_ontology(Annotation_File)

            print "Ontology compiled"

            #Retrieve OntologyDictionary
            Ontology_Dictionary = self.dfu.get_objects({'object_refs':["KBaseOntology/gene_ontology"]})['data'][0]['data']['term_hash']

            #Load Genome Object
            Genome_Result = self.dfu.get_objects({'object_refs':['Phytozome_Genomes/'+Genome_Name]})['data'][0]
            Genome_Object = Genome_Result['data']
            Genome_Meta = Genome_Result['info'][10]

            time_string = str(datetime.datetime.fromtimestamp(time.time()).strftime('%Y_%m_%d_%H_%M_%S'))

            #Iterate through Features
            for feature in Genome_Object['features']:
                feature["ontology_terms"]=dict()

                if(feature['id'] in Ontology):
                    ontology_terms = dict()                    
                    ontology_terms["GO"]=dict()
                    for Ontology_Term in Ontology[feature["id"]].keys():
                        if(Ontology_Term not in Ontology_Dictionary):
                            continue
                        if(Ontology_Term not in ontology_terms["GO"]):
                            OntologyEvidence=[{"method":"GFF_Fasta_Genome_to_KBaseGenomes_Genome",
                                               "timestamp":time_string,"method_version":"1.0"},
                                              {"method":"Phytozome annotation_info.txt",
                                               "timestamp":time_string,"method_version":"11"}]
                            OntologyData={"id":Ontology_Term,
                                          "ontology_ref":"KBaseOntology/gene_ontology",
                                          "term_name":Ontology_Dictionary[Ontology_Term]["name"],
                                          "term_lineage":[],
                                          "evidence":OntologyEvidence}
                            ontology_terms["GO"][Ontology_Term]=OntologyData
                    feature["ontology_terms"]=ontology_terms

            #Save Genome Object
            genome_string = simplejson.dumps(Genome_Object, sort_keys=True, indent=4, ensure_ascii=False)
            genome_file = open(self.scratch+'/'+Genome_Name+'.json', 'w+')
            genome_file.write(genome_string)
            genome_file.close()

            print "Saving: "+Genome_Name
            print "With Meta: ",Genome_Meta

#typedef structure {
#        int id;
#        list<ObjectSaveData> objects;
#    } SaveObjectsParams;

#            self.wsClient.save_objects({'workspace' : 'Phytozome_Genomes', 'objects' : {'type': 'KBaseGenomes.Genome',
#                                                                                        'data': Genome_Object,
#                                                                                        'name' : Genome_Name,
#                                                                                        'meta' : Genome_Meta}})
