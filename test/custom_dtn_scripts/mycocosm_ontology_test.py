# -*- coding: utf-8 -*-
import unittest
import os  # noqa: F401
import json  # noqa: F401
import time
import shutil
import re
import gzip
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

        cls.dtn_root = "/kb/module/data/"
        cls.dtn_phytozome_root = cls.dtn_root+"PhytozomeV11_[PhytozomeV11]/"
        cls.dtn_mycocosm_root = cls.dtn_root+"Fungi_[fungi]/"
        cls.dtn_custom_root = cls.dtn_root+"custom/"

    def compile_ontology(self, Annotation_File, Identifier_Column=1, Ontology_Column=9):
        annotations = dict()
        #hardcoded header for now
        annotation_header=[]

        input_file_handle = gzip.open(Annotation_File,'rb')
        current_line = input_file_handle.readline()
        while ( current_line != '' ):
            current_line=current_line.strip()

            if(current_line.startswith("#pacId")):
                #Store header
                annotation_header=current_line.split('\t')
            else:
                annotation_items=current_line.split('\t')

                #Skip empty lines
                if(len(annotation_items) <= 1 or len(annotation_items)<=Ontology_Column):
                    current_line = input_file_handle.readline()
                    continue

                #Skip empty ontology
                if(annotation_items[Ontology_Column]==""):
                    current_line = input_file_handle.readline()
                    continue

                annotation_dict=dict()
                for entry in annotation_items[Ontology_Column].split(","):
                    if(entry == ''):
                        continue
                    entry=entry.replace("GO:GO:","GO:")
                    annotation_dict[entry]=1
                annotations[annotation_items[Identifier_Column]]=annotation_dict

            current_line = input_file_handle.readline()
        return annotations

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

        Mycocosm_Log = open(self.dtn_root+"Mycocosm_Ontologies.log","w")

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
            time_string = str(datetime.datetime.fromtimestamp(time.time()).strftime('%Y_%m_%d_%H_%M_%S'))
            Mycocosm_Log.write(time_string+":"+Species+": Loading Ontologies\n")

            Folder = Species_Names_Dict[Species]['folder']
            Species_Dir = os.path.join(self.dtn_mycocosm_root,Folder,"Files")
            Annotation_Dir = os.path.join(Species_Dir,"Annotation","Filtered_Models___best__","Functional_Annotations","GO")
            Annotation_List = os.listdir(Annotation_Dir)
            Annotation_Files = []
            for File in Annotation_List:
                if('GO' in File or 'Gene' in File):
                    Annotation_Files.append(File)

            if(len(Annotation_Files)==0):
                print "No ontology files: "+Species
                continue

            if(len(Annotation_Files)>1):
                print "Too many ontology files: "+Species,Annotation_Files
                continue

            Annotation_File=os.path.join(Annotation_Dir,Annotation_Files[0])

            #Read Header to try and log possible identifiers
            ontologies_file_handle = gzip.open(Annotation_File, 'rb')
            current_line = ontologies_file_handle.readline()
            while ( current_line != '' ):
                current_line=current_line.strip()
                items = current_line.strip('\t')
                if(len(items) <= 1):
                    current_line = ontologies_file_handle.readline()
                    continue
                else:
                    break

            ontologies_headers = current_line.split('\t')
            Mycocosm_Log.write(time_string+":"+Species+": Ontologies Headers: "+"|".join(ontologies_headers)+"\n")
            Mycocosm_Log.write(time_string+":"+Species+": Utilizing Headers: "+"|".join([ontologies_headers[0],
                                                                                         ontologies_headers[4]])+"\n")


            Ontology = self.compile_ontology(Annotation_File,0,4)
            time_string = str(datetime.datetime.fromtimestamp(time.time()).strftime('%Y_%m_%d_%H_%M_%S'))
            Mycocosm_Log.write(time_string+":"+Species+": Ontologies Parsed\n")

            #Retrieve OntologyDictionary
            Ontology_Dictionary = self.dfu.get_objects({'object_refs':["KBaseOntology/gene_ontology"]})['data'][0]['data']['term_hash']

            #Load Genome Object
            Genome_Name = Species
            Genome_Result = self.dfu.get_objects({'object_refs':[self.wsName+'/'+Genome_Name]})['data'][0]
            Genome_Object = Genome_Result['data']
            Genome_Meta = Genome_Result['info'][10]
            Workspace_ID = Genome_Result['info'][6]

            #Iterate through Features
            Found_Ftrs_Count=list()
            for feature in Genome_Object['features']:
                feature["ontology_terms"]=dict()

                if(feature['id'] in Ontology):
                    Found_Ftrs_Count.append(feature['id'])
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

            time_string = str(datetime.datetime.fromtimestamp(time.time()).strftime('%Y_%m_%d_%H_%M_%S'))
            Mycocosm_Log.write(time_string+":"+Species+": "+str(len(Found_Ftrs_Count))+" Features Loaded with Ontologies\n")
            Mycocosm_Log.write(time_string+":"+Species+": Feature Sample IDs: "+"|".join(Found_Ftrs_Count[:2])+"\n")

            print "Saving: "+Genome_Name+" with ontology for "+str(len(Found_Ftrs_Count))+ " features"
            print "With Meta: ",Genome_Meta
            self.dfu.save_objects({'id':Workspace_ID, 'objects' : [ {'type': 'KBaseGenomes.Genome', 'data': Genome_Object,
                                                                     'meta' : Genome_Meta, 'name' : Genome_Name} ]})
            time_string = str(datetime.datetime.fromtimestamp(time.time()).strftime('%Y_%m_%d_%H_%M_%S'))
            Mycocosm_Log.write(time_string+":"+Species+": Saved with Ontologies\n")
