# -*- coding: utf-8 -*-
import unittest
import os  # noqa: F401
import json  # noqa: F401
import time
import shutil
import re
import datetime
import simplejson
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
        cls.dtn_custom_root = "/kb/module/data/custom/"

    def test_find_taxonomy(self):

        #Retrieve OntologyDictionary
        Taxonomy_Lookup = self.dfu.get_objects({'object_refs':["ReferenceTaxons/taxon_lookup"]})['data'][0]['data']['taxon_lookup']

    tax_id=0
    taxon_object_name = "unknown_taxon"
    display_sc_name = None
    if(full_organism[0:3] in taxon_lookup and full_organism in taxon_lookup[full_organism[0:3]]):
        tax_id=taxon_lookup[full_organism[0:3]][full_organism]
        taxon_object_name = "%s_taxon" % (str(tax_id))

    taxon_info = ws_client.get_objects([{"workspace": args.taxon_wsname,
                                         "name": taxon_object_name}])[0]

    taxon_ref = "%s/%s/%s" % (taxon_info['info'][6], taxon_info['info'][0], taxon_info['info'][4])
    display_sc_name = taxon_info['data']['scientific_name']
    taxonomy = taxon_info['data']['scientific_lineage']
    core_genome_name = "%s_%s_%s" % (args.organism,args.source,args.release)

        Custom_Annotation=self.dtn_custom_root+"GO_terms_Ptr_V2.2.txt"
        Custom_Genome="sunita:narrative_1500591494735/PtrichocarpaV2.2_annotated"
        Custom_Genome="Phytozome_Genomes/Ptrichocarpa_PhytozomeV11_v2.2"

        #Parse Annotation File
        Ontology = self.compile_ontology(Custom_Annotation,0,2)

        #Retrieve OntologyDictionary
        Ontology_Dictionary = self.dfu.get_objects({'object_refs':["KBaseOntology/gene_ontology"]})['data'][0]['data']['term_hash']

        #Load Genome Object
        Genome_Result = self.dfu.get_objects({'object_refs':[Custom_Genome]})['data'][0]
        Genome_Object = Genome_Result['data']
        Genome_Meta = Genome_Result['info'][10]
        Workspace_ID = Genome_Result['info'][6]

        time_string = str(datetime.datetime.fromtimestamp(time.time()).strftime('%Y_%m_%d_%H_%M_%S'))

        #Iterate through Features
        for feature in Genome_Object['features']:
            feature["ontology_terms"]=dict()
            
            if(feature['id'] in Ontology):
                ontology_terms = dict()                    
                ontology_terms["GO"]=dict()
                for Ontology_Term in Ontology[feature["id"]].keys():
                    if(Ontology_Term not in Ontology_Dictionary):
                        print "Missing: "+Ontology_Term
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
        genome_file = open(self.dtn_custom_root+'/Custom.json', 'w+')
        genome_file.write(genome_string)
        genome_file.close()

        print "Saving: "+Custom_Genome
        print "With Meta: ",Genome_Meta

        Genome_Name=Custom_Genome.split('/')[1]
        print "Saving as: "+Genome_Name
        self.dfu.save_objects({'id':Workspace_ID, 'objects' : [ {'type': 'KBaseGenomes.Genome', 'data': Genome_Object,
                                                                 'meta' : Genome_Meta, 'name' : Genome_Name} ]})
