# -*- coding: utf-8 -*-
import unittest
import os  # noqa: F401
import json  # noqa: F401
import time
import shutil
import re
import sys
import datetime
import collections
#import simplejson

from os import environ
try:
        from ConfigParser import ConfigParser  # py2
except:
	from configparser import ConfigParser  # py3

from pprint import pprint  # noqa: F401

from GenomeFileUtil.GenomeFileUtilImpl import GenomeFileUtil
from GenomeFileUtil.core.GenomeInterface import GenomeInterface
from GenomeFileUtil.GenomeFileUtilImpl import SDKConfig
from GenomeFileUtil.GenomeFileUtilServer import MethodContext
from GenomeFileUtil.core.FastaGFFToGenome import FastaGFFToGenome
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace as workspaceService

class FastaGFFToGenomeUploadTest(unittest.TestCase):

        @classmethod
        def setUpClass(cls):
                print('setting up class')
		token = environ.get('KB_AUTH_TOKEN', None)
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
		cls.wsClient = workspaceService(cls.wsURL, token=token)
		cls.serviceImpl = GenomeFileUtil(cls.cfg)

		cls.dfu = DataFileUtil(os.environ['SDK_CALLBACK_URL'], token=token)
		cls.scratch = cls.cfg['scratch']
		cls.shockURL = cls.cfg['shock-url']
		cls.gfu_cfg = SDKConfig(cls.cfg)
		cls.wsName = "Phytozome_Genomes"
#		 cls.wsName = "Athaliana_Test"
#		 cls.wsName = "sunita:narrative_1568652277639"
		cls.prepare_data()

#	 @classmethod
#	 def tearDownClass(cls):
#		 if hasattr(cls, 'wsName'):
#			 cls.wsClient.delete_workspace({'workspace': cls.wsName})
#			 print('Test workspace was deleted')

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
#		 cls.wsClient.create_workspace({'workspace': cls.wsName})

		cls.dtn_root = "/kb/module/genomes/Phytozome/"

	def compile_ontology(self, Annotation_File, Identifier_Column):
		annotations = dict()
		#hardcoded header for now
		annotation_header=[]
		ontology_column=9

		if not os.path.isfile(Annotation_File):
			return annotations

		with open(Annotation_File) as f:
			for line in f:
				line=line.strip()
				if(line.startswith("#pacId")):
					#Store header
					annotation_header=line.split('\t')
				else:
					annotation_items=line.split('\t')

					#Skip empty lines
					if(len(annotation_items) <= 1 or len(annotation_items)<=ontology_column):
						continue

					#Skip empty ontology
					if(annotation_items[ontology_column]==""):
						continue

					annotation_dict=dict()
					for entry in annotation_items[ontology_column].split(","):
						if(entry == ''):
							continue
						entry=entry.replace("GO:GO:","GO:")
						annotation_dict[entry]=1
					annotations[annotation_items[Identifier_Column]]=annotation_dict
		return annotations

	def compile_functions(self, Functions_File, Identifier_Column=0,Functions_Column=1,EC_Column=-1):
		functions = dict()

		if not os.path.isfile(Functions_File):
			return functions

		with open(Functions_File) as f:
			for line in f:
				line=line.strip()
				function_items=line.split('\t')

				if(len(function_items) <= Functions_Column):
					Functions_Column -= 1

				if(function_items[Functions_Column] == ""):
					continue

				Function = function_items[Functions_Column]

				if(EC_Column != -1):
					Function+=" (EC "+function_items[EC_Column]+")"

				if(function_items[Identifier_Column] not in functions):
					functions[function_items[Identifier_Column]]=dict()
				functions[function_items[Identifier_Column]][Function]=1

		return functions

	def compile_synonyms(self, Synonyms_File, Identifier_Column=0, Synonym_Column=1):
		synonyms=dict()

		if not os.path.isfile(Synonyms_File):
			return synonyms

		with open(Synonyms_File) as f:
			for line in f:
				line=line.strip()
				synonyms_items=line.split('\t')
				if(len(synonyms_items) <= Synonym_Column or synonyms_items[Synonym_Column] == ""):
					continue

				Synonym = synonyms_items[Synonym_Column]

				if(synonyms_items[Identifier_Column] not in synonyms):
					synonyms[synonyms_items[Identifier_Column]]=dict()
				synonyms[synonyms_items[Identifier_Column]][Synonym]=1
				
		return synonyms

	def test_phytozome_to_genome(self):

		#Read Species Names
		SN_File = os.path.join('data','Phytozome_Names.txt')
		Species_Names_Dict={}
		with open(SN_File) as f:
			for line in f:
				line=line.strip()
				array = line.split('\t')
				Species_Names_Dict[array[0]]=array[1]

		GM_File = os.path.join('data','Accepted_Phytozome_Versions_GeneModels.txt')
		Species_Dict=dict()
		with open(GM_File) as f:
			for line in f:
				line=line.strip('\r\n')
				(phytozome_release,species,phytozome_identifier,genemodel_version,tax_id)=line.split('\t')

				if(species not in Species_Dict):
					Species_Dict[species]=dict()
				Species_Dict[species][genemodel_version]={'release':phytozome_release,'identifier':phytozome_identifier,'tax_id':tax_id}

		Species_Versions_to_Skip=dict()
		with open(os.path.join('data','Phytozome_Upload_Summary.txt')) as f:
			for line in f:
				line=line.strip()
				array = line.split('\t')
		#		 array = array[0].split(':')
		#		 Species_Version = array[1].strip()
				Species_Version=array[0]
				Species_Versions_to_Skip[Species_Version]=1

		summary_file=open(os.path.join('data','Phytozome_Upload_Summary.txt'),'w')
		for species in sorted(Species_Dict):

			for version in sorted(Species_Dict[species]):

				Genome_Name = species+"_"+version

				if(Genome_Name in Species_Versions_to_Skip):
					continue

				phytozome_version=Species_Dict[species][version]['release']

				sp=species
				#Special exception: Zmays PH207
				if(species == "Zmays" and "PH207" in version):
					sp=species+version.split('_')[0]

				path = os.path.join(self.dtn_root,'Phytozome'+phytozome_version,sp)

				if(os.path.isdir(path) is not True):
					print("Path Not Found: "+path)
					continue

				has_assembly=0
				for files in os.listdir(path):
					if(files=="assembly"):
						has_assembly=1

				if(has_assembly == 0):
					for versions in os.listdir(path):
						if(versions == version):
							for files in os.listdir(os.path.join(path,versions)):
								if(files=="assembly"):
									path=os.path.join(path,versions)
									has_assembly=1

				#Assembly file retrieval, should only find one, if any
				assembly_file = os.listdir(os.path.join(path,'assembly'))[0]

				#Annotation file retrieval, at least one, maybe two or three
				gff_file = ""
				functions_file = ""
				ontology_file = ""
				names_file = ""
				for ann_file in os.listdir(os.path.join(path,'annotation')):
					if('gene' in ann_file):
						gff_file = ann_file
					elif('defline' in ann_file):
						functions_file = ann_file
					elif('info' in ann_file):
						ontology_file = ann_file
					elif('synonym' in ann_file):
						names_file = ann_file

				Fa_Path = os.path.join(path,'assembly',assembly_file)
				Gff_Path = os.path.join(path,'annotation',gff_file)

				phytozome_version = phytozome_version.split('_')[0]
				tax_id = Species_Dict[species][version]['tax_id']
				input_params = {'fasta_file': {'path': Fa_Path},
								'gff_file': {'path': Gff_Path},
								'genome_name': Genome_Name,
								'workspace_name': self.getWsName(),
								'source': 'JGI Phytozome '+phytozome_version,
								'source_id' : version,
								'type': 'Reference',
								'scientific_name': Species_Names_Dict[species],
								'taxon_id': tax_id,
								'genetic_code':1} 
				result = self.getImpl().fasta_gff_to_genome(self.getContext(), input_params)[0]
				print(result)
				
				#Load Genome Object
				Genome_Result = self.dfu.get_objects({'object_refs':[self.wsName+'/'+Genome_Name]})['data'][0]
				Genome_Object = Genome_Result['data']

				############################################################
				# Functions
				###########################################################

				Functions_Path = os.path.join(path,'annotation',functions_file)
				Functions = self.compile_functions(Functions_Path,0,2,1)
				print("Functions compiled")
				summary_file.write(Genome_Name+'\t'+functions_file+'\t'+str(len(list(Functions)))+'\n')
				
				Found_Count={'features':0,'mrnas':0,'cdss':0}
				if(len(list(Functions))>0):
					for key in Found_Count:
						print("Searching: "+key+"\t"+Genome_Object[key][0]['id'])
						for entity in Genome_Object[key]:
							if(entity['id'] in Functions):
								entity["functions"]=sorted(Functions[entity['id']].keys())
								Found_Count[key]+=1

					#If no features were annotated, and mrnas were annotated, use parent_gene to do transfer annotation
					parent_feature_functions = collections.defaultdict(dict)
					if(Found_Count['features']==0 and Found_Count['mrnas']!=0):
						#Create lookup dict
						parent_feature_index = dict([(f['id'], i) for i, f in enumerate(Genome_Object['features'])])
						for mrna in Genome_Object['mrnas']:
							if('functions' in mrna):
								parent_feature = parent_feature_index[mrna['parent_gene']]
								for function in mrna['functions']:
									parent_feature_functions[parent_feature][function]=1

					for index in parent_feature_functions:
						Genome_Object['features'][index]['functions']=sorted(parent_feature_functions[index].keys())
						Found_Count['features']+=1

					summary_file.write(Genome_Name+'\t'+functions_file+'\t'+str(Found_Count)+'\n')

				############################################################
				# Ontology
				###########################################################

				#Parse Annotation File
				Annotation_Path = os.path.join(path,'annotation',ontology_file)
				Feature_Ontology = self.compile_ontology(Annotation_Path,1)
				mRNA_Ontology = self.compile_ontology(Annotation_Path,2)
				print("Ontology compiled")
				summary_file.write(Genome_Name+'\t'+ontology_file+'\t'+str(len(Feature_Ontology.keys()))+'\n')
				summary_file.write(Genome_Name+'\t'+ontology_file+'\t'+str(len(mRNA_Ontology.keys()))+'\n')

				#Retrieve OntologyDictionary
				Ontology_Dictionary = self.dfu.get_objects({'object_refs':["KBaseOntology/gene_ontology"]})['data'][0]['data']['term_hash']
				time_string = str(datetime.datetime.fromtimestamp(time.time()).strftime('%Y_%m_%d_%H_%M_%S'))

				Found_Count={'features':0,'mrnas':0,'cdss':0}
				if(len(Feature_Ontology.keys())!=0 or len(mRNA_Ontology.keys())!=0):
					for key in Found_Count:
						print("Searching: "+key+"\t"+Genome_Object[key][0]['id'])
						for entity in Genome_Object[key]:
							if(entity['id'] in Feature_Ontology):
								ontology_terms = dict()				
								ontology_terms["GO"]=dict()
								for Ontology_Term in Feature_Ontology[entity["id"]].keys():
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
								entity["ontology_terms"]=ontology_terms
								Found_Count[key]+=1

							if(entity['id'] in mRNA_Ontology):
								ontology_terms = dict()				
								ontology_terms["GO"]=dict()
								for Ontology_Term in mRNA_Ontology[entity["id"]].keys():
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
								entity["ontology_terms"]=ontology_terms
								Found_Count[key]+=1

				summary_file.write(Genome_Name+'\t'+ontology_file+'\t'+str(Found_Count)+'\n')

				############################################################
				# Synonyms
				###########################################################

				Synonyms_Path = os.path.join(path,'annotation',names_file)
				Synonyms = self.compile_synonyms(Synonyms_Path,0,1)
				print("Synonyms compiled")
				summary_file.write(Genome_Name+'\t'+names_file+'\t'+str(len(list(Synonyms)))+'\n')

				Found_Count={'features':0,'mrnas':0,'cdss':0}
				if(len(list(Synonyms))>0):
					for key in Found_Count:
						print("Searching: "+key+"\t"+Genome_Object[key][0]['id'])
						for entity in Genome_Object[key]:
							if(entity['id'] in Synonyms):
								if("aliases" not in entity):
									entity["aliases"]=list()

								for synonym in sorted(Synonyms[entity['id']]):
									entity["aliases"].append(["JGI",synonym])

								Found_Count[key]+=1

					#If no features were annotated, and mrnas were annotated, use parent_gene to do transfer annotation
					parent_feature_synonyms = collections.defaultdict(dict)
					if(Found_Count['features']==0 and Found_Count['mrnas']!=0):
						#Create lookup dict
						parent_feature_index = dict([(f['id'], i) for i, f in enumerate(Genome_Object['features'])])
						for mrna in Genome_Object['mrnas']:
							if(entity['id'] in Synonyms):
								if("aliases" not in entity):
									entity["aliases"]=list()
								for synonym in sorted(Synonyms[entity['id']]):
									entity["aliases"].append(["JGI",synonym])

							if('aliases' in mrna):
								parent_feature = parent_feature_index[mrna['parent_gene']]
								for synonym in mrna['aliases']:
									parent_feature_synonyms[parent_feature][synonym[1]]=1

						for index in parent_feature_synonyms:
							if("aliases" not in Genome_Object['features'][index]):
								Genome_Object['features'][index]['aliases']=list()

							for synonym in sorted(parent_feature_synonyms[index].keys()):
								Genome_Object['features'][index]['aliases'].append(("JGI",synonym))

						Found_Count['features']+=1

				summary_file.write(Genome_Name+'\t'+names_file+'\t'+str(Found_Count)+'\n')

				############################################################
				# Saving
				###########################################################

				#Save Genome Object
				#genome_string = simplejson.dumps(Genome_Object, sort_keys=True, indent=4, ensure_ascii=False)
				#genome_file = open(self.scratch+'/'+Genome_Name+'.json', 'w+')
				#genome_file.write(genome_string)
				#genome_file.close()

				#Retaining metadata
				Genome_Meta = Genome_Result['info'][10]
				Workspace_ID = Genome_Result['info'][6]

				save_result = self.getImpl().save_one_genome(self.getContext(), {'workspace' : self.wsName,
																				 'name' : Genome_Name,
																				 'data' : Genome_Object,
																				 'upgrade' : 1})

				#Saving metadata
				Genome_Result = self.dfu.get_objects({'object_refs':[self.wsName+'/'+Genome_Name]})['data'][0]
				Genome_Object = Genome_Result['data']
				self.dfu.save_objects({'id':Workspace_ID, 'objects' : [ {'type': 'KBaseGenomes.Genome', 'data': Genome_Object,
																		 'meta' : Genome_Meta, 'name' : Genome_Name} ]})
