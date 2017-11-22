import os
import sys
import shutil
import uuid
import time
import json
import gzip
import copy
import itertools
import hashlib
import collections
import datetime
import re

# KBase imports
from DataFileUtil.DataFileUtilClient import DataFileUtil
from AssemblyUtil.AssemblyUtilClient import AssemblyUtil
from KBaseReport.KBaseReportClient import KBaseReport
from GenomeInterface import GenomeInterface

# 3rd party imports
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable
import Bio.SeqIO

codon_table = CodonTable.ambiguous_generic_by_name["Standard"]

def log(message, prefix_newline=False):
    """Logging function, provides a hook to suppress or redirect log messages."""
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))


class FastaGFFToGenome:

    def __init__(self, config):
        self.cfg = config
        self.au = AssemblyUtil(config.callbackURL)
        self.dfu = DataFileUtil(self.cfg.callbackURL)
        self.gi = GenomeInterface(self.cfg)
        self.taxon_wsname = self.cfg.raw['taxon-workspace-name']
        self.time_string = str(datetime.datetime.fromtimestamp(
            time.time()).strftime('%Y_%m_%d_%H_%M_%S'))
        yml_text = open('/kb/module/kbase.yml').read()
        self.version = re.search("module-version:\n\W+(.+)\n", yml_text).group(
            1)
        self.feature_dict = {}
        self.ontologies_present = collections.defaultdict(dict)
        self.skiped_features = collections.Counter()
        self.feature_counts = collections.Counter()

    def import_file(self, params):

        # 1) validate parameters
        self._validate_import_file_params(params)

        # 2) construct the input directory staging area
        input_directory = os.path.join(self.cfg.sharedFolder, 'fast_gff_upload_'+str(uuid.uuid4()))
        os.makedirs(input_directory)
        file_paths = self._stage_input(params, input_directory)

        # 3) extract out the parameters
        params = self._set_parsed_params(params)

        # 4) do the upload
        result = self.upload_genome(
            input_fasta_file=file_paths["fasta_file"],
            input_gff_file=file_paths["gff_file"],

            workspace_name=params['workspace_name'],
            core_genome_name=params['genome_name'],
            scientific_name=params['scientific_name'],
            taxon_reference=params['taxon_reference'],
            source=params['source'],
            genome_type=params['type'],
            release=params['release']
        )

        # 5) generate report
        output_data_ref = params['workspace_name']+"/"+params['genome_name']
        reportObj = {'objects_created': [{'ref': output_data_ref,
                                          'description': 'KBase Genome object'}],
                     'text_message': result['report_string']}

        reportClient = KBaseReport(os.environ['SDK_CALLBACK_URL'])
        report_info = reportClient.create({'report': reportObj,
                                           'workspace_name': params['workspace_name']})

        # 6) clear the temp directory
        shutil.rmtree(input_directory)

        # 7) return the result
        info = result['genome_info']
        details = {
            'genome_ref': str(info[6]) + '/' + str(info[0]) + '/' + str(info[4]),
            'genome_info': info,
            'report_name': report_info['name'],
            'report_ref':  report_info['ref']
        }

        return details

    def upload_genome(self, input_gff_file=None, input_fasta_file=None,
                      workspace_name=None, core_genome_name=None,
                      scientific_name="unknown_taxon", taxon_reference=None,
                      source=None, release=None, genome_type=None):

        # reading in GFF file
        features_by_contig = self._retrieve_gff_file(input_gff_file)

        # parse feature information
        fasta_contigs = Bio.SeqIO.parse(input_fasta_file, "fasta")
        for contig in fasta_contigs:
            molecule_type = str(contig.seq.alphabet).replace(
                'IUPACAmbiguous', '').strip('()')
            for feature in features_by_contig.get(contig.id, []):
                if feature['type'] == 'exon':
                    self._parse_exon(feature)
                else:
                    self._transform_feature(contig, feature)

        # save assembly file
        assembly_ref = self.au.save_assembly_from_fasta(
            {'file': {'path': input_fasta_file},
             'workspace_name': workspace_name,
             'assembly_name': core_genome_name+".assembly"})
        assembly_data = self.dfu.get_objects(
            {'object_refs': [assembly_ref],
             'ignore_errors': 0})['data'][0]['data']

        # generate genome info
        genome = self._gen_genome_info(core_genome_name, scientific_name,
                                       assembly_ref, source, assembly_data,
                                       input_gff_file, molecule_type)

        workspace_id = self.dfu.ws_name_to_id(workspace_name)
        genome_info = self.dfu.save_objects(
            {"id": workspace_id,
             "objects": [{"name": core_genome_name,
                          "type": "NewTempGenomes.Genome",
                          "data": genome}]})[0]
        report_string = ''

        return {'genome_info': genome_info, 'report_string': report_string}

    def _parse_exon(self, feature):
        parent_id = feature.get('Parent', '')
        if not parent_id:
            raise ValueError("An exon has no valid parent: {}".format(feature))
        parent = self.feature_dict[parent_id]
        if parent.get('exons', False):
            self.feature_dict[parent_id]['location'].append(
                self._location(feature))
        else:
            parent['exons'] = True
            self.feature_dict[parent_id]['location'] = [self._location(feature)]

    @staticmethod
    def _location(in_feature):
        return [
            in_feature['contig'],
            in_feature['start'],
            in_feature['strand'],
            in_feature['end'] - in_feature['start']
        ]

    def _get_ontology(self, feature):
        ontology = collections.defaultdict(dict)
        for key in ("GO_process", "GO_function", "GO_component"):
            if key in feature['attributes']:
                sp = feature['attributes'][key][0][3:].split(" - ")
                ontology['GO'][sp[0]] = [1]
                self.ontologies_present['GO'][sp[0]] = sp[1]
        # TODO: Support other ontologies
        return dict(ontology)

    def _validate_import_file_params(self, params):
        """
        validate_import_file_params:
                    validates params passed to FastaGFFToGenome.import_file method

        """

        # check for required parameters
        for p in ['workspace_name', 'genome_name', 'fasta_file', 'gff_file']:
            if p not in params:
                raise ValueError('"{}" parameter is required, but missing'.format(p))

        # one and only one of 'path', or 'shock_id' is required
        for key in ('fasta_file', 'gff_file'):
            file = params[key]
            if not isinstance(file, dict):
                raise ValueError('Required "{}" field must be a map/dict'.format(key))
            n_valid_fields = 0
            if 'path' in file and file['path'] is not None:
                n_valid_fields += 1
            if 'shock_id' in file and file['shock_id'] is not None:
                n_valid_fields += 1
            if 'ftp_url' in file and file['ftp_url'] is not None:
                n_valid_fields += 1
                raise ValueError('FTP link is currently not supported for FastaGFFToGenome')
            if n_valid_fields < 1:
                error_msg = 'Required "{}" field must include one source: '.format(key)
                error_msg += 'path | shock_id'
                raise ValueError(error_msg)
            if n_valid_fields > 1:
                error_msg = 'Required "{}" field has too many sources specified: '.format(key)
                error_msg += str(file.keys())
                raise ValueError(error_msg)

        # check for valid type param
        valid_types = ['Reference', 'User upload', 'Representative']
        if params.get('type') and params['type'] not in valid_types:
            error_msg = 'Entered value for type is not one of the valid entries of '
            error_msg += '[' + ''.join('"' + str(e) + '", ' for e in valid_types)[0: -2] + ']'
            raise ValueError(error_msg)

    def _set_parsed_params(self, params):
        log('Setting params')

        # default params
        default_params = {
            'taxon_wsname': self.cfg.raw['taxon-workspace-name'],
            'scientific_name': 'unknown_taxon',
            'taxon_reference': None,
            'source': 'User',
            'release': None,
            'type': 'User upload',
            'metadata': {}
        }

        for field in default_params:
            if field not in params:
                params[field] = default_params[field]

        log(json.dumps(params, indent=1))

        return params

    def _stage_input(self, params, input_directory):
        """
        stage_input: Setup the input_directory by fetching the files and uncompressing if needed

        """

        file_paths = dict()
        for key in ('fasta_file', 'gff_file'):
            file = params[key]
            file_path = None
            if 'path' in file and file['path'] is not None:
                local_file_path = file['path']
                file_path = os.path.join(input_directory, os.path.basename(local_file_path))
                log('Moving file from {} to {}'.format(local_file_path, file_path))
                shutil.copy2(local_file_path, file_path)

            if 'shock_id' in file and file['shock_id'] is not None:
                # handle shock file
                log('Downloading file from SHOCK node: {}-{}'.format(
                                                        self.cfg.sharedFolder, file['shock_id']))
                sys.stdout.flush()
                file_name = self.dfu.shock_to_file({'file_path': input_directory,
                                                    'shock_id': file['shock_id']
                                                    })['node_file_name']
                file_path = os.path.join(input_directory, file_name)

            # extract the file if it is compressed
            if file_path is not None:
                print("staged input file =" + file_path)
                sys.stdout.flush()
                dfUtil_result = self.dfu.unpack_file({'file_path': file_path})
                file_paths[key] = dfUtil_result['file_path']
            else:
                raise ValueError('No valid files could be extracted based on the input')

        return file_paths

    @staticmethod
    def _retrieve_fasta_file(input_fasta_file, core_genome_name,
                             scientific_name, source):
        """
        _retrieve_fasta_file: retrieve info from fasta_file
                              https://www.biostars.org/p/710/

        """
        #TODO: Use assembly utils for this?
        log("Reading FASTA file")

        assembly = {"contigs": {}, "dna_size": 0, "gc_content": 0, "md5": [], "base_counts": {}}
        contig_seq_start = 0

        input_file_handle = open(input_fasta_file, 'rb')

        # alternate header and sequence
        faiter = (x[1] for x in itertools.groupby(input_file_handle, lambda line: line[0] == ">"))
        for header in faiter:
            # drop the ">"
            header = header.next()[1:].strip()
            # join all sequence lines to one.
            seq = "".join(s.strip() for s in faiter.next())

            try:
                fasta_header, fasta_description = header.split(' ', 1)
            except:
                fasta_header = header
                fasta_description = None

            # Handle record
            seq = seq.upper()

            # Build contig objects for Assembly
            seq_count = dict(collections.Counter(seq))

            # to delete at end, but required for now
            contig_dict = {"sequence": seq}

            Ncount = 0
            if "N" in seq_count:
                Ncount = seq_count["N"]
            contig_dict["Ncount"] = Ncount

            for character in seq_count:
                if character in assembly["base_counts"]:
                    assembly["base_counts"][character] += seq_count[character]
                else:
                    assembly["base_counts"][character] = seq_count[character]

            contig_seq_length = len(seq)
            assembly["dna_size"] += contig_seq_length
            contig_gc_length = seq.count("G")
            contig_gc_length += seq.count("C")
            contig_dict["gc_content"] = float("{0:.2f}".format(float(contig_gc_length) /
                                              float(contig_seq_length)))
            assembly["gc_content"] += contig_gc_length
            contig_dict["contig_id"] = fasta_header
            contig_dict["name"] = fasta_header
            contig_dict["length"] = contig_seq_length
            contig_dict["md5"] = hashlib.md5(seq).hexdigest()
            assembly["md5"].append(contig_dict["md5"])

            if fasta_description is not None:
                contig_dict["description"] = fasta_description

            contig_dict["is_circular"] = "Unknown"
            contig_dict["start_position"] = contig_seq_start
            contig_dict["num_bytes"] = sys.getsizeof(contig_dict["sequence"])
            assembly["contigs"][fasta_header] = contig_dict

            # used for start of next sequence and total gc_content
            contig_seq_start += contig_seq_length

        assembly["gc_content"] = float("{0:.2f}".format(float(assembly["gc_content"]) /
                                       float(contig_seq_start)))
        assembly["md5"] = hashlib.md5(",".join(assembly["md5"])).hexdigest()
        assembly["assembly_id"] = core_genome_name+"_assembly"
        assembly["name"] = scientific_name
        assembly["external_source"] = source
        assembly["external_source_id"] = os.path.basename(input_fasta_file)
        assembly["external_source_origination_date"] = str(os.stat(input_fasta_file).st_ctime)
        assembly["num_contigs"] = len(assembly["contigs"].keys())
        assembly["type"] = "Unknown"
        assembly["notes"] = "Note MD5s are generated from uppercasing the sequences"

        return assembly

    def _retrieve_gff_file(self, input_gff_file):
        """
        _retrieve_gff_file: retrieve info from gff_file
    
        """
        log("Reading GFF file")
    
        feature_list = dict()
        is_phytozome = 0
        is_patric = 0

        gff_file_handle = open(input_gff_file, 'rb')
        current_line = gff_file_handle.readline()
        line_count = 0

        while (current_line != ''):
            current_line = current_line.strip()

            if(current_line.isspace() or current_line == "" or current_line.startswith("#")):
                pass
            else:
                #Split line
                (contig_id, source_id, feature_type, start, end,
                 score, strand, phase, attributes) = current_line.split('\t')

                #Checking to see if Phytozome
                if("phytozome" in source_id or "Phytozome" in source_id):
                    is_phytozome=1

                #Checking to see if Phytozome
                if("PATRIC" in source_id):
                    is_patric=1

                #PATRIC prepends their contig ids with some gibberish
                if(is_patric and "|" in contig_id):
                    contig_id = contig_id.split("|",1)[1]

                #Features grouped by contigs first
                if(contig_id not in feature_list):
                    feature_list[contig_id] = list()

                #Populating basic feature object
                ftr = {'contig': contig_id, 'source': source_id,
                       'type': feature_type, 'start': int(start),
                       'end': int(end), 'score': score, 'strand': strand,
                       'phase': phase, 'attributes': collections.defaultdict(list)}

                #Populating with attribute key-value pair
                #This is where the feature id is from
                for attribute in attributes.split(";"):
                    attribute = attribute.strip()

                    #Sometimes empty string
                    if not attribute:
                        continue

                    #Use of 1 to limit split as '=' character can also be made available later
                    #Sometimes lack of "=", assume spaces instead
                    if("=" in attribute):
                        key, value = attribute.split("=", 1)
                        ftr['attributes'][key].append(value)
                    elif(" " in attribute):
                        key, value = attribute.split(" ", 1)
                        ftr['attributes'][key].append(value)
                    else:
                        log("Warning: attribute "+attribute+" cannot be separated into key,value pair")
                if "ID" in ftr['attributes']:
                    ftr['ID'] = ftr['attributes']['ID'][0]
                if "Parent" in ftr['attributes']:
                    ftr['Parent'] = ftr['attributes']['Parent'][0]

                feature_list[contig_id].append(ftr)

            current_line = gff_file_handle.readline()

        gff_file_handle.close()

        #Some GFF/GTF files don't use "ID" so we go through the possibilities        
        feature_list = self._add_missing_identifiers(feature_list)

        #Most bacterial files have only CDSs
        #In order to work with prokaryotic and eukaryotic gene structure synonymously
        #Here we add feature dictionaries representing the parent gene and mRNAs
        feature_list = self._add_missing_parents(feature_list)

        #Phytozome has the annoying habit of editing their identifiers so we fix them
        if(is_phytozome):
            self._update_phytozome_features(feature_list)

        #All identifiers need to be checked so that they follow the same general rules
        #Rules are listed within the function itself
        feature_list = self._update_identifiers(feature_list)

        #If phytozome, the edited files need to be re-printed as GFF so that it works better with RNA-Seq pipeline
        if(is_phytozome):
            self._print_phytozome_gff(input_gff_file,feature_list)

        return feature_list

    def _add_missing_identifiers(self,feature_list):

        #General rule is to iterate through a range of possibilities if "ID" is missing
        for contig in feature_list.keys():
            for i in range(len(feature_list[contig])):
                if("ID" not in feature_list[contig][i]):
                    for key in ("transcriptId", "proteinId", "PACid", "pacid", "Parent"):
                        if(key in feature_list[contig][i]):
                            feature_list[contig][i]['ID']=feature_list[contig][i][key]
                            break

                    #If the process fails, throw an error
                    for ftr_type in ("gene", "mRNA", "CDS"):
                        if(ftr_type not in feature_list[contig][i]):
                            continue

                        if("ID" not in feature_list[contig][i]):
                            log("Error: Cannot find unique ID to utilize in GFF attributes: "+ \
                                    feature_list[contig][i]['contig']+"."+ \
                                    feature_list[contig][i]['source']+"."+ \
                                    feature_list[contig][i]['type']+": "+ \
                                    feature_list[contig][i]['attributes'])
        return feature_list

    @staticmethod
    def _generate_feature_hierarchy(feature_list):

        feature_hierarchy = { contig : {} for contig in feature_list }

        #Need to remember mRNA/gene links for CDSs
        mRNA_gene_dict = {}
        exon_list_position_dict = {}

        for contig in feature_list:
            for i in range(len(feature_list[contig])):
                ftr = feature_list[contig][i]
            
                if("gene" in ftr["type"]):
                    feature_hierarchy[contig][ftr["ID"]] = {
                        "utrs": [], "mrnas": [], "cdss": [], "index": i}

                if("UTR" in ftr["type"]):
                    feature_hierarchy[contig][mRNA_gene_dict[ftr["Parent"]]]["utrs"].append(
                        {"id": ftr["ID"], "index": i})

                if("RNA" in ftr["type"]):
                    feature_hierarchy[contig][ftr["Parent"]]["mrnas"].append(
                        {"id": ftr["ID"], "index": i, "cdss": []})
                    mRNA_gene_dict[ftr["ID"]] = ftr["Parent"]
                    exon_list_position_dict[ftr["ID"]] = len(
                        feature_hierarchy[contig][ftr["Parent"]]["mrnas"])-1

                if("CDS" in ftr["type"]):
                    feature_hierarchy[contig][mRNA_gene_dict[ftr["Parent"]]]["mrnas"]\
                        [exon_list_position_dict[ftr["Parent"]]]["cdss"].append(
                        {"id": ftr["ID"], "index": i})

                #TODO: add other feature types
                                                                                      
        return feature_hierarchy

    @staticmethod
    def _add_missing_parents(feature_list):

        #General rules is if CDS or RNA missing parent, add them
        for contig in feature_list.keys():
            ftrs = feature_list[contig]
            new_ftrs = []
            for i in range(len(ftrs)):
                if("Parent" not in ftrs[i]):
                    #Assuming parent doesn't exist at all, so create de novo instead of trying to find it
                    if("RNA" in ftrs[i]["type"] or "CDS" in ftrs[i]["type"]):
                        new_gene_ftr = copy.deepcopy(ftrs[i])
                        new_gene_ftr["type"] = "gene"
                        ftrs[i]["Parent"]=new_gene_ftr["ID"]
                        new_ftrs.append(new_gene_ftr)

                    if("CDS" in ftrs[i]["type"]):
                        new_rna_ftr = copy.deepcopy(ftrs[i])
                        new_rna_ftr["type"] = "mRNA"
                        new_ftrs.append(new_rna_ftr)
                        ftrs[i]["Parent"]=new_rna_ftr["ID"]

                new_ftrs.append(ftrs[i])
            feature_list[contig]=new_ftrs
        return feature_list

    @staticmethod
    def _update_phytozome_features(feature_list):

        #General rule is to use the "Name" field where possible
        #And update parent attribute correspondingly
        for contig in feature_list.keys():
            feature_position_dict = {}
            for i in range(len(feature_list[contig])):

                #Maintain old_id for reference
                #Sometimes ID isn't available, so use PACid
                old_id = None
                for key in ("ID", "PACid", "pacid"):
                    if(key in feature_list[contig][i]):
                        old_id = feature_list[contig][i][key]
                        break
                if(old_id is None):
                    #This should be an error
                    print ("Cannot find unique ID, PACid, or pacid in GFF "
                           "attributes: " + feature_list[contig][i][contig])
                    continue

                #Retain old_id
                feature_position_dict[old_id]=i

                #In Phytozome, gene and mRNA have "Name" field, CDS do not
                if("Name" in feature_list[contig][i]):
                    feature_list[contig][i]["ID"] = feature_list[contig][i]["Name"]

                if("Parent" in feature_list[contig][i]):
                    #Update Parent to match new ID of parent ftr
                    feature_list[contig][i]["Parent"] = feature_list[contig][feature_position_dict[feature_list[contig][i]["Parent"]]]["ID"]

        return feature_list

    @staticmethod
    def _update_identifiers(feature_list):

        #General rules:
        #1) Genes keep identifier
        #2) RNAs keep identifier only if its different from gene, otherwise append ".mRNA"
        #3) CDS always uses RNA identifier with ".CDS" appended
        #4) CDS appended with an incremented digit

        CDS_count_dict = dict()
        mRNA_parent_dict = dict()

        for contig in feature_list.keys():
            for ftr in feature_list[contig]:
                if("Parent" in ftr):

                    #Retain old_id of parents
                    old_id = ftr["ID"]

                    if(ftr["ID"] == ftr["Parent"] or "CDS" in ftr["type"]):
                        ftr["ID"] = ftr["Parent"]+"."+ftr["type"]

                    #link old to new ids for mRNA to use with CDS
                    if("RNA" in ftr["type"]):
                        mRNA_parent_dict[old_id]=ftr["ID"]

                    if("CDS" in ftr["type"]):
                        #Increment CDS identifier
                        if(ftr["ID"] not in CDS_count_dict):
                            CDS_count_dict[ftr["ID"]] = 1
                        else:
                            CDS_count_dict[ftr["ID"]] += 1
                        ftr["ID"]=ftr["ID"]+"."+str(CDS_count_dict[ftr["ID"]])

                        #Recall new mRNA id for parent
                        ftr["Parent"]=mRNA_parent_dict[ftr["Parent"]]

        return feature_list

    @staticmethod
    def _print_phytozome_gff(input_gff_file, feature_list):

        #Write modified feature ids to new file
        input_gff_file = input_gff_file.replace("gene", "edited_gene")+".gz"
        try:
            print "Printing to new file: "+input_gff_file
            gff_file_handle = gzip.open(input_gff_file, 'wb')
        except:
            print "Failed to open"

        for contig in sorted(feature_list.iterkeys()):
            for ftr in feature_list[contig]:

                #Re-build attributes
                attributes_dict = {}
                for attribute in ftr["attributes"].split(";"):
                    attribute=attribute.strip()

                    #Sometimes empty string
                    if(attribute == ""):
                        continue

                    #Use of 1 to limit split as '=' character can also be made available later
                    #Sometimes lack of "=", assume spaces instead
                    if("=" in attribute):
                        key, value = attribute.split("=", 1)
                    elif(" " in attribute):
                        key, value = attribute.split(" ", 1)
                    else:
                        log("Warning: attribute "+attribute+" cannot be separated into key,value pair")

                    if(ftr[key] != value):
                        value = ftr[key]
                    attributes_dict[key]=value

                ftr["attributes"]=";".join(key+"="+attributes_dict[key] for key in attributes_dict.keys())

                new_line = "\t".join( str(ftr[key]) for key in ['contig', 'source', 'type', 'start', 'end',
                                                                'score', 'strand', 'phase', 'attributes'])
                gff_file_handle.write(new_line)
        gff_file_handle.close()
        return

    def _retrieve_genome_feature_list(self, feature_list, feature_hierarchy, assembly):

        genome_features_list = list()
        genome_mrnas_list = list()
        genome_cdss_list = list()
        genome_translation_issues = list()

        for contig in feature_hierarchy:
            for gene in feature_hierarchy[contig]:

                #We only iterate through the gene objects
                #And then for each gene object, retrieve the necessary mRNA and CDS objects indirectly

                ftr = feature_list[contig][feature_hierarchy[contig][gene]["index"]]
                contig_sequence = assembly["contigs"][ftr["contig"]]["sequence"]
                gene_ftr = self._convert_ftr_object(ftr, contig_sequence) #reverse-complementation for negative strands done here

                #Add non-optional terms
                gene_ftr["mrnas"] = list()
                gene_ftr["cdss"] = list()
                gene_ftr["ontology_terms"] = dict()

                #Retaining longest sequences for gene feature
                longest_protein_length = 0
                longest_protein_sequence = ""
                for mRNA in feature_hierarchy[contig][gene]["mrnas"]:

                    ########################################################
                    # Construct mRNA Ftr
                    ########################################################
                    ftr = feature_list[contig][mRNA["index"]]
                    contig_sequence = assembly["contigs"][ftr["contig"]]["sequence"]
                    mRNA_ftr = self._convert_ftr_object(ftr, contig_sequence) #reverse-complementation for negative strands done here

                    #Modify mrna object for use in mrna array
                    #Objects will be un-used until further notice
                    mRNA_ftr['parent_gene'] = gene_ftr['id']

                    #If there are CDS, then New CDS ID without incrementation as they were aggregated
                    if(len(mRNA['cdss'])>0):
                        mRNA_ftr['cds'] = mRNA_ftr['id']+".CDS"
                    else:
                        mRNA_ftr['cds'] = ""

                    #Add to mrnas array
                    genome_mrnas_list.append(mRNA_ftr)

                    #Add ids to gene_ftr arrays
                    gene_ftr["mrnas"].append(mRNA_ftr["id"])

                    ########################################################
                    # Construct transcript, protein sequence, UTR, CDS locations
                    ########################################################

                    #At time of writing, all of this aggregation should probably be done in a single function
                    cds_exons_locations_array = list()
                    cds_cdna_sequence = str()
                    protein_sequence = str()
                    if(len(mRNA["cdss"])>0):
                        (cds_exons_locations_array, cds_cdna_sequence, protein_sequence) = \
                            self._cds_aggregation_translation(mRNA["cdss"],feature_list[contig],assembly,genome_translation_issues)
                    
                    UTRs=list()
                    if("utrs" in feature_hierarchy[contig][gene] and len(feature_hierarchy[contig][gene]["utrs"])>0):
                        for UTR in feature_hierarchy[contig][gene]["utrs"]:
                            ftr = feature_list[contig][UTR["index"]]
                            if("Parent" in ftr and ftr["Parent"] == mRNA_ftr["id"]):
                                UTRs.append(ftr)

                    mrna_exons_locations_array = copy.deepcopy(cds_exons_locations_array)
                    mrna_transcript_sequence = str(cds_cdna_sequence)
                    if(len(UTRs)>0):
                        (mrna_exons_locations_array, mrna_transcript_sequence) = \
                            self._utr_aggregation(UTRs,assembly,mrna_exons_locations_array,cds_cdna_sequence)

                    #Update sequence and locations
                    mRNA_ftr["dna_sequence"]=mrna_transcript_sequence
                    mRNA_ftr["dna_sequence_length"]=len(mrna_transcript_sequence)
                    mRNA_ftr["location"]=mrna_exons_locations_array
                    mRNA_ftr["md5"] = hashlib.md5(mRNA_ftr["dna_sequence"]).hexdigest()

                    #Remove DNA
                    del mRNA_ftr["dna_sequence"]
                    del mRNA_ftr["dna_sequence_length"]

                    #Skip CDS if not present
                    if(len(mRNA["cdss"])==0):
                        continue

                    #Remove asterix representing stop codon if present
                    if(len(protein_sequence)>0 and protein_sequence[-1] == '*'):
                        protein_sequence = protein_sequence[:-1]

                    #Save longest sequence
                    if(len(protein_sequence) > longest_protein_length):
                        longest_protein_length = len(protein_sequence)
                        longest_protein_sequence = protein_sequence

                    ########################################################
                    # Construct CDS Ftr
                    ########################################################
                    CDS_ftr = dict()
                    CDS_ftr['type']='CDS'

                    #New CDS ID without incrementation as they were aggregated
                    CDS_ftr['id'] = mRNA_ftr['id']+'.CDS'

                    #Add gene/mrna links
                    CDS_ftr['parent_gene']=gene_ftr['id']
                    CDS_ftr['parent_mrna']=mRNA_ftr['id']

                    #Update sequence and locations
                    CDS_ftr["dna_sequence"]=cds_cdna_sequence
                    CDS_ftr["dna_sequence_length"]=len(cds_cdna_sequence)
                    CDS_ftr["location"]=cds_exons_locations_array
                    CDS_ftr["md5"] = hashlib.md5(CDS_ftr["dna_sequence"]).hexdigest()

                    #Add protein
                    CDS_ftr["protein_translation"] = str(protein_sequence).upper()
                    CDS_ftr["protein_translation_length"] = len(CDS_ftr["protein_translation"])
                    CDS_ftr["protein_md5"] = hashlib.md5(CDS_ftr["protein_translation"]).hexdigest()

                    #Add empty non-optional fields for populating in future
                    CDS_ftr["ontology_terms"] = dict()
                    if "aliases" not in CDS_ftr:
                        CDS_ftr["aliases"] = list()
                    if "function" not in CDS_ftr:
                        CDS_ftr["function"] = []
                    if 'product' in CDS_ftr:
                        CDS_ftr["function"].append(CDS_ftr['product'])

                    #Add to cdss array
                    genome_cdss_list.append(CDS_ftr)

                    #Add ids to gene_ftr arrays
                    gene_ftr["cdss"].append(CDS_ftr["id"])

                gene_ftr["protein_translation"] = longest_protein_sequence
                gene_ftr["protein_translation_length"] = longest_protein_length
                genome_features_list.append(gene_ftr)

        msg = "Genome features processed: {} genes, {} RNAs, and {} CDSs\n".format(len(genome_features_list),len(genome_mrnas_list),len(genome_cdss_list))
        msg += "{} mRNA(s) had errors during translation".format(len(genome_translation_issues))
        log(msg)

        return genome_features_list, genome_mrnas_list, genome_cdss_list

    def _transform_feature(self, contig, in_feature):
        """Converts a feature from the gff ftr format into the appropriate
        format for a genome object """
        def _aliases(feat):
            keys = ('locus_tag', 'old_locus_tag', 'protein_id',
                    'transcript_id', 'gene', 'EC_number')
            alias_list = []
            for key in keys:
                if key in feat['attributes']:
                    alias_list.extend([(key, val) for val in feat['attributes'][key]])
            return alias_list

        def _propagate_cds_props_to_gene(cds, gene):
            # Put longest protein_translation to gene
            if "protein_translation" not in gene or (
                len(gene["protein_translation"]) <
                    len(cds["protein_translation"])):
                gene["protein_translation"] = cds["protein_translation"]
                gene["protein_translation_length"] = len(
                    cds["protein_translation"])
            # Merge cds list attributes with gene
            for key in ('function', 'aliases', 'db_xref'):
                if cds.get(key, []):
                    gene[key] = cds.get(key, []) + gene.get(key, [])
            # Merge cds["ontology_terms"] -> gene["ontology_terms"]
            terms2 = cds.get("ontology_terms")
            if terms2 is not None:
                terms = gene.get("ontology_terms")
                if terms is None:
                    gene["ontology_terms"] = terms2
                else:
                    for source in terms2:
                        if source in terms:
                            terms[source].update(terms2[source])
                        else:
                            terms[source] = terms2[source]

        # UTRs don't have any information not implied by the exon
        excluded_features = ('five_prime_UTR', 'three_prime_UTR')
        if in_feature['type'] in excluded_features:
            self.skiped_features[in_feature['type']] += 1
            return

        # if the feature ID is duplicated we only need to update the location
        if in_feature['ID'] in self.feature_dict:
            self.feature_dict[in_feature['ID']]['location'].append(
                self._location(in_feature))

        self.feature_counts[in_feature['type']] += 1
        feat_seq = contig.seq[in_feature['start']-1:in_feature['end']-1]
        if in_feature['strand'] == '-':
            feat_seq = feat_seq.reverse_complement()

        # The following is common to all the feature types
        out_feat = {
            "id": in_feature['ID'],
            "type": in_feature['type'],
            "location": [self._location(in_feature)],
            "dna_sequence": str(feat_seq),
            "dna_sequence_length": len(feat_seq),
            "md5": hashlib.md5(str(feat_seq)).hexdigest(),
            "function": in_feature['attributes'].get('function', []),
            'warnings': []
        }
        # TODO: Flags, inference_data
        # add optional fields
        if 'note' in in_feature['attributes']:
            out_feat['note'] = in_feature['attributes']["note"][0]
        ont = self._get_ontology(in_feature)
        if ont:
            out_feat['ontology_terms'] = ont
        aliases = _aliases(in_feature)
        if aliases:
            out_feat['aliases'] = aliases
        if 'db_xref' in in_feature['attributes']:
            out_feat['db_xref'] = [tuple(x.split(":")) for x in
                                   in_feature['attributes']['db_xref']]
        if 'product' in in_feature['attributes']:
            out_feat['function'] += ["product:" + x for x in
                                     in_feature['attributes']["product"]]
        parent_id = in_feature.get('Parent', '')
        if parent_id and parent_id not in self.feature_dict:
            raise ValueError("Parent ID: {} was not found in feature ID list.")

        # add type specific features
        if in_feature['type'] == 'gene':
            out_feat['mrnas'] = []
            out_feat['cdss'] = []

        elif in_feature['type'] == 'CDS':
            prot_seq = str(feat_seq.translate()).strip("*")
            out_feat.update({
                "protein_translation": prot_seq,
                "protein_md5": hashlib.md5(prot_seq).hexdigest(),
                "protein_translation_length": len(prot_seq),
                'parent_gene': "",
            })
            if parent_id:
                # TODO: add location checks and warnings
                parent = self.feature_dict[parent_id]
                if 'cdss' in parent:  # parent must be a gene
                    parent['cdss'].append(in_feature['ID'])
                    out_feat['parent_gene'] = parent_id
                else:  # parent must be mRNA
                    parent['cds'] = in_feature['ID']
                    out_feat['parent_mrna'] = parent_id
                    parent_gene = self.feature_dict[parent['parent_gene']]
                    parent_gene['cdss'].append(in_feature['ID'])
                    out_feat['parent_gene'] = parent['parent_gene']
                    _propagate_cds_props_to_gene(out_feat, parent_gene)

        elif in_feature['type'] == 'mRNA':
            if parent_id:
                # TODO: add location checks and warnings
                parent = self.feature_dict[parent_id]
                if 'cdss' in parent:  # parent must be a gene
                    parent['mrnas'].append(in_feature['ID'])
                    out_feat['parent_gene'] = parent_id

        else:
            out_feat["type"] = in_feature.type
            if parent_id:
                # TODO: add location checks and warnings
                parent = self.feature_dict[parent_id]
                if 'children' not in parent:
                    parent['children'] = []
                parent['children'].append(out_feat['id'])
                out_feat['parent_gene'] = parent_id

        self.feature_dict[out_feat['id']] = out_feat

    def _gen_genome_info(self, core_genome_name, scientific_name, assembly_ref,
                         source, assembly, input_gff_file, molecule_type):
        """
        _gen_genome_info: generate genome info

        """
        genome = dict()
        genome["id"] = core_genome_name
        genome["scientific_name"] = scientific_name
        genome["assembly_ref"] = assembly_ref
        genome['molecule_type'] = molecule_type
        genome["features"] = []
        genome["cdss"] = []
        genome["mrnas"] = []
        genome['non_coding_features'] = []
        genome["gc_content"] = assembly["gc_content"]
        genome["dna_size"] = assembly["dna_size"]
        genome['md5'] = assembly['md5']
        genome['contig_ids'], genome['contig_lengths'] = zip(
            *[(k, v['length']) for k, v in assembly['contigs'].items()])
        genome['num_contigs'] = len(assembly['contigs'])
        genome["ontology_events"] = [{
            "method": "GenomeFileUtils Genbank uploader from annotations",
            "method_version": self.version,
            "timestamp": self.time_string,
            # TODO: remove this hardcoding
            "id": "GO",
            "ontology_ref": "KBaseOntology/gene_ontology"
        }]
        genome['ontologies_present'] = dict(self.ontologies_present)
        genome['feature_counts'] = dict(self.feature_counts)

        genome['taxonomy'], genome['taxon_ref'], genome['domain'], \
            genome["genetic_code"] = self.gi.retrieve_taxon(self.taxon_wsname,
                                                            genome['scientific_name'])
        genome['source'], genome['genome_tiers'] = self.gi.determine_tier(
            source)

        #gff_file_to_shock = self.dfu.file_to_shock(
        #    {'file_path': input_gff_file, 'make_handle': 1, 'pack': "gzip"})
        #genome['gff_handle_ref'] = gff_file_to_shock['handle']['hid']

        # sort features into their respective arrays
        for feature in self.feature_dict.values():
            if 'exons' in feature:
                del feature['exons']
            if feature['type'] == 'CDS':
                del feature['type']
                genome['cdss'].append(feature)
            elif feature['type'] == 'mRNA':
                del feature['type']
                genome['mrnas'].append(feature)
            elif feature['type'] == 'gene':
                del feature['type']
                if genome['cdss']:
                    self.feature_counts["protein_encoding_gene"] += 1
                    genome['features'].append(feature)
                else:
                    del genome['mrnas'], genome['cdss']
                    self.feature_counts["non-protein_encoding_gene"] += 1
                    genome['non_coding_features'].append(feature)

        return genome

    def _convert_ftr_object(self, old_ftr, contig):
        new_ftr = dict()
        new_ftr["id"] = old_ftr["ID"]

        dna_sequence = Seq(contig[old_ftr["start"]-1:old_ftr["end"]], IUPAC.ambiguous_dna)

        # reverse complement
        if(old_ftr["strand"] == "-"):
            dna_sequence = dna_sequence.reverse_complement()
            old_start = old_ftr["start"]
            old_ftr["start"] = old_ftr["end"]
            old_ftr["end"]=old_start

        new_ftr["dna_sequence"] = str(dna_sequence).upper()
        new_ftr["dna_sequence_length"] = len(dna_sequence)
        new_ftr["md5"] = hashlib.md5(str(dna_sequence)).hexdigest()
        new_ftr["location"] = [[old_ftr["contig"], old_ftr["start"], 
                                old_ftr["strand"], len(dna_sequence)]]
        new_ftr["type"] = old_ftr["type"]

        new_ftr["aliases"] = list()
        for key in ("transcriptId", "proteinId", "PACid", "pacid"):
            if key in old_ftr.keys():
                new_ftr["aliases"].append(key+":"+old_ftr[key])

        return new_ftr

    def _utr_aggregation(self, utr_list, assembly, exons, exon_sequence):

        #create copies of locations and transcript
        utrs_exons = list(exons)
        utr_exon_sequence = exon_sequence

        five_prime_dna_sequence = ""
        three_prime_dna_sequence = ""
        five_prime_locations = list()
        three_prime_locations = list()

        for UTR in (utr_list):
            contig_sequence = assembly["contigs"][UTR["contig"]]["sequence"]
            UTR_ftr = self._convert_ftr_object(UTR, contig_sequence)  #reverse-complementation for negative strands done here

            #aggregate sequences and locations
            if("five_prime" in UTR_ftr["id"]):
                five_prime_dna_sequence += UTR_ftr["dna_sequence"]
                five_prime_locations.append(UTR_ftr["location"][0])
            if("three_prime" in UTR_ftr["id"]):
                three_prime_dna_sequence += UTR_ftr["dna_sequence"]
                three_prime_locations.append(UTR_ftr["location"][0])

        #Handle five_prime UTRs
        if(len(five_prime_locations)>0):

            #Sort UTRs by "start" (reverse-complement UTRs in Phytozome appear to be incorrectly ordered in the GFF file
            five_prime_locations = sorted(five_prime_locations, key=lambda x: x[1])

            #Merge last UTR with CDS if "next" to each other
            if(five_prime_locations[-1][1]+five_prime_locations[-1][3] == utrs_exons[0][1]):

                #Remove last UTR
                last_five_prime_location = five_prime_locations[-1]
                five_prime_locations = five_prime_locations[:-1]

                #"Add" last UTR to first exon
                utrs_exons[0][1]=last_five_prime_location[1]
                utrs_exons[0][3]+=last_five_prime_location[3]
                        
            #Prepend other UTRs if available
            if(len(five_prime_locations)>0):
                utrs_exons = five_prime_locations + utrs_exons

        utr_exon_sequence = five_prime_dna_sequence+utr_exon_sequence

        #Handle three_prime UTRs
        if(len(three_prime_locations)>0):

            #Sort UTRs by "start" (reverse-complement UTRs in Phytozome appear to be incorrectly ordered in the GFF file
            three_prime_locations = sorted(three_prime_locations, key=lambda x: x[1])

            #Merge first UTR with CDS if "next to each other
            if(utrs_exons[-1][1]+utrs_exons[-1][3] == three_prime_locations[0][1]):

                #Remove first UTR
                first_three_prime_location = three_prime_locations[0]
                three_prime_locations = three_prime_locations[1:]

                #"Add" first UTR to last exon
                utrs_exons[-1][3]+=first_three_prime_location[3]

        #Append other UTRs if available
        if(len(three_prime_locations)>0):
            utrs_exons = utrs_exons + three_prime_locations

        utr_exon_sequence += three_prime_dna_sequence

        return (utrs_exons, utr_exon_sequence)

    def _cds_aggregation_translation(self, cds_list, feature_list, assembly, issues):

        dna_sequence = ""
        locations = list()

        # collect phases, and lengths of exons
        # right now, this is only for the purpose of error reporting
        phases = list()
        exons = list()

        #Saving parent mRNA identifier
        Parent_mRNA = cds_list[0]["id"]
        for CDS in (cds_list):
            ftr = feature_list[CDS["index"]]
            phases.append(ftr["phase"])
            Parent_mRNA=ftr["Parent"]

            contig_sequence = assembly["contigs"][ftr["contig"]]["sequence"]
            CDS_ftr = self._convert_ftr_object(ftr, contig_sequence) #reverse-complementation for negative strands done here
            exons.append(len(CDS_ftr["dna_sequence"]))

            # Remove base(s) according to phase, but only for first CDS
            if(CDS == cds_list[0] and int(ftr["phase"]) != 0):
                log("Adjusting phase for first CDS: "+CDS["id"])
                CDS_ftr["dna_sequence"] = CDS_ftr["dna_sequence"][int(ftr["phase"]):]

            #aggregate sequences and locations
            dna_sequence += CDS_ftr["dna_sequence"]
            locations.append(CDS_ftr["location"][0])

        # translate sequence
        dna_sequence_obj = Seq(dna_sequence, IUPAC.ambiguous_dna)
        rna_sequence = dna_sequence_obj.transcribe()

        # incomplete gene model with no start codon
        if str(rna_sequence.upper())[:3] not in codon_table.start_codons:
            msg = "Missing start codon for {}. Possibly incomplete gene model.".format(Parent_mRNA)
            log(msg)

        # You should never have this problem, needs to be reported rather than "fixed"
        codon_count = len(str(rna_sequence)) % 3
        if codon_count != 0:
            msg = "Number of bases for RNA sequence for {} ".format(Parent_mRNA)
            msg += "is not divisible by 3. "
            msg += "The resulting protein may well be mis-translated."
            log(msg)
            issues.append(Parent_mRNA)

        protein_sequence = Seq("")
        try:
            protein_sequence = rna_sequence.translate()
        except CodonTable.TranslationError as te:
            log("TranslationError")

        return (locations,dna_sequence.upper(),str(protein_sequence).upper())
