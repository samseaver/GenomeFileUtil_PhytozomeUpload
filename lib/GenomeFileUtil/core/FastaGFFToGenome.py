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

# KBase imports
from DataFileUtil.DataFileUtilClient import DataFileUtil
from AssemblyUtil.AssemblyUtilClient import AssemblyUtil
from KBaseReport.KBaseReportClient import KBaseReport

# 3rd party imports
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable

codon_table = CodonTable.ambiguous_generic_by_name["Standard"]


def log(message, prefix_newline=False):
    """Logging function, provides a hook to suppress or redirect log messages."""
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))


class FastaGFFToGenome:

    def __init__(self, config):
        self.cfg = config
        self.dfu = DataFileUtil(self.cfg.callbackURL)

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
            shock_service_url=self.cfg.shockURL,
            handle_service_url=self.cfg.handleURL,
            workspace_service_url=self.cfg.workspaceURL,
            callback_url=self.cfg.callbackURL,

            input_fasta_file=file_paths["fasta_file"],
            input_gff_file=file_paths["gff_file"],

            workspace_name=params['workspace_name'],
            core_genome_name=params['genome_name'],
            scientific_name=params['scientific_name'],
            taxon_wsname=params['taxon_wsname'],
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

    def upload_genome(self, shock_service_url=None, handle_service_url=None,
                      workspace_service_url=None, callback_url=None, input_gff_file=None,
                      input_fasta_file=None, workspace_name=None, core_genome_name=None,
                      scientific_name="unknown_taxon", taxon_wsname='ReferenceTaxons',
                      taxon_reference=None, source=None, release=None, genome_type=None):

        # retrieve taxon
        taxonomy, taxon_reference = self._retrieve_taxon(taxon_reference,
                                                         taxon_wsname, scientific_name)

        # reading in Fasta file
        assembly = self._retrieve_fasta_file(input_fasta_file, core_genome_name,
                                             scientific_name, source)

        if taxon_reference is not None:
            assembly["taxon_ref"] = taxon_reference

        # reading in GFF file
        feature_list = self._retrieve_gff_file(input_gff_file)

        # compile links between features
        feature_hierarchy = self._generate_feature_hierarchy(feature_list)

        # retrieve genome feature list
        (genome_features_list,
         genome_mrnas_list,
         genome_cdss_list) = self._retrieve_genome_feature_list(feature_list, feature_hierarchy, assembly)

        # remove sequences before loading
        for contig in assembly["contigs"]:
            del assembly["contigs"][contig]["sequence"]

        aUtil = AssemblyUtil(callback_url)
        assembly_ref = aUtil.save_assembly_from_fasta(
                                                {'file':
                                                    {'path': input_fasta_file,
                                                     'assembly_name': assembly['assembly_id']},
                                                 'workspace_name': workspace_name,
                                                 'assembly_name': assembly['assembly_id']})

        # generate genome info
        genome = self._gen_genome_info(core_genome_name, scientific_name, assembly_ref,
                                       genome_features_list, genome_cdss_list, genome_mrnas_list,
                                       source, assembly, taxon_reference, taxonomy, input_gff_file)

        workspace_id = self.dfu.ws_name_to_id(workspace_name)
        genome_info = self.dfu.save_objects({"id": workspace_id,
                                             "objects": [{"name": core_genome_name,
                                                          "type": "KBaseGenomes.Genome",
                                                          "data": genome}]})[0]
        report_string = ''

        return {'genome_info': genome_info, 'report_string': report_string}

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

    def _retrieve_taxon(self, taxon_reference, taxon_wsname, scientific_name):
        """
        _retrieve_taxon: retrieve taxonomy and taxon_reference

        """
        taxon_id = -1
        taxon_object_name = "unknown_taxon"

        # retrieve lookup object if scientific name provided
        if(taxon_reference is None and scientific_name is not "unknown_taxon"):
            # retrieve taxon lookup object then find taxon id
            taxon_lookup = self.dfu.get_objects(
                                    {'object_refs': [taxon_wsname+"/taxon_lookup"],
                                     'ignore_errors': 0})['data'][0]['data']['taxon_lookup']

            if(scientific_name[0:3] in taxon_lookup and
               scientific_name in taxon_lookup[scientific_name[0:3]]):
                taxon_id = taxon_lookup[scientific_name[0:3]][scientific_name]
                taxon_object_name = "{}_taxon".format(str(taxon_id))

        # retrieve Taxon object
        taxon_info = {}
        if(taxon_reference is None):
            taxon_info = self.dfu.get_objects({'object_refs': [taxon_wsname+"/"+taxon_object_name],
                                               'ignore_errors': 0})['data'][0]
            taxon_reference = "{}/{}/{}".format(taxon_info['info'][6],
                                                taxon_info['info'][0],
                                                taxon_info['info'][4])
        else:
            taxon_info = self.dfu.get_objects({"object_refs": [taxon_reference],
                                              'ignore_errors': 0})['data'][0]

        taxonomy = taxon_info['data']['scientific_lineage']

        return taxonomy, taxon_reference

    def _retrieve_fasta_file(self, input_fasta_file, core_genome_name,
                             scientific_name, source):
        """
        _retrieve_fasta_file: retrieve info from fasta_file
                              https://www.biostars.org/p/710/

        """
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

                #Only allow certain types for the time being
                if("RNA" not in feature_type and "CDS" not in feature_type and "gene" not in feature_type):
                    current_line = gff_file_handle.readline()
                    continue

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
                ftr = {'contig': contig_id, 'source': source_id, 'type': feature_type, 'start': int(start), 'end': int(end),
                       'score': score, 'strand': strand, 'phase': phase, 'attributes': attributes}

                #Populating with attribute key-value pair
                #This is where the feature id is from
                for attribute in attributes.split(";"):
                    #Sometimes empty string or lack of "="
                    if(attribute == "" or "=" not in attribute):
                        continue

                    #Use of 1 to limit split as '=' character can also be made available later
                    key, value = attribute.split("=", 1)
                    ftr[key] = value

                feature_list[contig_id].append(ftr)

            current_line = gff_file_handle.readline()

        gff_file_handle.close()

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

    def _generate_feature_hierarchy(self,feature_list):

        feature_hierarchy = { contig : {} for contig in feature_list }

        #Need to remember mRNA/gene links for CDSs
        mRNA_gene_dict = {}
        mRNA_position_dict = {}

        for contig in feature_list:
            for i in range(len(feature_list[contig])):
                ftr = feature_list[contig][i]
            
                if("gene" in ftr["type"]):
                    feature_hierarchy[contig][ftr["ID"]]={"mrnas":[],"cdss":[],"index":i}

                if("RNA" in ftr["type"]):
                    feature_hierarchy[contig][ftr["Parent"]]["mrnas"].append( { "id" : ftr["ID"], "index" : i, "cdss" : [] } )
                    mRNA_gene_dict[ftr["ID"]]=ftr["Parent"]
                    mRNA_position_dict[ftr["ID"]]=len(feature_hierarchy[contig][ftr["Parent"]]["mrnas"])-1
                if("CDS" in ftr["type"]):
                    feature_hierarchy[contig][mRNA_gene_dict[ftr["Parent"]]]["mrnas"][mRNA_position_dict[ftr["Parent"]]]["cdss"].append( { "id": ftr["ID"], "index" : i } )
                                                                                      
        return feature_hierarchy

    def _add_missing_parents(self,feature_list):

        #General rules is if CDS or RNA missing parent, add them
        for contig in feature_list.keys():
            ftrs = feature_list[contig]
            new_ftrs = []
            for i in range(len(ftrs)):
                if("Parent" not in ftrs[i]):
                    #Assuming parent doesn't exist at all, so create de novo instead of trying to find it
                    if("RNA" in ftrs[i]["type"] or "CDS" in ftrs[i]["type"]):
                        new_gene_ftr = copy.copy(ftrs[i])
                        new_gene_ftr["type"] = "gene"
                        ftrs[i]["Parent"]=new_gene_ftr["ID"]
                        new_ftrs.append(new_gene_ftr)

                    if("CDS" in ftrs[i]["type"]):
                        new_rna_ftr = copy.copy(ftrs[i])
                        new_rna_ftr["type"] = "mRNA"
                        new_ftrs.append(new_rna_ftr)
                        ftrs[i]["Parent"]=new_rna_ftr["ID"]

                new_ftrs.append(ftrs[i])
            feature_list[contig]=new_ftrs
        return feature_list

    def _update_phytozome_features(self,feature_list):

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
                    print ("Cannot find unique ID, PACid, or pacid in GFF attributes: "+feature_list[contig][i][contig],feature_list[contig][i][source],feature_list[contig][i][attributes])
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

    def _update_identifiers(self,feature_list):

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

    def _print_phytozome_gff(self, input_gff_file, feature_list):

        #Write modified feature ids to new file
        input_gff_file = input_gff_file.replace("gene", "edited_gene")+".gz"
        gff_file_handle = gzip.open(input_gff_file, 'wb')
        for contig in sorted(feature_list.iterkeys()):
            for ftr in feature_list[contig]:

                #Re-build attributes
                attributes_dict = {}
                for attribute in ftr["attributes"].split(";"):
                    #Sometimes empty string or lack of "="
                    if(attribute == "" or "=" not in attribute):
                        continue

                    #Use of 1 to limit split as '=' character can also be made available later
                    key, value = attribute.split("=", 1)
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
                gene_ftr = self._convert_ftr_object(ftr, contig_sequence)

                #Add non-optional terms
                gene_ftr["mrnas"] = list()
                gene_ftr["cdss"] = list()
                gene_ftr["ontology_terms"] = dict()

                #Retaining longest sequences for gene feature
                longest_protein_length = 0
                longest_protein_sequence = ""
                for mRNA in feature_hierarchy[contig][gene]["mrnas"]:
                    ftr = feature_list[contig][mRNA["index"]]
                    contig_sequence = assembly["contigs"][ftr["contig"]]["sequence"]
                    mRNA_ftr = self._convert_ftr_object(ftr, contig_sequence)

                    #Modify mrna object for use in mrna array
                    #Objects will be un-used until further notice
                    mRNA_ftr['parent_gene'] = gene_ftr['id']

                    #If there are CDS, then New CDS ID without incrementation as they were aggregated
                    if(len(mRNA['cdss'])>0):
                        mRNA_ftr['cds'] = mRNA_ftr['id']+".CDS"
                    else:
                        mRNA_ftr['cds'] = ""

                    #Remove DNA
                    del mRNA_ftr["dna_sequence"]
                    del mRNA_ftr["dna_sequence_length"]

                    #Add to mrnas array
                    genome_mrnas_list.append(mRNA_ftr)

                    #Add ids to gene_ftr arrays
                    gene_ftr["mrnas"].append(mRNA_ftr["id"])

                    #Checks to see if there's any CDS first, possible there isn't
                    if(len(mRNA['cdss'])==0):
                        continue

                    protein_sequence = self._cds_aggregation_translation(mRNA["cdss"],feature_list[contig],assembly,genome_translation_issues)
                    if(len(protein_sequence) > longest_protein_length):
                        longest_protein_length = len(protein_sequence)
                        longest_protein_sequence = protein_sequence

                    #Modify mrna object for use in cds array
                    #Objects will be un-used until further notice
                    CDS_ftr = copy.deepcopy(mRNA_ftr)

                    #New CDS ID without incrementation as they were aggregated
                    CDS_ftr['id'] = mRNA_ftr['id']+'.CDS'

                    #Remove CDS link and add mrna link
                    del CDS_ftr['cds']
                    CDS_ftr['parent_mrna']=mRNA_ftr['id']

                    #Add protein
                    CDS_ftr["protein_translation"] = str(protein_sequence).upper()
                    CDS_ftr["protein_translation_length"] = len(CDS_ftr["protein_translation"])
                    CDS_ftr["md5"] = hashlib.md5(CDS_ftr["protein_translation"]).hexdigest()

                    #Add empty non-optional fields for populating in future
                    CDS_ftr["ontology_terms"] = dict()
                    if("aliases" not in CDS_ftr):
                        CDS_ftr["aliases"] = list()
                    if("function" not in CDS_ftr):
                        CDS_ftr["function"] = ""

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

    def _gen_genome_info(self, core_genome_name, scientific_name, assembly_ref,
                         genome_features_list, genome_cdss_list, genome_mrnas_list,
                         source, assembly, taxon_reference, taxonomy, input_gff_file):
        """
        _gen_genome_info: generate genome info

        """
        genome = dict()
        genome["id"] = core_genome_name
        genome["scientific_name"] = scientific_name
        genome["assembly_ref"] = assembly_ref
        genome["features"] = genome_features_list
        genome["cdss"] = genome_cdss_list
        genome["mrnas"] = genome_mrnas_list
        genome["source"] = source
        genome["domain"] = "Eukaryota"
        genome["genetic_code"] = 1
        genome["gc_content"] = assembly["gc_content"]
        genome["dna_size"] = assembly["dna_size"]

        if taxon_reference is not None:
            genome["taxon_ref"] = taxon_reference
            genome["taxonomy"] = taxonomy

        gff_file_to_shock = self.dfu.file_to_shock({'file_path': input_gff_file,
                                                    'make_handle': 1, 'pack': "gzip"})
        gff_handle_ref = gff_file_to_shock['handle']['hid']

        genome['gff_handle_ref'] = gff_handle_ref

        return genome

    def _convert_ftr_object(self, old_ftr, contig):
        new_ftr = dict()
        new_ftr["id"] = old_ftr["ID"]

        dna_sequence = Seq(contig[old_ftr["start"]-1:old_ftr["end"]], IUPAC.ambiguous_dna)

        # reverse complement
        if(old_ftr["strand"] == "-"):
            dna_sequence = dna_sequence.reverse_complement()
            old_ftr["start"] = old_ftr["end"]

        new_ftr["dna_sequence"] = str(dna_sequence).upper()
        new_ftr["dna_sequence_length"] = len(dna_sequence)
        new_ftr["md5"] = hashlib.md5(str(dna_sequence)).hexdigest()
        new_ftr["location"] = [[old_ftr["contig"], old_ftr["start"], 
                                old_ftr["strand"], len(dna_sequence)]]
        new_ftr["type"]=old_ftr["type"]
        return new_ftr

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
            CDS_ftr = self._convert_ftr_object(ftr, contig_sequence)
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
            log("TranslationError for: "+feature_object["id"], phases, exons, " : "+str(te))

        return str(protein_sequence).upper()
