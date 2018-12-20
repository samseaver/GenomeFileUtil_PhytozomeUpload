import os
import sys
import shutil
import uuid
import time
import json
import hashlib
import collections
import datetime
import re
import copy
import urllib.parse as parse

# KBase imports
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.AssemblyUtilClient import AssemblyUtil
from GenomeFileUtil.core import GenomeUtils
from GenomeFileUtil.core.GenomeUtils import is_parent, warnings, check_full_contig_length_or_multi_strand_feature
from GenomeFileUtil.core.GenomeUtils import propagate_cds_props_to_gene, load_ontology_mappings
from GenomeFileUtil.core.GenomeInterface import GenomeInterface

# 3rd party imports
from Bio.Data import CodonTable
from Bio.Data.CodonTable import TranslationError
import Bio.SeqIO
from Bio.Seq import Seq

codon_table = CodonTable.ambiguous_generic_by_name["Standard"]
strand_table = str.maketrans("1?.", "+++")
MAX_MISC_FEATURE_SIZE = 10000


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
        self.version = re.search("module-version:\n\W+(.+)\n", yml_text).group(1)
        self.ont_mappings = load_ontology_mappings('/kb/module/data')
        self.code_table = 11
        self.skip_types = ('exon', 'five_prime_UTR', 'three_prime_UTR',
                           'start_codon', 'stop_codon', 'region', 'chromosome', 'scaffold')
        self.spoof_gene_count = 0
        self.is_phytozome = False
        self.strict = True
        self.generate_genes = False
        self.warnings = []
        self.feature_dict = collections.OrderedDict()
        self.cdss = set()
        self.ontologies_present = collections.defaultdict(dict)
        self.ontology_events = list()
        self.skiped_features = collections.Counter()
        self.feature_counts = collections.Counter()

    def warn(self, message):
        self.warnings.append(message)

    def generate_genome_json(self, params):
        # 1) validate parameters
        self._validate_import_file_params(params)
        self.code_table = params.get('genetic_code', 11)
        # 2) construct the input directory staging area
        input_directory = os.path.join(self.cfg.sharedFolder,
                                       'fast_gff_upload_' + str(uuid.uuid4()))
        os.makedirs(input_directory)
        file_paths = self._stage_input(params, input_directory)
        # 3) extract out the parameters
        params = self._set_parsed_params(params)
        if params.get('generate_missing_genes'):
            self.generate_genes = True

        # 4) do the upload
        genome = self._gen_genome_json(
            input_fasta_file=file_paths["fasta_file"],
            input_gff_file=file_paths["gff_file"],
            workspace_name=params['workspace_name'],
            core_genome_name=params['genome_name'],
            scientific_name=params['scientific_name'],
            source=params['source'],
            source_id=params['source_id'],
            release=params['release'],
        )
        if params.get('genetic_code'):
            genome["genetic_code"] = params['genetic_code']

        return genome, input_directory

    def import_file(self, params):

        genome, input_directory = self.generate_genome_json(params)

        json.dump(genome, open("{}/{}.json".format(self.cfg.sharedFolder,
                                                   genome['id']), 'w'),
                  indent=4)
        result = self.gi.save_one_genome({
            'workspace': params['workspace_name'],
            'name': params['genome_name'],
            'data': genome,
            "meta": params.get('metadata', {}),
        })
        report_string = 'A genome with {} contigs and the following feature ' \
                        'types was imported: {}'\
            .format(len(genome['contig_ids']), "\n".join(
                [k + ": " + str(v) for k, v in genome['feature_counts'].items()]))
        log(report_string)

        # 5) clear the temp directory
        shutil.rmtree(input_directory)

        # 6) return the result
        info = result['info']
        details = {
            'genome_ref': str(info[6]) + '/' + str(info[0]) + '/' + str(info[4]),
            'genome_info': info
        }

        return details

    def _gen_genome_json(self, input_gff_file=None, input_fasta_file=None,
                        workspace_name=None, core_genome_name=None,
                        scientific_name="unknown_taxon", source=None, source_id=None,
                        release=None):

        # reading in GFF file
        features_by_contig = self._retrieve_gff_file(input_gff_file)
        contig_ids = set()

        # parse feature information
        fasta_contigs = Bio.SeqIO.parse(input_fasta_file, "fasta")
        for contig in fasta_contigs:
            molecule_type = str(contig.seq.alphabet).replace(
                'IUPACAmbiguous', '').strip('()')
            contig_ids.add(contig.id)
            for feature in features_by_contig.get(contig.id, []):
                self._transform_feature(contig, feature)

        for cid in set(features_by_contig.keys()) - contig_ids:
            self.warn("Sequence name {} does not match a sequence id in the "
                      "FASTA file. {} features will not be imported."
                      .format(cid, len(features_by_contig[cid])))
            if self.strict:
                raise ValueError("Every feature sequence id must match a fasta sequence id")
        self._process_cdss()

        # save assembly file
        assembly_ref = self.au.save_assembly_from_fasta(
            {'file': {'path': input_fasta_file},
             'workspace_name': workspace_name,
             'assembly_name': core_genome_name + ".assembly"})
        assembly_data = self.dfu.get_objects(
            {'object_refs': [assembly_ref],
             'ignore_errors': 0})['data'][0]['data']

        # generate genome info
        genome = self._gen_genome_info(core_genome_name, scientific_name,
                                       assembly_ref, source, source_id, assembly_data,
                                       input_gff_file, molecule_type)
        genome['release'] = release
        if self.spoof_gene_count > 0:
            genome['warnings'] = genome.get('warnings', []) + \
                                    [warnings['spoofed_genome'].format(self.spoof_gene_count)]
            genome['suspect'] = 1

        return genome

    @staticmethod
    def _location(in_feature):
        in_feature['strand'] = in_feature['strand'].replace(
            "-1", "-").translate(strand_table)
        if in_feature['strand'] == '+':
            start = in_feature['start']
        elif in_feature['strand'] == '-':
            start = in_feature['end']
        else:
            raise ValueError('Invalid feature strand: {}'
                             .format(in_feature['strand']))
        return [
            in_feature['contig'],
            start,
            in_feature['strand'],
            in_feature['end'] - in_feature['start'] + 1
        ]

    @staticmethod
    def _validate_import_file_params(params):
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
                error_msg += str(list(file.keys()))
                raise ValueError(error_msg)
        if params.get('genetic_code'):
            if not (isinstance(params['genetic_code'], int) and 0 < params['genetic_code'] < 32):
                raise ValueError("Invalid genetic code specified: {}".format(params))

    def _set_parsed_params(self, params):
        log('Setting params')

        default_params = {
            'taxon_wsname': self.cfg.raw['taxon-workspace-name'],
            'scientific_name': 'unknown_taxon',
            'taxon_reference': None,
            'source': 'User',
            'release': None,
            'metadata': {},
            'source_id': 'unknown',
        }
        default_params.update(params)
        log(json.dumps(default_params, indent=1))
        return default_params

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
                log("staged input file =" + file_path)
                sys.stdout.flush()
                dfUtil_result = self.dfu.unpack_file({'file_path': file_path})
                file_paths[key] = dfUtil_result['file_path']
            else:
                raise ValueError('No valid files could be extracted based on the input')

        return file_paths

    def _retrieve_gff_file(self, input_gff_file):
        """
        _retrieve_gff_file: retrieve info from gff_file
    
        """
        log("Reading GFF file")
    
        feature_list = collections.defaultdict(list)
        is_patric = 0

        gff_file_handle = open(input_gff_file)
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
                if "phytozome" in source_id.lower():
                    self.is_phytozome = True

                #Checking to see if Phytozome
                if "PATRIC" in source_id:
                    is_patric = True

                #PATRIC prepends their contig ids with some gibberish
                if is_patric and "|" in contig_id:
                    contig_id = contig_id.split("|", 1)[1]

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
                        ftr['attributes'][key.lower()].append(parse.unquote(value.strip('"')))
                    elif(" " in attribute):
                        key, value = attribute.split(" ", 1)
                        ftr['attributes'][key.lower()].append(parse.unquote(value.strip('"')))
                    else:
                        pass
                        #log("Warning: attribute "+attribute+" cannot be separated into key,value pair")

                ftr['attributes']['raw'] = attributes
                if "id" in ftr['attributes']:
                    ftr['ID'] = ftr['attributes']['id'][0]
                if "parent" in ftr['attributes']:
                    ftr['Parent'] = ftr['attributes']['parent'][0]

                feature_list[contig_id].append(ftr)

            current_line = gff_file_handle.readline()

        gff_file_handle.close()

        #Some GFF/GTF files don't use "ID" so we go through the possibilities        
        feature_list = self._add_missing_identifiers(feature_list)

        #Most bacterial files have only CDSs
        #In order to work with prokaryotic and eukaryotic gene structure synonymously
        #Here we add feature dictionaries representing the parent gene and mRNAs
        #feature_list = self._add_missing_parents(feature_list)

        #Phytozome has the annoying habit of editing their identifiers so we fix them
        if self.is_phytozome:
            self._update_phytozome_features(feature_list)

        #All identifiers need to be checked so that they follow the same general rules
        #Rules are listed within the function itself
        feature_list = self._update_identifiers(feature_list)

        return feature_list

    def _add_missing_identifiers(self, feature_list):
        log("Adding missing identifiers")
        #General rule is to iterate through a range of possibilities if "ID" is missing
        for contig in feature_list:
            for i, feat in enumerate(feature_list[contig]):
                if "ID" not in feature_list[contig][i]:
                    for key in ("transcriptid", "proteinid", "pacid",
                                "parent", "name", 'transcript_id'):
                        if key in feature_list[contig][i]['attributes']:
                            feature_list[contig][i]['ID'] = feature_list[
                                contig][i]['attributes'][key][0]
                            break
                    if feat['type'] not in self.skip_types:
                        self.feature_counts[feat['type']] += 1

                    #If the process fails, throw an error
                    if "ID" not in feature_list[contig][i]:
                        feat['ID'] = "{}_{}".format(feat['type'],
                                                    self.feature_counts[feat['type']])
                        #log("Warning: Could find unique ID to utilize in GFF attributes: {}. "
                        #    "ID '{}' has been assigned".format(feat['attributes'], feat['ID']))
        return feature_list

    def _add_missing_parents(self, feature_list):

        #General rules is if CDS or RNA missing parent, add them
        for contig in feature_list:
            ftrs = feature_list[contig]
            new_ftrs = []
            for i in range(len(ftrs)):
                if ftrs[i]["type"] in self.skip_types:
                    continue
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
        for contig in feature_list:
            feature_position_dict = {}
            for i in range(len(feature_list[contig])):

                #Maintain old_id for reference
                #Sometimes ID isn't available, so use PACid
                old_id = None
                for key in ("id", "pacid"):
                    if(key in feature_list[contig][i]['attributes']):
                        old_id = feature_list[contig][i]['attributes'][key][0]
                        break
                if(old_id is None):
                    #This should be an error
                    #log("Cannot find unique ID, PACid, or pacid in GFF "
                    #    "attributes: " + feature_list[contig][i][contig])
                    continue

                #Retain old_id
                feature_position_dict[old_id]=i

                # Clip off the increment on CDS IDs so fragments of the same
                # CDS share the same ID
                if "CDS" in feature_list[contig][i]["ID"]:
                    feature_list[contig][i]["ID"] = feature_list[contig][i]["ID"].rsplit('.', 1)[0]

                #In Phytozome, gene and mRNA have "Name" field, CDS do not
                if("name" in feature_list[contig][i]['attributes']):
                    feature_list[contig][i]["ID"] = feature_list[contig][i]['attributes']['name'][0]

                if("Parent" in feature_list[contig][i]):
                    #Update Parent to match new ID of parent ftr
                    feature_list[contig][i]["Parent"] = feature_list[contig][feature_position_dict[feature_list[contig][i]["Parent"]]]["ID"]

        return feature_list

    def _update_identifiers(self, feature_list):

        #General rules:
        #1) Genes keep identifier
        #2) RNAs keep identifier only if its different from gene, otherwise append ".mRNA"
        #3) CDS always uses RNA identifier with ".CDS" appended

        mRNA_parent_dict = dict()

        for contig in feature_list:
            for ftr in feature_list[contig]:
                if ftr["type"] in self.skip_types:
                    continue
                if("Parent" in ftr):
                    #Retain old_id of parents
                    old_id = ftr["ID"]

                    if(ftr["ID"] == ftr["Parent"] or "CDS" in ftr["type"]):
                        ftr["ID"] = ftr["Parent"]+"."+ftr["type"]

                    #link old to new ids for mRNA to use with CDS
                    if("RNA" in ftr["type"]):
                        mRNA_parent_dict[old_id]=ftr["ID"]

        return feature_list

    def _check_location_order(self,locations):
        """If order looks good return None.  
           If out of order return warning
           If on multiple strands return warning"""
        strand = None
        last_start = 0
        for location in locations:
            if strand == None:
                strand = location[2]
            elif strand != location[2]:
                return warnings["both_strand_coordinates"]
        if strand == "-":
            locations = reversed(locations)
        for location in locations:
            if last_start > location[1]:
                return warnings["out_of_order"]
            else:
                last_start = location[1]
        return None        

    def _create_ontology_event(self, ontology_type):
        """Creates the ontology_event if necessary
        Returns the index of the ontology event back."""
        if ontology_type not in self.ont_mappings:
            raise ValueError("{} is not a supported ontology".format(ontology_type))

        if "event_index" not in self.ont_mappings[ontology_type]:
            self.ont_mappings[ontology_type]['event_index'] = len(self.ontology_events)
            if ontology_type == "GO":
                ontology_ref = "KBaseOntology/gene_ontology"
            elif ontology_type == "PO":
                ontology_ref = "KBaseOntology/plant_ontology"
            else:
                ontology_ref = f"KBaseOntology/{ontology_type.lower()}_ontology"
            self.ontology_events.append({
                "method": "GenomeFileUtils Genbank uploader from annotations",
                "method_version": self.version,
                "timestamp": self.time_string,
                "id": ontology_type,
                "ontology_ref": ontology_ref
            })

        return self.ont_mappings[ontology_type]['event_index']

    def _get_ontology_db_xrefs(self, feature):
        """Splits the ontology info from the other db_xrefs"""
        ontology = collections.defaultdict(dict)
        db_xrefs = []
        # these are keys are formatted strangely and require special parsing
        for key in ("go_process", "go_function", "go_component"):
            ontology_event_index = self._create_ontology_event("GO")
            for term in feature.get(key, []):
                sp = term.split(" - ")
                ontology['GO'][sp[0]] = [ontology_event_index]
                self.ontologies_present['GO'][sp[0]] = self.ont_mappings['GO'].get(sp[0], '')

        # CATH terms are not distinct from EC numbers so myst be found by key
        for term in feature.get('cath_funfam', []) + feature.get('cath', []):
            for ref in term.split(','):
                ontology['CATH'][ref] = [self._create_ontology_event("CATH")]
                self.ontologies_present['CATH'][ref] = self.ont_mappings['CATH'].get(ref, '')

        search_keys = ['ontology_term', 'db_xref', 'dbxref', 'product_source', 'tigrfam', 'pfam',
                       'cog', 'go', 'po', 'ko']
        ont_terms = []
        # flatten out into list of values
        for key in search_keys:
            if key in feature:
                ont_terms += [x for y in feature[key] for x in y.split(',')]

        for ref in ont_terms:
            if ref.startswith('GO:'):
                ontology['GO'][ref] = [self._create_ontology_event("GO")]
                self.ontologies_present['GO'][ref] = self.ont_mappings['GO'].get(ref, '')
            elif ref.startswith('PO:'):
                ontology['PO'][ref] = [self._create_ontology_event("PO")]
                self.ontologies_present['PO'][ref] = self.ont_mappings['PO'].get(ref, '')
            elif ref.startswith('KO:'):
                ontology['KO'][ref] = [self._create_ontology_event("KO")]
                self.ontologies_present['KO'][ref] = self.ont_mappings['KO'].get(ref, '')
            elif ref.startswith('COG'):
                ontology['COG'][ref] = [self._create_ontology_event("COG")]
                self.ontologies_present['COG'][ref] = self.ont_mappings['COG'].get(ref, '')
            elif ref.startswith('PF'):
                ontology['PFAM'][ref] = [self._create_ontology_event("PFAM")]
                self.ontologies_present['PFAM'][ref] = self.ont_mappings['PFAM'].get(ref, '')
            elif ref.startswith('TIGR'):
                ontology['TIGRFAM'][ref] = [self._create_ontology_event("TIGRFAM")]
                self.ontologies_present['TIGRFAM'][ref] = self.ont_mappings['TIGRFAM'].get(ref, '')
            else:
                db_xrefs.append(tuple(ref.split(":", 1)))
        return dict(ontology), db_xrefs

    def _transform_feature(self, contig, in_feature):
        """Converts a feature from the gff ftr format into the appropriate
        format for a genome object """
        def _aliases(feat):
            keys = ('locus_tag', 'old_locus_tag', 'protein_id',
                    'transcript_id', 'gene', 'ec_number', 'gene_synonym')
            alias_list = []
            for key in keys:
                if key in feat['attributes']:
                    alias_list.extend([(key, val) for val in feat['attributes'][key]])
            return alias_list

        if in_feature['start'] < 1 or in_feature['end'] > len(contig):
            self.warn("Feature with invalid location for specified "
                      "contig: " + str(in_feature))
            if self.strict:
                raise ValueError("Features must be completely contained within the Contig in the "
                                 "Fasta file. Feature: " + str(in_feature))
            return

        feat_seq = contig.seq[in_feature['start']-1:in_feature['end']].upper()
        if in_feature['strand'] in {'-', '-1'}:
            feat_seq = feat_seq.reverse_complement()

        # if the feature ID is duplicated (CDS or transpliced gene) we only
        # need to update the location and dna_sequence
        if in_feature.get('ID') in self.feature_dict:
            existing = self.feature_dict[in_feature['ID']]
            existing['location'].append(self._location(in_feature))
            existing['dna_sequence'] = existing.get('dna_sequence', '') + str(feat_seq)
            existing['dna_sequence_length'] = len(existing['dna_sequence'])
            return

        # The following is common to all the feature types
        out_feat = {
            "id": in_feature.get('ID'),
            "type": in_feature['type'],
            "location": [self._location(in_feature)],
            "dna_sequence": str(feat_seq),
            "dna_sequence_length": len(feat_seq),
            "md5": hashlib.md5(str(feat_seq).encode('utf8')).hexdigest(),
        }

        # add optional fields
        if 'note' in in_feature['attributes']:
            out_feat['note'] = in_feature['attributes']["note"][0]
        ont, db_xrefs = self._get_ontology_db_xrefs(in_feature['attributes'])
        if ont:
            out_feat['ontology_terms'] = ont
        aliases = _aliases(in_feature)
        if aliases:
            out_feat['aliases'] = aliases
        if db_xrefs:
            out_feat['db_xrefs'] = db_xrefs
        if 'product' in in_feature['attributes']:
            out_feat['functions'] = in_feature['attributes']["product"]
        if 'product_name' in in_feature['attributes']:
            if "functions" in out_feat:
                out_feat['functions'].extend(in_feature['attributes']["product_name"])
            else:
                out_feat['functions'] = in_feature['attributes']["product_name"]
        if 'function' in in_feature['attributes']:
            out_feat['functional_descriptions'] = in_feature['attributes']["function"]
        if 'inference' in in_feature['attributes']:
            GenomeUtils.parse_inferences(in_feature['attributes']['inference'])
        if 'trans-splicing' in in_feature['attributes'].get('exception', []):
            out_feat['flags'] = out_feat.get('flags',[]) + ['trans_splicing']
        if 'pseudo' in in_feature['attributes'].get('exception', []):
            out_feat['flags'] = out_feat.get('flags',[]) + ['pseudo']
        if 'ribosomal-slippage' in in_feature['attributes'].get('exception', []):
            out_feat['flags'] = out_feat.get('flags',[]) + ['ribosomal_slippage']
        parent_id = in_feature.get('Parent', '')
        if parent_id and parent_id not in self.feature_dict:
            raise ValueError("Parent ID: {} was not found in feature ID list.".format(parent_id))

        # if the feature is a exon or UTR, it will only be used to update the
        # location and sequence of it's parent, we add the info to it parent
        # feature but not the feature dict
        if in_feature['type'] in self.skip_types:
            if parent_id and in_feature['type'] in {'exon', 'five_prime_UTR', 'three_prime_UTR'}:
                parent = self.feature_dict[parent_id]
                if in_feature['type'] not in parent:
                    parent[in_feature['type']] = []
                parent[in_feature['type']].append(out_feat)
            return

        # add type specific features
        elif 'gene' in in_feature['type']:
            out_feat['protein_translation_length'] = 0
            out_feat['cdss'] = []

        elif in_feature['type'] == 'CDS':
            if parent_id:
                parent = self.feature_dict[parent_id]
                if 'cdss' in parent:  # parent must be a gene
                    if not is_parent(parent,out_feat):
                        parent["warnings"] = parent.get('warnings', [] ) + [
                            warnings["genes_CDS_child_fails_location_validation"].format(out_feat["id"])]
                        out_feat["warnings"] = out_feat.get('warnings', []) + [
                            warnings["CDS_fail_child_of_gene_coordinate_validation"].format(parent_id)]      
                    parent['cdss'].append(in_feature['ID'])
                    out_feat['parent_gene'] = parent_id             
                else:  # parent must be mRNA
                    if not is_parent(parent,out_feat):
                        parent["warnings"] = parent.get('warnings', [] ) + [
                            warnings["mRNA_fail_parent_coordinate_validation"].format(out_feat["id"])]
                        out_feat["warnings"] = out_feat.get('warnings', []) + [
                            warnings["CDS_fail_child_of_mRNA_coordinate_validation"].format(parent_id)]      
                    parent['cds'] = in_feature['ID']
                    out_feat['parent_mrna'] = parent_id
                    parent_gene = self.feature_dict[parent['parent_gene']]
                    parent_gene['cdss'].append(in_feature['ID'])
                    out_feat['parent_gene'] = parent['parent_gene']
            # keep track of CDSs for post processing
            self.cdss.add(out_feat['id'])

        elif in_feature['type'] == 'mRNA':
            if parent_id:
                parent = self.feature_dict[parent_id]
                if 'mrnas' not in parent:                         
                    parent['mrnas'] = []
                if 'cdss' in parent:  # parent must be a gene
                    parent['mrnas'].append(in_feature['ID'])
                    out_feat['parent_gene'] = parent_id
                if not is_parent(parent,out_feat):
                    parent["warnings"] = parent.get('warnings', [] ) + [
                        warnings["genes_mRNA_child_fails_location_validation"].format(out_feat["id"])]
                    out_feat["warnings"] = out_feat.get('warnings', []) + [
                        warnings["mRNAs_parent_gene_fails_location_validation"].format(parent_id)]    

        else:
            out_feat["type"] = in_feature['type']
            # this prevents big misc_features from blowing up the genome size
            if out_feat['dna_sequence_length'] > MAX_MISC_FEATURE_SIZE:
                del out_feat['dna_sequence']
            if parent_id:
                parent = self.feature_dict[parent_id]
                if 'children' not in parent:
                    parent['children'] = []
                parent['children'].append(out_feat['id'])
                out_feat['parent_gene'] = parent_id
                if not is_parent(parent,out_feat):
                    parent["warnings"] = parent.get('warnings', [] ) + [
                        warnings["generic_parents_child_fails_location_validation"].format(out_feat["id"])]
                    out_feat["warnings"] = out_feat.get('warnings', []) + [
                        warnings["generic_childs_parent_fails_location_validation"].format(parent_id)]    

        self.feature_dict[out_feat['id']] = out_feat

    def _process_cdss(self):
        """Because CDSs can have multiple fragments, it's necessary to go
        back over them to calculate a final protein sequence"""
        for cds_id in self.cdss:
            cds = self.feature_dict[cds_id]
            try:
                prot_seq = str(Seq(cds['dna_sequence']).translate(
                            self.code_table, cds=True).strip("*"))
            except TranslationError as e:
                cds['warnings'] = cds.get('warnings', []) + [str(e)]
                prot_seq = ""

            cds.update({
                "protein_translation": prot_seq,
                "protein_md5": hashlib.md5(prot_seq.encode('utf8')).hexdigest(),
                "protein_translation_length": len(prot_seq),
            })
            if 'parent_gene' in cds:
                parent_gene = self.feature_dict[cds['parent_gene']]
                # no propigation for now
                propagate_cds_props_to_gene(cds, parent_gene)
            elif self.generate_genes:
                spoof = copy.copy(cds)
                spoof['type'] = 'gene'
                spoof['id'] = cds['id']+"_gene"
                spoof['cdss'] = [cds['id']]
                spoof['warnings'] = [warnings['spoofed_gene'].format(cds['id'])]
                self.feature_dict[spoof['id']] = spoof
                cds['parent_gene'] = spoof['id']
                self.spoof_gene_count += 1
            else:
                raise ValueError(warnings['no_spoof'])

            self.feature_dict[cds['id']] = cds

    def _update_from_exons(self, feature):
        """This function updates the sequence and location of a feature based
            on it's UTRs, CDSs and exon information"""
        # note that start and end here are in direction of translation
        def start(loc):
            return loc[0][1]

        def end(loc):
            if loc[-1][2] == "+":
                return loc[-1][1] + loc[-1][3] + 1
            else:
                return loc[-1][1] - loc[-1][3] - 1

        if 'exon' in feature:
            # update the feature with the exon locations and sequences
            feature['location'] = [x['location'][0] for x in feature['exon']]
            feature['dna_sequence'] = "".join(
                x['dna_sequence'] for x in feature['exon'])
            feature['dna_sequence_length'] = len(feature['dna_sequence'])

        # construct feature location from utrs and cdss if present
        elif 'cds' in feature:
            cds = [copy.deepcopy(self.feature_dict[feature['cds']])]
            locs = []
            seq = ""
            for frag in feature.get('five_prime_UTR', []) + cds + \
                    feature.get('three_prime_UTR', []):

                # merge into last location if adjacent
                if locs and abs(end(locs) - start(frag['location'])) == 1:
                    # extend the location length by the length of the first
                    # location in the fragment
                    first = frag['location'].pop(0)
                    locs[-1][3] += first[3]

                locs.extend(frag['location'])
                seq += frag['dna_sequence']

            feature['location'] = locs
            feature['dna_sequence'] = seq
            feature['dna_sequence_length'] = len(seq)

        # remove these properties as they are no longer needed
        for x in ['five_prime_UTR', 'three_prime_UTR', 'exon']:
            feature.pop(x, None)

        else:
            ValueError('Feature {} must contain either exon or cds data to '
                       'construct an accurate location and sequence'.format(
                        feature['id']))

    def _gen_genome_info(self, core_genome_name, scientific_name, assembly_ref,
                         source, source_id, assembly, input_gff_file, molecule_type):
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
        genome['ontologies_present'] = dict(self.ontologies_present)
        genome['ontology_events'] = self.ontology_events
        genome['taxonomy'], genome['taxon_ref'], genome['domain'], \
            genome["genetic_code"] = self.gi.retrieve_taxon(self.taxon_wsname,
                                                            genome['scientific_name'])
        genome['source'], genome['genome_tiers'] = self.gi.determine_tier(
            source)
        genome['source_id'] = source_id

        # Phytozome gff files are not compatible with the RNASeq Pipeline
        # so it's better to build from the object than cache the file
        if self.is_phytozome:
            gff_file_to_shock = self.dfu.file_to_shock(
                {'file_path': input_gff_file, 'make_handle': 1, 'pack': "gzip"})
            genome['gff_handle_ref'] = gff_file_to_shock['handle']['hid']

        for feature in self.feature_dict.values():
            self.feature_counts[feature['type']] += 1
            if 'exon' in feature or feature['type'] == 'mRNA':
                self._update_from_exons(feature)

            # Test if location order is in order.
            is_transpliced = "flags" in feature and "trans_splicing" in feature["flags"]
            if not is_transpliced and len(feature["location"]) > 1:
                # Check the order only if not trans_spliced and has more than 1 location.
                location_warning = self._check_location_order(feature["location"])
                if location_warning is not None:
                    feature["warnings"] = feature.get('warnings', []) + [location_warning]

            contig_len = genome["contig_lengths"][genome["contig_ids"].index(feature["location"][0][0])]
            feature = check_full_contig_length_or_multi_strand_feature(
                feature, is_transpliced,contig_len, self.skip_types)

            # sort features into their respective arrays
            if feature['type'] == 'CDS':
                del feature['type']
                genome['cdss'].append(feature)
            elif feature['type'] == 'mRNA':
                del feature['type']
                genome['mrnas'].append(feature)
            elif feature['type'] == 'gene':
                # remove duplicates that may arise from CDS info propagation
                for key in ('functions', 'aliases', 'db_xrefs'):
                    if key in feature:
                        feature[key] = list(set(feature[key]))
                if feature['cdss']:
                    del feature['type']
                    self.feature_counts["protein_encoding_gene"] += 1
                    genome['features'].append(feature)
                else:
                    feature.pop('mrnas', None)
                    feature.pop('cdss', None)
                    feature.pop('protein_translation_length', None)
                    self.feature_counts["non_coding_features"] += 1
                    genome['non_coding_features'].append(feature)
            else:
                genome['non_coding_features'].append(feature)

        if self.warnings:
            genome['warnings'] = self.warnings
        genome['feature_counts'] = dict(self.feature_counts)
        return genome
