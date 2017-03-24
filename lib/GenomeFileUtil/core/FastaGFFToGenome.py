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
            # provenance=ctx['provenance'],
            # usermeta=parsed_params['metadata']
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

        # retrieve feature identifiers
        (features_identifiers_dict,
         features_identifiers_list,
         features_identifiers_count) = self._retrieve_feature_identifiers(feature_list)

        (updated_features_identifiers_dict,
         updated_features_list,
         updated_features_identifiers_count) = self._update_feature_identifiers(
                                                                features_identifiers_dict,
                                                                features_identifiers_list,
                                                                features_identifiers_count)

        # retrieve genome feature list
        (genome_features_list,
         genome_cdss_list,
         genome_mrnas_list) = self._retrieve_genome_feature_list(
                                                                updated_features_identifiers_dict,
                                                                updated_features_list,
                                                                updated_features_identifiers_count,
                                                                assembly)

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
        if 'type' in params and params['type'] not in valid_types:
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

        header = list()
        feature_list = dict()
        original_CDS_count = dict()
        original_feature_ids = dict()

        gff_file_handle = open(input_gff_file, 'rb')
        current_line = gff_file_handle.readline()
        gff_object = dict()
        while (current_line != ''):
            current_line = current_line.strip()

            if(current_line.startswith("##") or current_line.startswith("#!")):
                header.append(current_line)
                if('headers' not in gff_object):
                    gff_object['headers'] = list()
                gff_object['headers'].append(current_line)
            else:
                if('features' not in gff_object):
                    gff_object['features'] = list()

                (contig_id, source_id, feature_type, start, end,
                 score, strand, phase, attributes) = current_line.split('\t')
                attributes_dict = dict()
                for attribute in attributes.split(";"):
                    if(attribute == "" or "=" not in attribute):
                        continue
                    key, value = attribute.split("=", 1)
                    attributes_dict[key] = value

                # ID should be transferred from Name or Parent
                old_id = None
                for key in ("ID", "PACid", "pacid"):
                    if(key in attributes_dict):
                        old_id = attributes_dict[key]
                        break
                if(old_id is None):
                    print ("Cannot find unique ID, PACid, or pacid in GFF attributes: "+attributes)
                    continue

                if("Name" in attributes_dict):
                    attributes_dict["ID"] = attributes_dict["Name"]
                else:
                    attributes_dict["ID"] = original_feature_ids[
                                                    attributes_dict["Parent"]]+"."+feature_type

                    # if CDS have to increment
                    if(feature_type == "CDS"):
                        if(attributes_dict["ID"] not in original_CDS_count):
                            original_CDS_count[attributes_dict["ID"]] = 1
                        else:
                            original_CDS_count[attributes_dict["ID"]] += 1

                        attributes_dict["ID"] += "."+str(original_CDS_count[attributes_dict["ID"]])

                # Update parent
                if("Parent" in attributes_dict):
                    attributes_dict["Parent"] = original_feature_ids[attributes_dict["Parent"]]

                original_feature_ids[old_id] = attributes_dict["ID"]

                # recreate line for GFF
                partial_line, attributes = current_line.rsplit('\t', 1)
                new_line = partial_line + "\t" + ";".join(key+"="+attributes_dict[key] for key in
                                                          attributes_dict.keys())
                gff_object['features'].append(new_line)

                # if(contig_id not in assembly["contigs"]):
                #     log("Missing contig: "+contig_id)

                if(contig_id not in feature_list):
                    feature_list[contig_id] = list()

                feature = {'type': feature_type, 'start': int(start), 'end': int(end),
                           'score': score, 'strand': strand, 'phase': phase}
                for attribute in attributes.split(";"):
                    if(attribute == "" or "=" not in attribute):
                        continue
                    key, value = attribute.split("=", 1)
                    feature[key] = value

                # Append contig identifier
                feature["contig"] = contig_id
                feature_list[contig_id].append(feature)

            current_line = gff_file_handle.readline()
        gff_file_handle.close()

        # Writing updated lines to gff_file_handle
        input_gff_file = input_gff_file.replace("gene", "edited_gene")
        gff_file_handle = gzip.open(input_gff_file, 'wb')
        if('headers' in gff_object):
            gff_file_handle.write("\n".join(gff_object["headers"]))
        gff_file_handle.write("\n".join(gff_object["features"]))
        gff_file_handle.close()

        return feature_list

    def _retrieve_feature_identifiers(self, feature_list):
        features_identifiers_dict = dict()
        features_identifiers_list = list()
        features_identifiers_count = dict()
        features_parents_dict = dict()
        features_name_id_dict = dict()
        CDS_count = dict()
        for contig in sorted(feature_list):
            for feature in feature_list[contig]:
                # We're only considering gene, mRNA, and CDS for brevity's sake
                if(feature["type"] not in ("gene", "mRNA", "CDS")):
                    continue

                # gene and mRNA always have name, CDS do not
                if("Name" not in feature):
                    feature["Name"] = None

                # Update parent following name/id switch
                if("Parent" in feature and feature["Parent"] in features_name_id_dict):
                    feature["Parent"] = features_name_id_dict[feature["Parent"]]

                # ID should be transferred to Name, but need to maintain parent
                if(feature["Name"] is not None):
                    features_name_id_dict[feature["ID"]] = feature["Name"]
                    feature["ID"] = feature["Name"]
                else:
                    feature["ID"] = feature["Parent"]+"."+feature["type"]
                    # if CDS have to increment
                    if(feature["type"] == "CDS"):
                        if(feature["ID"] not in CDS_count):
                            CDS_count[feature["ID"]] = 1
                        else:
                            CDS_count[feature["ID"]] += 1

                        feature["ID"] += "."+str(CDS_count[feature["ID"]])

                # Collect
                if(feature["type"] == "gene"):
                    features_identifiers_dict[feature["ID"]] = dict()
                if(feature["type"] == "mRNA"):
                    features_identifiers_dict[feature["Parent"]][feature["ID"]] = dict()
                    features_parents_dict[feature["ID"]] = feature["Parent"]
                if(feature["type"] == "CDS"):
                    features_identifiers_dict[features_parents_dict[feature[
                                                "Parent"]]][feature["Parent"]][feature["ID"]] = 1

                features_identifiers_list.append(feature)
                features_identifiers_count[feature["ID"]] = len(features_identifiers_list)-1

        return features_identifiers_dict, features_identifiers_list, features_identifiers_count

    def _update_feature_identifiers(self, features_identifiers_dict,
                                    features_identifiers_list, features_identifiers_count):

        updated_features_identifiers_dict = dict()
        updated_features_list = list()
        updated_features_identifiers_count = dict()
        updated_features_parents_dict = dict()
        updated_CDS_count = dict()
        for gene in sorted(features_identifiers_dict):

            # retrieve original object
            gene_ftr = features_identifiers_list[features_identifiers_count[gene]]

            # store gene
            updated_features_identifiers_dict[gene_ftr["ID"]] = dict()
            updated_features_list.append(gene_ftr)
            updated_features_identifiers_count[gene_ftr["ID"]] = len(updated_features_list)-1

            for mRNA in sorted(features_identifiers_dict[gene], key=lambda x:
                               features_identifiers_count[x]):
                # retrieve feature
                mRNA_ftr = features_identifiers_list[features_identifiers_count[mRNA]]

                if("PAC" in mRNA[0:3]):
                    if("Name" in mRNA_ftr):
                        mRNA_ftr["ID"] = mRNA_ftr["Name"]

                updated_features_identifiers_dict[gene_ftr["ID"]][mRNA_ftr["ID"]] = dict()
                updated_features_parents_dict[mRNA_ftr["ID"]] = mRNA_ftr["Parent"]

                updated_features_list.append(mRNA_ftr)
                updated_features_identifiers_count[mRNA_ftr["ID"]] = len(updated_features_list)-1

                for CDS in sorted(features_identifiers_dict[gene][mRNA], key=lambda x:
                                  features_identifiers_count[x]):
                    # retrieve feature
                    CDS_ftr = features_identifiers_list[features_identifiers_count[CDS]]

                    if("PAC" in CDS[0:3]):
                        CDS_ftr["ID"] = mRNA_ftr["ID"]+".CDS"

                        if(CDS_ftr["ID"] not in updated_CDS_count):
                            updated_CDS_count[CDS_ftr["ID"]] = 1
                        else:
                            updated_CDS_count[CDS_ftr["ID"]] += 1

                        CDS_ftr["ID"] += "."+str(updated_CDS_count[CDS_ftr["ID"]])
                        CDS_ftr["Parent"] = mRNA_ftr["ID"]

                    updated_features_identifiers_dict[gene_ftr["ID"]][
                                                            mRNA_ftr["ID"]][CDS_ftr["ID"]] = 1
                    updated_features_parents_dict[CDS_ftr["ID"]] = CDS_ftr["Parent"]

                    updated_features_list.append(CDS_ftr)
                    updated_features_identifiers_count[CDS_ftr["ID"]] = len(updated_features_list)-1

        return (updated_features_identifiers_dict,
                updated_features_list,
                updated_features_identifiers_count)

    def _retrieve_genome_feature_list(self, updated_features_identifiers_dict,
                                      updated_features_list, updated_features_identifiers_count,
                                      assembly):

        genome_features_list = list()
        genome_mrnas_list = list()
        genome_cdss_list = list()
        for gene in sorted(updated_features_identifiers_dict):
            # retrieve updated object
            gene_ftr = updated_features_list[updated_features_identifiers_count[gene]]

            gene_object = self._convert_ftr_object(gene_ftr,
                                                   assembly["contigs"][gene_ftr[
                                                                    "contig"]]["sequence"])
            gene_object["type"] = "gene"

            # New terms, TODO, move to end of gene loop
            gene_object["cdss"] = list()
            gene_object["mrnas"] = list()

            # use function of longest CDS for gene
            longest_protein_length = 0
            longest_protein_sequence = ""
            for mRNA in sorted(updated_features_identifiers_dict[gene], key=lambda x:
                               updated_features_identifiers_count[x]):
                # retrieve updated object
                mRNA_ftr = updated_features_list[updated_features_identifiers_count[mRNA]]

                feature_object = self._convert_ftr_object(mRNA_ftr,
                                                          assembly["contigs"][mRNA_ftr[
                                                                        "contig"]]["sequence"])
                feature_object['parent_gene'] = gene_object['id']

                mrna_object = copy.deepcopy(feature_object)
                cds_object = copy.deepcopy(feature_object)

                cds_object['id'] = mrna_object['id']+".CDS"
                mrna_object['cds'] = cds_object['id']

                cds_object['parent_mrna'] = mrna_object['id']

                del mrna_object["dna_sequence"]
                del mrna_object["dna_sequence_length"]

                cds_object["ontology_terms"] = dict()

                gene_object["mrnas"].append(mrna_object["id"])
                gene_object["cdss"].append(cds_object["id"])

                # CDS aggregation needs to be done to build protein sequence and list of locations
                CDS_list = sorted(updated_features_identifiers_dict[gene][mRNA], key=lambda x:
                                  updated_features_identifiers_count[x])

                dna_sequence = ""
                locations = list()

                # collect phases, and lengths of exons
                # right now, this is only for the purpose of error reporting
                phases = list()
                exons = list()

                for CDS in (CDS_list):
                    # retrieve updated partial CDS
                    add_ftr = updated_features_list[updated_features_identifiers_count[CDS]]
                    phases.append(add_ftr["phase"])

                    add_ftr_obj = self._convert_ftr_object(add_ftr,
                                                           assembly["contigs"][
                                                                    add_ftr["contig"]]["sequence"])
                    exons.append(len(add_ftr_obj["dna_sequence"]))

                    # Remove base(s) according to phase, but only for first CDS
                    if(CDS == CDS_list[0] and int(add_ftr["phase"]) != 0):
                        log("Adjusting phase for first CDS: "+CDS)
                        add_ftr_obj["dna_sequence"] = add_ftr_obj["dna_sequence"][int(
                                                                                add_ftr["phase"]):]

                    dna_sequence += add_ftr_obj["dna_sequence"]
                    locations.append(add_ftr_obj["location"][0])

                # translate sequence
                dna_sequence_obj = Seq(dna_sequence, IUPAC.ambiguous_dna)
                rna_sequence = dna_sequence_obj.transcribe()

                # incomplete gene model with no start codon
                if str(rna_sequence.upper())[:3] not in codon_table.start_codons:
                    msg = "Missing start codon for {} Assuming incomplete gene model.".format(
                                                                            feature_object["id"])
                    log(msg)

                # You should never have this problem, needs to be reported rather than "fixed"
                codon_count = len(str(rna_sequence)) % 3
                if codon_count != 0:
                    msg = "Number of bases for RNA sequence for {} ".format(feature_object["id"])
                    msg += "is not divisible by 3. "
                    msg += "The resulting protein may well be mis-translated."
                    log(msg)

                protein_sequence = Seq("")
                try:
                    protein_sequence = rna_sequence.translate()
                except CodonTable.TranslationError as te:
                    log("TranslationError for: "+feature_object["id"], phases, exons, " : "+str(te))

                cds_object["protein_translation"] = str(protein_sequence).upper()
                cds_object["protein_translation_length"] = len(cds_object["protein_translation"])
                cds_object["md5"] = hashlib.md5(cds_object["protein_translation"]).hexdigest()

                if(cds_object["protein_translation_length"] > longest_protein_length):
                    longest_protein_length = cds_object["protein_translation_length"]
                    longest_protein_sequence = cds_object["protein_translation"]

                del cds_object["dna_sequence"]
                del cds_object["dna_sequence_length"]
                if("aliases" not in cds_object):
                    cds_object["aliases"] = list()
                if("function" not in cds_object):
                    cds_object["function"] = ""

                # End of mRNA loop
                genome_mrnas_list.append(mrna_object)
                genome_cdss_list.append(cds_object)

            # End of gene loop
            gene_object["ontology_terms"] = dict()
            gene_object["protein_translation"] = longest_protein_sequence
            gene_object["protein_translation_length"] = longest_protein_length
            genome_features_list.append(gene_object)

        return genome_features_list, genome_cdss_list, genome_mrnas_list

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
        new_ftr["location"] = [[old_ftr["contig"], old_ftr["start"], old_ftr["strand"],
                                len(dna_sequence)]]
        return new_ftr
