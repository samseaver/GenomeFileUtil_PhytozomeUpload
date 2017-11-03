
import datetime
import time
import os
import re
import sys
import shutil
import uuid
import hashlib
import itertools
from collections import Counter, defaultdict

import Bio.SeqIO
import Bio.SeqUtils

from AssemblyUtil.AssemblyUtilClient import AssemblyUtil
from DataFileUtil.DataFileUtilClient import DataFileUtil
from GenomeFileUtil.core.GenbankUploaderScript import upload_genome
from GenomeInterface import GenomeInterface
from KBaseReport.KBaseReportClient import KBaseReport


class GenbankToGenome:
    def __init__(self, config):
        self.cfg = config
        self.gi = GenomeInterface(config)
        self.dfu = DataFileUtil(config.callbackURL)
        self.aUtil = AssemblyUtil(config.callbackURL)
        self.report_client = KBaseReport(config.callbackURL)
        self._messages = []
        self.time_string = str(datetime.datetime.fromtimestamp(
            time.time()).strftime('%Y_%m_%d_%H_%M_%S'))
        yml_text = open('/kb/module/kbase.yml').read()
        self.version = re.search("module-version:\n\W+(.+)\n", yml_text).group(1)
        self.ontologies_present = defaultdict(dict)
        self.default_params = {
            'source': 'Genbank',
            'taxon_wsname': self.cfg.raw['taxon-workspace-name'],
            'taxon_lookup_obj_name': self.cfg.raw['taxon-lookup-object-name'],
            'taxon_reference': None,

            'ontology_wsname': self.cfg.raw['ontology-workspace-name'],
            'ontology_GO_obj_name': self.cfg.raw['ontology-gene-ontology-obj-name'],
            'ontology_PO_obj_name': self.cfg.raw['ontology-plant-ontology-obj-name'],

            'release': None,
            'genetic_code': 11,
            'generate_ids_if_needed': 0,
            'exclude_ontologies': 0,
            'type': 'User upload',
            'metadata': {}
        }

    def log(self, message):
        self._messages.append(message)
        print message

    @property
    def messages(self):
        return "\n".join(self._messages)

    def import_file(self, ctx, params):

        # 1) validate parameters and extract defaults
        self.validate_params(params)

        # 2) construct the input directory staging area
        input_directory = self.stage_input(params)

        # 3) extract out the parameters
        parsed_params = self.default_params
        parsed_params.update(params)

        # add any optional parameters
        optional_param_fields_to_check = [
                'source',
                'taxon_wsname',
                'taxon_reference',
                'release',
                'genetic_code',
                'generate_ids_if_needed',
                'exclude_ontologies',
                'type',
                'metadata'
            ]

        for field in optional_param_fields_to_check:
            if field in params:
                parsed_params[field] = params[field]

        # 4) Do the upload
        result = upload_genome(
                logger=None,
        
                shock_service_url = self.cfg.shockURL,
                handle_service_url = self.cfg.handleURL,
                workspace_service_url = self.cfg.workspaceURL,
                callback_url = self.cfg.callbackURL,

                input_directory=input_directory,
        
                workspace_name   = parsed_params['workspace_name'],
                core_genome_name = parsed_params['genome_name'],

                source           = parsed_params['source'],

                taxon_wsname     = parsed_params['taxon_wsname'],
                taxon_lookup_obj_name = parsed_params['taxon_lookup_obj_name'],
                taxon_reference = parsed_params['taxon_reference'],

                release          = parsed_params['release'],
                genetic_code     = parsed_params['genetic_code'],
                type             = parsed_params['type'],
                generate_ids_if_needed = parsed_params['generate_ids_if_needed'],

                exclude_ontologies = parsed_params['exclude_ontologies'],
                ontology_wsname = parsed_params['ontology_wsname'],
                ontology_GO_obj_name = parsed_params['ontology_GO_obj_name'],
                ontology_PO_obj_name = parsed_params['ontology_PO_obj_name'],


                provenance = ctx['provenance'],
                usermeta = parsed_params['metadata']
            )

        # 5) clear the temp directory
        shutil.rmtree(input_directory)

        # 6) return the result
        info = result['genome_info']
        details = {
            'genome_ref': str(info[6]) + '/' + str(info[0]) + '/' + str(info[4]),
            'genome_info': info,
            'report_name': result['report_name'],
            'report_ref': result['report_ref']
        }

        return details

    def refactored_import(self, ctx, params):
        # 1) validate parameters and extract defaults
        self.validate_params(params)

        # 2) construct the input directory staging area
        input_directory = self.stage_input(params)

        # 3) update default params
        self.default_params.update(params)
        params = self.default_params

        # 4) Do the upload
        files = self._find_input_files(input_directory)
        consolidated_file = self._join_files_skip_empty_lines(files)
        genome = self.parse_genbank(consolidated_file, params)
        print genome
        result = self.gi.save_one_genome({
            'workspace': params['workspace_name'],
            'name': params['genome_name'],
            'data': genome,
            "meta": params['metadata'],
        })
        ref = "{}/{}/{}".format(result['info'][6], result['info'][0],
                                result['info'][4])
        self.log("Genome saved to {}".format(ref))

        # 5) Generate Report
        reportObj = {'objects_created': [
            {'ref': ref, 'description': 'KBase Genome object'}],
                     'text_message': self.messages}
        report_info = self.report_client.create({
            'report': reportObj,
            'workspace_name': params['workspace_name']
        })

        # 6) clear the temp directory
        shutil.rmtree(input_directory)

        # 7) return the result
        info = result['info']
        details = {
            'genome_ref': ref,
            'genome_info': info,
            'report_name': report_info['name'],
            'report_ref': report_info['ref']
        }

        return details

    def validate_params(self, params):
        if 'workspace_name' not in params:
            raise ValueError('required "workspace_name" field was not defined')
        if 'genome_name' not in params:
            raise ValueError('required "genome_name" field was not defined')
        if 'file' not in params:
            raise ValueError('required "file" field was not defined')

        # one and only one of 'path', 'shock_id', or 'ftp_url' is required
        file = params['file']
        if not isinstance(file, dict):
            raise ValueError('required "file" field must be a map/dict')
        n_valid_fields = 0
        if 'path' in file and file['path'] is not None:
            n_valid_fields += 1
        if 'shock_id' in file and file['shock_id'] is not None:
            n_valid_fields += 1
        if 'ftp_url' in file and file['ftp_url'] is not None:
            n_valid_fields += 1
        if n_valid_fields < 1:
            raise ValueError('required "file" field must include one source: '
                             'path | shock_id | ftp_url')
        if n_valid_fields > 1:
            raise ValueError('required "file" field has too many sources '
                             'specified: ' + str(file.keys()))

        valid_types = ['Reference', 'User upload', 'Representative']
        if 'type' in params and params['type'] not in valid_types:
            raise ValueError('Entered value for type is not one of the valid '
                             'entries: {}'.format(", ".join(valid_types)))

    def stage_input(self, params):
        ''' Setup the input_directory by fetching the files and uncompressing if needed. '''

        # construct the input directory where we stage files
        input_directory =  os.path.join(self.cfg.sharedFolder, 'genome-upload-staging-'+str(uuid.uuid4()))
        os.makedirs(input_directory)

        # at this point, the 'file' input is validated, so we don't have to catch any special cases
        # we expect one and only one of path, shock_id, or ftp_url

        # determine how to get the file: if it is from shock, download it.  If it
        # is just sitting there, then use it.  Move the file to the staging input directory
        file = params['file']
        genbank_file_path = None
        if 'path' in file and file['path'] is not None:
            # copy the local file to the input staging directory
            # (NOTE: could just move it, but then this method would have the side effect of moving your
            # file which another SDK module might have an open handle on)
            local_file_path = file['path']
            genbank_file_path = os.path.join(input_directory, os.path.basename(local_file_path))
            shutil.copy2(local_file_path, genbank_file_path)

        if 'shock_id' in file and file['shock_id'] is not None:
            # handle shock file
            print('Downloading file from SHOCK node: {} - {}'.format(
                self.cfg.shockURL, file['shock_id']))
            sys.stdout.flush()
            file_name = self.dfu.shock_to_file({
                                    'file_path': input_directory,
                                    'shock_id': file['shock_id']
                                })['node_file_name']
            genbank_file_path = os.path.join(input_directory, file_name)

        if 'ftp_url' in file and file['ftp_url'] is not None:
            print('Downloading file from: ' + str(file['ftp_url']))
            local_file_path = self.dfu.download_web_file({
                'file_url': file['ftp_url'],
                'download_type': 'FTP'
            })['copy_file_path']
            genbank_file_path = os.path.join(input_directory,
                                             os.path.basename(local_file_path))
            shutil.copy2(local_file_path, genbank_file_path)

        # extract the file if it is compressed
        if genbank_file_path is not None:
            print("staged input file =" + genbank_file_path)
            self.dfu.unpack_file({'file_path': genbank_file_path})

        else:
            raise ValueError('No valid files could be extracted based on the input')

        return input_directory

    def parse_genbank(self, file_path, params):
        print("Saving original file to shock")
        shock_res = self.dfu.file_to_shock({
            'file_path': file_path,
            'make_handle': 1,
            'pack': 'gzip',
        })
        # Write and save assembly file
        contigs = Bio.SeqIO.parse(file_path, "genbank")
        assembly_ref = self._save_assembly(contigs, params)
        assembly_data = self.dfu.get_objects(
            {'object_refs': [assembly_ref],
             'ignore_errors': 0})['data'][0]['data']
        genome = {
            "id": params['genome_name'],
            "source": params['source'],
            "type": params['type'],
            "original_source_file_name": os.path.basename(file_path),
            "genetic_code": params['genetic_code'],
            "assembly_ref": assembly_ref,
            "gc_content": assembly_data['gc_content'],
            "dna_size": assembly_data['dna_size'],
            "md5": assembly_data['md5'],
            "genbank_handle_ref": shock_res['handle']['hid'],
            "publications": set(),
            "domain": "Unknown",
            "ontology_events": [{
                "method": "GenomeFileUtils Genbank uploader from annotations",
                "method_version": self.version,
                "timestamp": self.time_string
            }],
            "contig_ids": [],
            "contig_lengths": [],
            "features": [],
            "non_coding_features": [],
            "cdss": [],
            'mrnas': [],
        }
        dates = []
        # Parse data from genbank file
        contigs = Bio.SeqIO.parse(file_path, "genbank")
        for record in contigs:
            if 'date' in record.annotations:
                dates.append(time.strptime(record.annotations['date'],
                                           "%d-%b-%Y"))
            genome['contig_ids'].append(record.id)
            genome['contig_lengths'].append(len(record))
            for k, v in self._parse_features(record).items():
                genome[k].extend(v)
            genome["publications"] |= self._get_pubs(record)

            if "source_id" in genome:
                continue  # only do the following once
            genome["source_id"] = record.id.split('.')[0]
            genome['scientific_name'] = record.annotations.get('organism',
                                                               'unknown_taxon')
            tax_info = self.gi.retrieve_taxon(params['taxon_wsname'],
                                              genome['scientific_name'])
            if tax_info:
                genome['taxonomy'], genome['taxon_ref'], genome['domain'] = tax_info
            genome['notes'] = record.annotations.get('comment', "")
        genome['num_contigs'] = len(genome['contig_ids'])
        dates.sort()
        if dates:
            genome['external_source_origination_date'] = time.strftime(
                "%d-%b-%Y",dates[0])
            if dates[0] != dates[-1]:
                genome['external_source_origination_date'] += " _ " + \
                    time.strftime("%d-%b-%Y", dates[-1])
        genome['ontology_present'] = dict(self.ontologies_present)
        # can't serialize a set
        genome['publications'] = list(genome['publications'])
        return genome

    def _save_assembly(self, contigs, params):
        print("Saving sequence as Assembly object")
        assembly_id = "{}_assembly".format(params['genome_name'])
        fasta_file = "{}/{}_assembly.fasta".format(
            self.cfg.sharedFolder, params['genome_name'], self.time_string)
        Bio.SeqIO.write(contigs, fasta_file, "fasta")
        assembly_ref = self.aUtil.save_assembly_from_fasta(
            {'file': {'path': fasta_file},
             'workspace_name': params['workspace_name'],
             'assembly_name': assembly_id})
        self.log("Assembly saved to {}".format(assembly_ref))
        return assembly_ref

    def _find_input_files(self, input_directory):
        print("Scanning for Genbank Format files.")
        valid_extensions = [".gbff", ".gbk", ".gb", ".genbank", ".dat", ".gbf"]

        files = os.listdir(os.path.abspath(input_directory))
        self.log("Genbank Files : " + ", ".join(files))
        genbank_files = [x for x in files if
                         os.path.splitext(x)[-1] in valid_extensions]

        if len(genbank_files) == 0:
            raise Exception(
                "The input directory does not have any files with one of the "
                "following extensions %s." % (",".join(valid_extensions)))

        self.log("Found {} genbank files".format(len(genbank_files)))

        input_files = []
        for genbank_file in genbank_files:
            input_files.append(os.path.join(input_directory, genbank_file))

        return input_files

    def _join_files_skip_empty_lines(self, input_files):
            """ Applies strip to each line of each input file.
            Args:
                input_files: Paths to input files in Genbank format.
            Returns:
                Path to resulting file (currenly it's the same file as input).
            """
            if len(input_files) == 0:
                raise ValueError("NO GENBANK FILE")
            temp_dir = os.path.join(os.path.dirname(input_files[0]), "combined")
            if not os.path.exists(temp_dir):
                os.makedirs(temp_dir)
            ret_file = os.path.join(temp_dir, os.path.basename(input_files[0]))

            # take in Genbank file and remove all empty lines from it.
            with open(ret_file, 'w', buffering=2 ** 20) as f_out:
                for input_file in input_files:
                    with open(input_file, 'r') as f_in:
                        for line in f_in:
                            line = line.rstrip('\r\n')
                            if line.strip():
                                f_out.write(line + '\n')
            return ret_file

    def _get_pubs(self, record):
        pub_list = []
        for in_pub in record.annotations.get('references', []):
            out_pub = [
                0,  # pmid
                "",  # source
                in_pub.title,
                "",  # web address
                "",  # date
                in_pub.authors,
                in_pub.journal,
            ]
            date_match = re.match("\((\d{4})\)", in_pub.journal)
            if date_match:
                out_pub[4] = date_match.group(1)
            if in_pub.pubmed_id:
                out_pub[0:4] = [
                    int(in_pub.pubmed_id),
                    "PubMed",
                    in_pub.title,
                    "http://www.ncbi.nlm.nih.gov/pubmed/{}".format(
                        in_pub.pubmed_id)]
            pub_list.append(tuple(out_pub))
        self.log("Parsed {} publication records".format(len(pub_list)))
        return set(pub_list)

    def _get_ontology(self, feature):
        ontology = defaultdict(dict)
        for key in ("GO_process", "GO_function", "GO_component"):
            if key in feature.qualifiers:
                sp = feature.qualifiers[key][0].split(" - ")
                ontology['GO'][sp[0]] = sp[1]
                self.ontologies_present['GO'][sp[0]] = sp[1]
        return dict(ontology)

    def _parse_features(self, record):
        def _location(seq, feat):
            strand_trans = ("", "+", "-")
            loc = []
            for part in feat.location.parts:
                if part.strand >= 0:
                    begin = int(part.start)
                else:
                    begin = int(part.end)
                loc.append((
                        record.id,
                        begin + 1,
                        strand_trans[part.strand],
                        len(part)))
            return loc

        def _aliases(feat):
            keys = ('locus_tag', 'old_locus_tag', 'db_xref', 'protein_id')
            aliases = []
            for k, v in feat.qualifiers.items():
                if k in keys:
                    aliases.extend(v)
            return aliases

        def _is_parent(feat1, feat2):
            def _contains(loc1, loc2):
                if loc1[2] != loc2[2]:  # different strands
                    return False
                elif loc1[2] == "+":
                    return loc2[1] >= loc1[1] and (
                    loc2[1] + loc2[3] <= loc1[1] + loc1[3])
                else:
                    return loc2[1] <= loc1[1] and (
                    loc2[1] - loc2[3] >= loc1[1] - loc1[3])

            # if any location in the CDS in contained by any location of the gene
            return any([_contains(x, y) for x, y in
                        itertools.product(feat1['location'], feat2['location'])])

        def _propagate_cds_props_to_gene(cds, gene):
            # Check gene function
            if "function" not in gene or gene["function"] is None or len(
                    gene["function"]) == 0:
                gene["function"] = cds.get("function", "")
            # Put longest protein_translation to gene
            if "protein_translation" not in gene or (
                len(gene["protein_translation"]) <
                    len(cds["protein_translation"])):
                gene["protein_translation"] = cds["protein_translation"]
                gene["protein_translation_length"] = len(
                    cds["protein_translation"])
            # Merge cds["aliases"] -> gene["aliases"]
            gene["aliases"] = list(set(cds['aliases']+gene['aliases']))
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

        skiped_features = Counter()
        noncoding_types = Counter()
        excluded_features = ('source', 'exon')
        genes, cdss, mrnas, noncoding = {}, {}, {}, []
        for in_feature in record.features:
            if in_feature.type in excluded_features:
                skiped_features[in_feature.type] += 1
                continue
            feat_seq = in_feature.extract(record)
            _id = in_feature.qualifiers.get("locus_tag", [""])[0]
            if not _id:
                _id = in_feature.qualifiers.get("gene", [""])[0]
            out_feature = {
                "id": "_".join([in_feature.type, _id]),
                "location": _location(feat_seq, in_feature),
                "dna_sequence": str(feat_seq.seq),
                "dna_sequence_length": len(feat_seq),
                "md5": hashlib.md5(str(feat_seq)).hexdigest(),
                "function": in_feature.qualifiers.get("product", [""])[0],
                "ontology_terms": {},  # self._get_ontology(in_feature),
                "note": in_feature.qualifiers.get("note", ""),
                'aliases': _aliases(in_feature),
            }

            if in_feature.type == 'CDS':
                out_feature.update({
                    "protein_translation": in_feature.qualifiers.get(
                        "translation", [""])[0],
                    "protein_md5": hashlib.md5(str(feat_seq)).hexdigest(),
                    "parent_gene": "",
                    "parent_mrna": "",
                })
                out_feature["protein_md5"] = hashlib.md5(
                    str(out_feature['protein_translation'])).hexdigest(),
                out_feature['protein_translation_length'] = len(
                    out_feature['protein_translation'])
                if _id in genes:
                    out_feature['id'] += "_" + str(len(genes[_id]['cdss'])+1)
                    genes[_id]['cdss'].append(out_feature['id'])
                    _propagate_cds_props_to_gene(out_feature, genes[_id])
                    out_feature['parent_gene'] = _id
                    if not _is_parent(genes[_id], out_feature):
                        self.log("{} is annotated as the parent gene of {} "
                                 "but coordinates do not match".format(
                                    _id, out_feature['id']))

                mrna_id = "mRNA" + out_feature['id'][3:]
                if mrna_id in mrnas:
                    if not _is_parent(mrnas[mrna_id], out_feature):
                        self.log("{} is annotated as the parent transcript of "
                                 "{} but coordinates do not match".format(
                                    mrna_id, out_feature['id']))
                    mrnas[mrna_id]['cds'] = out_feature['id']
                    out_feature['parent_mrna'] = mrna_id
                cdss[out_feature['id']] = out_feature

            elif in_feature.type == 'gene':
                out_feature.update({
                    "id": _id,
                    "type": 'gene',
                    "protein_translation": "",
                    "protein_translation_length": 0,
                    "mrnas": [],
                    'cdss': [],
                })
                genes[_id] = out_feature

            elif in_feature.type == 'mRNA':
                out_feature.update({
                    "parent_gene": "",
                    "cds": "",
                })
                # not in this feature type
                del out_feature['function'], out_feature['ontology_terms']
                if _id in genes:
                    out_feature['id'] += "_" + str(len(genes[_id]['mrnas'])+1)
                    genes[_id]['mrnas'].append(out_feature['id'])
                    out_feature['parent_gene'] = _id
                    if not _is_parent(genes[_id], out_feature):
                        self.log("{} is annotated as the parent gene of {} "
                                 "but coordinates do not match".format(
                            _id, out_feature['id']))
                mrnas[out_feature['id']] = out_feature

            else:
                noncoding_types[in_feature.type] += 1
                out_feature["type"] = in_feature.type
                out_feature['id'] += "_" + str(noncoding_types[in_feature.type])
                noncoding.append(out_feature)

        coding = []
        for g in genes.values():
            if len(g['cdss']):
                coding.append(g)
            else:
                del g['protein_translation'], g['protein_translation_length'],\
                    g['mrnas'], g['cdss']
                noncoding.append(g)
        self.log("Features skipped\n{}\n".format("\n".join([
            "{}: {}".format(k, v) for k, v in skiped_features.items()])))
        self.log("Noncoding Features\n{}\n".format("\n".join([
            "{}: {}".format(k, v) for k, v in noncoding_types.items()])))

        return {'features': coding, 'non_coding_features': noncoding,
                'cdss': cdss.values(), 'mrnas': mrnas.values()}
