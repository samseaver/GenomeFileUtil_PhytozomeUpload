
import copy
import datetime
import hashlib
import json
import os
import re
import shutil
import sys
import time
import uuid
from collections import Counter, defaultdict, OrderedDict

import Bio.SeqIO
import Bio.SeqUtils
from Bio import Seq
from Bio.Data.CodonTable import TranslationError
from Bio.SeqFeature import ExactPosition

from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.DataFileUtilClient import DataFileUtil
from GenomeFileUtil.core.GenomeInterface import GenomeInterface
from GenomeFileUtil.core.GenomeUtils import is_parent, propagate_cds_props_to_gene, warnings
from GenomeFileUtil.core.GenomeUtils import parse_inferences, load_ontology_mappings
from installed_clients.WorkspaceClient import Workspace

MAX_MISC_FEATURE_SIZE = 10000
MAX_PARENT_LOOKUPS = 5


class GenbankToGenome:
    def __init__(self, config):
        self.cfg = config
        self.gi = GenomeInterface(config)
        self.dfu = DataFileUtil(config.callbackURL)
        self.aUtil = AssemblyUtil(config.callbackURL)
        self.ws = Workspace(config.workspaceURL)
        self._messages = []
        self.time_string = str(datetime.datetime.fromtimestamp(
            time.time()).strftime('%Y_%m_%d_%H_%M_%S'))
        yml_text = open('/kb/module/kbase.yml').read()
        self.version = re.search("module-version:\n\W+(.+)\n", yml_text).group(1)
        self.generate_parents = False
        self.generate_ids = False
        self.genes = OrderedDict()
        self.mrnas = OrderedDict()
        self.cdss = OrderedDict()
        self.noncoding = []
        self.ontologies_present = defaultdict(dict)
        self.ontology_events = list()
        self.skiped_features = Counter()
        self.feature_counts = Counter()
        self.orphan_types = Counter()
        self.contig_seq = {}
        self.circ_contigs = set()
        self.features_spaning_zero = set()
        self.genome_warnings = []
        self.genome_suspect = False
        self.defects = Counter()
        self.spoofed_genes = 0
        self.excluded_features = ('source', 'exon')
        self.ont_mappings = load_ontology_mappings('/kb/module/data')
        self.code_table = 11
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
            'metadata': {}
        }

    def log(self, message):
        self._messages.append(message)
        print('{0:.2f}'.format(time.time()) + ': ' + str(message))

    @property
    def messages(self):
        return "\n".join(self._messages)

    def refactored_import(self, ctx, params):
        # 1) validate parameters and extract defaults
        self.validate_params(params)

        # 2) construct the input directory staging area
        input_directory = self.stage_input(params)

        # 3) update default params
        self.default_params.update(params)
        params = self.default_params
        self.generate_parents = params.get('generate_missing_genes')
        self.generate_ids = params.get('generate_ids_if_needed')
        if params.get('genetic_code'):
            self.code_table = params['genetic_code']

        # 4) Do the upload
        files = self._find_input_files(input_directory)
        consolidated_file = self._join_files_skip_empty_lines(files)
        genome = self.parse_genbank(consolidated_file, params)
        if params.get('genetic_code'):
            genome["genetic_code"] = params['genetic_code']
        ###
        # DEBUGGING CODE INSTRUCTIONS TO BE KEPT FOR DETERMINING ISSUE WITH A PARTICULAR FILE
        # THAT FAILS TYPESPEC CHECKING - THIS ALLOW YOU TO LOOK AT THE JSON BEFORE SAVED:
        # Turn this on :
        #  1) uncomment the TWO LINES for json printing lines below the ###
        #  2) move skip/utility_test/problem_genome_test.py into the test dir 
        #  3) change the file location in the test_problem_genome_for_json test
        #  4) add your test file to the test/data dir
        #  5) run kb-sdk test as normal
        #  6) after test completes go and find the json file at test_local/workdir/tmp/ProblemGenome.json
        ###
        #with open(self.cfg.sharedFolder+'/ProblemGenome.json', 'w') as outfile:
        #    json.dump(genome, outfile, indent=4)
        #    json.dump(genome, outfile)

        result = self.gi.save_one_genome({
            'workspace': params['workspace_name'],
            'name': params['genome_name'],
            'data': genome,
            "meta": params['metadata'],
        })
        ref = "{}/{}/{}".format(result['info'][6], result['info'][0],
                                result['info'][4])
        self.log("Genome saved to {}".format(ref))

        # 5) clear the temp directory
        shutil.rmtree(input_directory)

        # 6) return the result
        info = result['info']
        details = {
            'genome_ref': ref,
            'genome_info': info
        }

        return details

    @staticmethod
    def validate_params(params):
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
                             'specified: ' + str(list(file.keys())))
        if params.get('genetic_code'):
            if not (isinstance(params['genetic_code'], int) and 0 < params['genetic_code'] < 32):
                raise ValueError("Invalid genetic code specified: {}".format(params))

    def stage_input(self, params):
        """ Setup the input_directory by fetching the files and uncompressing if needed. """

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
            self.log('Downloading file from SHOCK node: {} - {}'.format(
                self.cfg.shockURL, file['shock_id']))
            sys.stdout.flush()
            file_name = self.dfu.shock_to_file({
                                    'file_path': input_directory,
                                    'shock_id': file['shock_id']
                                })['node_file_name']
            genbank_file_path = os.path.join(input_directory, file_name)

        if 'ftp_url' in file and file['ftp_url'] is not None:
            self.log('Downloading file from: ' + str(file['ftp_url']))
            local_file_path = self.dfu.download_web_file({
                'file_url': file['ftp_url'],
                'download_type': 'FTP'
            })['copy_file_path']
            genbank_file_path = os.path.join(input_directory,
                                             os.path.basename(local_file_path))
            shutil.copy2(local_file_path, genbank_file_path)

        # extract the file if it is compressed
        if genbank_file_path is not None:
            self.log("staged input file =" + genbank_file_path)
            self.dfu.unpack_file({'file_path': genbank_file_path})

        else:
            raise ValueError('No valid files could be extracted based on the input')

        return input_directory

    def parse_genbank(self, file_path, params):
        self.log("Saving original file to shock")
        shock_res = self.dfu.file_to_shock({
            'file_path': file_path,
            'make_handle': 1,
            'pack': 'gzip',
        })
        # Write and save assembly file
        assembly_ref = self._save_assembly(file_path, params)
        assembly_data = self.dfu.get_objects(
            {'object_refs': [assembly_ref],
             'ignore_errors': 0})['data'][0]['data']
        genome = {
            "id": params['genome_name'],
            "original_source_file_name": os.path.basename(file_path),
            "assembly_ref": assembly_ref,
            "gc_content": assembly_data['gc_content'],
            "dna_size": assembly_data['dna_size'],
            "md5": assembly_data['md5'],
            "genbank_handle_ref": shock_res['handle']['hid'],
            "publications": set(),
            "contig_ids": [],
            "contig_lengths": [],
        }
        genome['source'], genome['genome_tiers'] = \
            self.gi.determine_tier(params['source'])

        if params.get('genome_type'):
            genome['genome_type'] = params['genome_type']

        dates = []
        # Parse data from genbank file
        contigs = Bio.SeqIO.parse(file_path, "genbank")
        for record in contigs:
            r_annot = record.annotations
            self.log("parsing contig: " + record.id)
            if 'date' in r_annot:
                dates.append(time.strptime(r_annot['date'], "%d-%b-%Y"))
            genome['contig_ids'].append(record.id)
            genome['contig_lengths'].append(len(record))
            genome["publications"] |= self._get_pubs(r_annot)
            organism = r_annot.get('organism', 'Unknown Organism')
            if 'scientific_name' not in genome:
                genome['scientific_name'] = organism
            elif genome['scientific_name'] != organism:
                warn = "Multiple organism in provided files: {}, {}".format(
                    genome['scientific_name'], organism)
                genome['warnings'] = genome.get('warnings', []) + [warn]

            # only do the following once(on the first contig)
            if "source_id" not in genome:
                genome["source_id"] = record.id.split('.')[0]
                genome['taxonomy'], genome['taxon_ref'], genome['domain'], \
                genome['genetic_code'] = self.gi.retrieve_taxon(
                    params['taxon_wsname'], genome['scientific_name'])
                self.code_table = genome['genetic_code']
                genome["molecule_type"] = r_annot.get('molecule_type', 'DNA')
                genome['notes'] = r_annot.get('comment', "").replace('\\n', '\n')

            self._parse_features(record, params['source'])

        genome.update(self.get_feature_lists())

        genome['num_contigs'] = len(genome['contig_ids'])
        # add dates
        dates.sort()
        if dates:
            genome['external_source_origination_date'] = time.strftime(
                "%d-%b-%Y", dates[0])
            if dates[0] != dates[-1]:
                genome['external_source_origination_date'] += " _ " + \
                    time.strftime("%d-%b-%Y", dates[-1])

        if self.ontologies_present:
            genome['ontologies_present'] = dict(self.ontologies_present)
            genome["ontology_events"] = self.ontology_events
        genome['feature_counts'] = dict(self.feature_counts)
        # can't serialize a set
        genome['publications'] = list(genome['publications'])

        if len(genome['cdss']) and (self.defects['cds_seq_not_matching'] /
                                    float(len(genome['cdss'])) > 0.02):
            self.genome_warnings.append(warnings["genome_inc_translation"].format(
                self.defects['cds_seq_not_matching'], len(genome['cdss'])))
            self.genome_suspect = 1

        if self.defects['bad_parent_loc']:
            self.genome_warnings.append("There were {} parent/child "
                "relationships that were not able to be determined. Some of "
                "these may have splice variants that may be valid "
                "relationships.".format(self.defects['bad_parent_loc']))

        if self.defects['spoofed_genes']:
            self.genome_warnings.append(warnings['spoofed_genome'].format(
                self.defects['spoofed_genes']))
            genome['suspect'] = 1

        if self.defects['not_trans_spliced']:
            self.genome_warnings.append(warnings['genome_not_trans_spliced']
                                        .format(self.defects['not_trans_spliced']))
            genome['suspect'] = 1

        if self.genome_warnings:
            genome['warnings'] = self.genome_warnings
        if self.genome_suspect:
            genome['suspect'] = 1
        self.log("Feature Counts: {}".format(genome['feature_counts']))
        return genome

    def _save_assembly(self, genbank_file, params):
        """Convert genbank file to fasta and sve as assembly"""
        contigs = Bio.SeqIO.parse(genbank_file, "genbank")
        assembly_id = "{}_assembly".format(params['genome_name'])
        fasta_file = "{}/{}_assembly.fasta".format(
            self.cfg.sharedFolder, params['genome_name'], self.time_string)

        out_contigs = []
        extra_info = defaultdict(dict)
        for in_contig in contigs:
            if in_contig.annotations.get('topology', "") == 'circular':
                extra_info[in_contig.id]['is_circ'] = 1
                self.circ_contigs.add(in_contig.id)
            elif in_contig.annotations.get('topology', "") == 'linear':
                extra_info[in_contig.id]['is_circ'] = 0
            out_contigs.append(in_contig)
            self.contig_seq[in_contig.id] = in_contig.seq.upper()

        assembly_ref = params.get("use_existing_assembly")
        if assembly_ref:
            if not re.match("\d+\/\d+\/\d+", assembly_ref):
                raise ValueError("Assembly ref: {} is not a valid format. Must"
                                 " be in numerical <ws>/<object>/<version>"
                                 " format.".format(assembly_ref))
            ret = self.dfu.get_objects(
                {'object_refs': [assembly_ref]}
            )['data'][0]
            if "KBaseGenomeAnnotations.Assembly" not in ret['info'][2]:
                raise ValueError("{} is not a reference to an assembly"
                                 .format(assembly_ref))
            unmatched_ids = list()
            unmatched_ids_md5s = list()
            for current_contig in self.contig_seq.keys():
                current_contig_md5 = hashlib.md5(
                    str(self.contig_seq[current_contig]).encode('utf8')
                ).hexdigest()
                if current_contig in ret['data']['contigs']:
                    if current_contig_md5 != ret['data']['contigs'][current_contig]['md5']:
                        unmatched_ids_md5s.append(current_contig)
                else:
                    unmatched_ids.append(current_contig)
            if len(unmatched_ids) > 0:
                raise ValueError(warnings['assembly_ref_extra_contigs'].format(", ".join(unmatched_ids)))
            if len(unmatched_ids_md5s) > 0:
                raise ValueError(warnings["assembly_ref_diff_seq"].format(", ".join(unmatched_ids_md5s)))                
            self.log("Using supplied assembly: {}".format(assembly_ref))
            return assembly_ref
        self.log("Saving sequence as Assembly object")
        Bio.SeqIO.write(out_contigs, fasta_file, "fasta")
        assembly_ref = self.aUtil.save_assembly_from_fasta(
            {'file': {'path': fasta_file},
             'workspace_name': params['workspace_name'],
             'assembly_name': assembly_id,
             'contig_info': extra_info})
        self.log("Assembly saved to {}".format(assembly_ref))
        return assembly_ref

    def _find_input_files(self, input_directory):
        self.log("Scanning for Genbank Format files.")
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

    def _get_pubs(self, r_annotations):
        """Get a contig's publications"""
        pub_list = []
        for in_pub in r_annotations.get('references', []):
            # don't add blank pubs
            if not in_pub.authors:
                continue
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

    def _parse_features(self, record, source):
        def _get_id(feat, tags=None):
            """Assign a id to a feature based on the first tag that exists"""
            _id = ""
            if not tags:
                tags = ['locus_tag', 'kbase_id']
            for t in tags:
                _id = feat.qualifiers.get(t, [""])[0]
                if _id:
                    break

            if not _id:
                if feat.type == 'gene':
                    if not self.generate_ids:
                        raise ValueError("Unable to find a valid id for genes "
                                         "among these tags: {}. Correct the "
                                         "file or rerun with generate_ids"
                                         .format(", ".join(tags)))
                    self.orphan_types['gene'] += 1
                    _id = "gene_{}".format(self.orphan_types['gene'])
                if 'rna' in feat.type.lower() or feat.type in {'CDS', 'sig_peptide',
                                                               'five_prime_UTR', 'three_prime_UTR'}:
                    _id = "gene_{}".format(self.orphan_types['gene'])

            return _id

        def _location(feat):
            """Convert to KBase style location objects"""
            strand_trans = ("", "+", "-")
            loc = []
            for part in feat.location.parts:
                contig_id = part.ref if part.ref else record.id
                if part.strand >= 0:
                    begin = int(part.start) + 1
                else:
                    begin = int(part.end)
                loc.append((
                        contig_id,
                        begin,
                        strand_trans[part.strand],
                        len(part)))
            return loc

        def _warn(message):
            if message not in out_feat.get('warnings', []):
                out_feat['warnings'] = out_feat.get('warnings', []) + [message]

        def _check_suspect_location(parent=None):
            if 'trans_splicing' in out_feat.get('flags', []):
                return

            if out_feat['location'] == sorted(out_feat['location'],
                    reverse=(in_feature.location.strand == -1)):
                return

            if record.id in self.circ_contigs and \
                    in_feature.location.start == 0 \
                    and in_feature.location.end == len(record):
                self.features_spaning_zero.add(out_feat['id'])
                return

            if parent and parent['id'] in self.features_spaning_zero:
                return

            _warn(warnings['not_trans_spliced'])
            self.defects['not_trans_spliced'] += 1

        for in_feature in record.features:
            if in_feature.type in self.excluded_features:
                self.skiped_features[in_feature.type] += 1
                continue
            feat_seq = self._get_seq(in_feature, record.id)
            if source == "Ensembl":
                _id = _get_id(in_feature, ['gene', 'locus_tag'])
            else:
                _id = _get_id(in_feature)
            # The following is common to all the feature types
            out_feat = {
                "id": "_".join([_id, in_feature.type]),
                "location": _location(in_feature),
                "dna_sequence": str(feat_seq),
                "dna_sequence_length": len(feat_seq),
                "md5": hashlib.md5(str(feat_seq).encode('utf8')).hexdigest(),
            }
            if not _id:
                out_feat['id'] = in_feature.type

            # validate input feature
            # note that end is the larger number regardless of strand
            if int(in_feature.location.end) > len(record):
                self.genome_warnings.append(
                    warnings["coordinates_off_end"].format(out_feat['id']))
                self.genome_suspect = 1
                continue

            for piece in in_feature.location.parts:
                if not isinstance(piece.start, ExactPosition) \
                        or not isinstance(piece.end, ExactPosition):
                    _warn(warnings["non_exact_coordinates"])

            self.feature_counts[in_feature.type] += 1

            # add optional fields
            if 'note' in in_feature.qualifiers:
                out_feat['note'] = in_feature.qualifiers["note"][0]

            out_feat.update(self._get_aliases_flags_functions(in_feature))

            ont, db_xrefs = self._get_ontology_db_xrefs(in_feature)
            if ont:
                out_feat['ontology_terms'] = ont
            if db_xrefs:
                out_feat['db_xrefs'] = db_xrefs

            if 'inference' in in_feature.qualifiers:
                out_feat['inference_data'] = parse_inferences(
                    in_feature.qualifiers['inference'])

            _check_suspect_location(self.genes.get(_id))

            # add type specific features
            if in_feature.type == 'CDS':
                self.process_cds(_id, feat_seq, in_feature, out_feat)

            elif in_feature.type == 'gene':
                self.process_gene(_id, out_feat)

            elif in_feature.type == 'mRNA':
                self.process_mrna(_id, out_feat)

            else:
                self.noncoding.append(self.process_noncoding(_id, in_feature.type, out_feat))

    def get_feature_lists(self):
        """sort genes into their final arrays"""
        coding = []
        for g in self.genes.values():
            if len(g['cdss']):
                if g['mrnas'] and len(g['mrnas']) != len(g['cdss']):
                    msg = "The length of the mrna and cdss arrays are not equal"
                    g['warnings'] = g.get('warnings', []) + [msg]

                # remove duplicates that may arise from CDS info propagation
                for key in ('functions', 'aliases', 'db_xrefs'):
                    if key in g:
                        g[key] = list(set(g[key]))
                if not g['mrnas']:
                    del g['mrnas']
                del g['type']
                coding.append(g)
                self.feature_counts["protein_encoding_gene"] += 1
            else:
                del g['mrnas'], g['cdss']
                self.noncoding.append(g)
                self.feature_counts["non_coding_features"] += 1
        return {'features': coding, 'non_coding_features': self.noncoding,
                'cdss': list(self.cdss.values()), 'mrnas': list(self.mrnas.values())}

    def _get_seq(self, feat, contig):
        """Extract the DNA sequence for a feature"""
        seq = []
        for part in feat.location.parts:
            strand = part.strand
            # handle trans-splicing across contigs
            if part.ref:
                part_contig = part.ref
            else:
                part_contig = contig

            if strand >= 0:
                seq.append(str(self.contig_seq[part_contig]
                               [part.start:part.end]))
            else:
                seq.append(str(self.contig_seq[part_contig]
                               [part.start:part.end].reverse_complement()))
        return "".join(seq)

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
        ontology = defaultdict(dict)
        db_xrefs = []
        for key in ("GO_process", "GO_function", "GO_component"):
            ontology_event_index = self._create_ontology_event("GO")
            for term in feature.qualifiers.get(key, []):
                sp = term.split(" - ")
                ontology['GO'][sp[0]] = [ontology_event_index]
                self.ontologies_present['GO'][sp[0]] = self.ont_mappings['GO'].get(sp[0], '')

        for ref in feature.qualifiers.get('db_xref', []):
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

        return dict(ontology), sorted(db_xrefs)

    @staticmethod
    def _get_aliases_flags_functions(feat):
        """Get the values for aliases flags and features from qualifiers"""
        alias_keys = {'locus_tag', 'old_locus_tag', 'protein_id',
                      'transcript_id', 'gene', 'EC_number', 'gene_synonym'}
        result = defaultdict(list)
        for key, val_list in feat.qualifiers.items():
            if key in alias_keys:
                result['aliases'].extend([(key, val) for val in val_list])
            # flags have no other information associated with them
            if val_list == ['']:
                result['flags'].append(key)
            if key == 'function':
                result['functional_descriptions'].extend(val_list[0].split('; '))
            if key == 'product':
                result['functions'] = val_list

        return result

    def _find_parent_gene(self, potential_id, feature):
        if potential_id in self.genes:
            lookup_attempts = 0
            while lookup_attempts < MAX_PARENT_LOOKUPS:
                if is_parent(self.genes[potential_id], feature):
                    return potential_id

                lookup_attempts += 1
                try:
                    potential_id = list(self.genes.keys())[-(lookup_attempts + 1)]
                except IndexError:
                    break  # no more genes that could match exist

            self.defects['bad_parent_loc'] += 1
        return None

    def process_gene(self, _id, out_feat):
        out_feat.update({
            "id": _id,
            "type": 'gene',
            "mrnas": [],
            'cdss': [],
        })
        if _id in self.genes:
            raise ValueError("Duplicate gene ID: {}".format(_id))
        self.genes[_id] = out_feat

    def process_noncoding(self, gene_id, feat_type, out_feat):
        out_feat["type"] = feat_type

        # this prevents big misc_features from blowing up the genome size
        if out_feat['dna_sequence_length'] > MAX_MISC_FEATURE_SIZE:
            del out_feat['dna_sequence']

        gene_id = self._find_parent_gene(gene_id, out_feat)
        if gene_id:
            if 'children' not in self.genes[gene_id]:
                self.genes[gene_id]['children'] = []
            out_feat['id'] += "_" + str(len(self.genes[gene_id]['children']) + 1)
            self.genes[gene_id]['children'].append(out_feat['id'])
            out_feat['parent_gene'] = gene_id
        else:
            self.orphan_types[feat_type] += 1
            out_feat['id'] += "_" + str(self.orphan_types[feat_type])

        return out_feat

    def process_mrna(self, gene_id, out_feat):
        if gene_id not in self.genes and self.generate_parents:
                self.process_gene(gene_id, copy.copy(out_feat))

        gene_id = self._find_parent_gene(gene_id, out_feat)
        if gene_id:
            out_feat['id'] = "_".join((gene_id, "mRNA", str(len(self.genes[gene_id]['mrnas']) + 1)))
            self.genes[gene_id]['mrnas'].append(out_feat['id'])
            out_feat['parent_gene'] = gene_id
        else:
            self.orphan_types['mrna'] += 1
            out_feat['id'] = "mRNA_{}".format(self.orphan_types['mrna'])
            out_feat['warnings'] = out_feat.get('warnings', []) + [
                'Unable to find parent gene for ' + str(out_feat['id'])]

        self.mrnas[out_feat['id']] = out_feat

    def process_cds(self, gene_id, feat_seq, in_feature, out_feat):
        # Associate CDS with parents
        if gene_id not in self.genes:
            if not self.generate_parents:
                self.log("Expected gene id: {}".format(gene_id))
                raise ValueError(warnings['no_spoof'])
            new_feat = copy.copy(out_feat)
            new_feat['id'] = gene_id
            new_feat['warnings'] = [warnings['spoofed_gene']]
            self.feature_counts['gene'] += 1
            self.defects['spoofed_genes'] += 1
            self.process_gene(gene_id, new_feat)

        gene_id = self._find_parent_gene(gene_id, out_feat)
        if gene_id:
            out_feat['id'] = "_".join((gene_id, "CDS", str(len(self.genes[gene_id]['cdss']) + 1)))
            self.genes[gene_id]['cdss'].append(out_feat['id'])
            out_feat['parent_gene'] = gene_id
        else:
            self.orphan_types['cds'] += 1
            out_feat['id'] = "CDS_{}".format(self.orphan_types['cds'])
            out_feat['warnings'] = out_feat.get('warnings', []) + [
                'Unable to find parent gene for ' + str(out_feat['id'])]

        # there is a 1 to 1 relationship of mRNA to CDS so XXX_mRNA_1 will match XXX_CDS_1
        mrna_id = out_feat["id"].replace('CDS', 'mRNA')
        if mrna_id in self.mrnas:
            if not is_parent(self.mrnas[mrna_id], out_feat):
                out_feat['warnings'] = out_feat.get('warnings', []) + [
                    warnings['cds_mrna_cds'].format(mrna_id)]
                self.mrnas[mrna_id]['warnings'] = self.mrnas[mrna_id].get(
                    'warnings', []) + [warnings['cds_mrna_mrna']]
                self.defects['bad_parent_loc'] += 1
            else:
                out_feat['parent_mrna'] = mrna_id
                self.mrnas[mrna_id]['cds'] = out_feat['id']

        # process protein
        prot_seq = in_feature.qualifiers.get("translation", [""])[0]

        # allow a little slack to account for frameshift and stop codon
        if prot_seq and abs(len(prot_seq) * 3 - len(feat_seq)) > 4:
            out_feat['warnings'] = out_feat.get('warnings', []) + [
                warnings["inconsistent_CDS_length"].format(len(feat_seq),
                                                           len(prot_seq))]
            self.genome_warnings.append(
                warnings['genome_inc_CDS_length'].format(
                    out_feat['id'], len(feat_seq), len(prot_seq)))
            self.genome_suspect = 1

        try:
            if prot_seq and prot_seq != Seq.translate(
                    feat_seq, self.code_table, cds=True).strip("*"):
                out_feat['warnings'] = out_feat.get('warnings', []) + [
                    warnings["inconsistent_translation"]]
                self.defects['cds_seq_not_matching'] += 1

        except TranslationError as e:
            out_feat['warnings'] = out_feat.get('warnings', []) + [
                "Unable to verify protein sequence:" + str(e)]

        if not prot_seq:
            try:
                prot_seq = Seq.translate(
                        feat_seq, self.code_table, cds=True).strip("*")
                out_feat['warnings'] = out_feat.get('warnings', []) + [
                        warnings["no_translation_supplied"]]
            except TranslationError as e:
                out_feat['warnings'] = out_feat.get('warnings', []) + [
                    warnings["no_translation_supplied"] + str(e)]

        out_feat.update({
            "protein_translation": prot_seq,
            "protein_md5": hashlib.md5(prot_seq.encode('utf8')).hexdigest(),
            "protein_translation_length": len(prot_seq),
        })

        if out_feat.get('parent_gene'):
            propagate_cds_props_to_gene(out_feat, self.genes[gene_id])

        self.cdss[out_feat['id']] = out_feat
