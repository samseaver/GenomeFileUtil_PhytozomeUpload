
import copy
import datetime
import time
import os
import re
import sys
import shutil
import uuid
import hashlib
import json
from collections import Counter, defaultdict, OrderedDict

import Bio.SeqIO
import Bio.SeqUtils
from Bio import Seq
from Bio.Data.CodonTable import TranslationError
from Bio.SeqFeature import ExactPosition

from AssemblyUtil.AssemblyUtilClient import AssemblyUtil
from DataFileUtil.DataFileUtilClient import DataFileUtil
from GenomeInterface import GenomeInterface
from GenomeUtils import is_parent, propagate_cds_props_to_gene
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
        self.generate_parents = False
        self.ontologies_present = defaultdict(dict)
        self.skiped_features = Counter()
        self.feature_counts = Counter()
        self.orphan_types = Counter()
        self.contig_seq = {}
        self.circ_contigs = set()
        self.features_spaning_zero = set()
        self.genome_warnings = []
        self.genome_suspect = False
        self.cds_seq_not_matching = 0
        self.bad_parent_loc = 0
        self.excluded_features = ('source', 'exon')
        self.go_mapping = json.load(
            open('/kb/module/data/go_ontology_mapping.json'))
        self.po_mapping = json.load(
            open('/kb/module/data/go_ontology_mapping.json'))
        self.ontology_events = [
            {
                "method": "GenomeFileUtils Genbank uploader from annotations",
                "method_version": self.version,
                "timestamp": self.time_string,
                "id": "GO",
                "ontology_ref": "KBaseOntology/gene_ontology"
            },
            {
                "method": "GenomeFileUtils Genbank uploader from annotations",
                "method_version": self.version,
                "timestamp": self.time_string,
                "id": "PO",
                "ontology_ref": "KBaseOntology/plant_ontology"
            }
        ]
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
            'exclude_ontologies': 0,
            'type': 'User upload',
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

        # 4) Do the upload
        files = self._find_input_files(input_directory)
        consolidated_file = self._join_files_skip_empty_lines(files)
        genome = self.parse_genbank(consolidated_file, params)
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
        assembly_ref = self._save_assembly(file_path, params)
        assembly_data = self.dfu.get_objects(
            {'object_refs': [assembly_ref],
             'ignore_errors': 0})['data'][0]['data']
        genome = {
            "id": params['genome_name'],
            "type": params['type'],
            "original_source_file_name": os.path.basename(file_path),
            "assembly_ref": assembly_ref,
            "gc_content": assembly_data['gc_content'],
            "dna_size": assembly_data['dna_size'],
            "md5": assembly_data['md5'],
            "genbank_handle_ref": shock_res['handle']['hid'],
            "publications": set(),
            "contig_ids": [],
            "contig_lengths": [],
            "features": [],
            "non_coding_features": [],
            "cdss": [],
            'mrnas': [],
        }
        genome['source'], genome['genome_tiers'] = \
            self.gi.determine_tier(params['source'])

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

            # extends the feature, cdss, mrna & non_coding_genes arrays
            for k, v in self._parse_features(record, params['source']).items():
                genome[k].extend(v)

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

        if len(genome['cdss']) and self.cds_seq_not_matching / float(len(
                genome['cdss'])) > 0.01:
            self.genome_warnings.append("SUSPECT This Genome has a high "
                "proportion ({} out of {}) CDS features that do not "
                "translate the supplied translation".format(
                self.cds_seq_not_matching, len(genome['cdss'])))
            self.genome_suspect = 1
        if self.bad_parent_loc:
            self.genome_warnings.append("There were {} parent/child "
                "relationships that were not able to be determined. Some of "
                "these may have splice variants that may be valid "
                "relationships.".format(self.bad_parent_loc))

        if self.genome_warnings:
            genome['warnings'] = self.genome_warnings
        if self.genome_suspect:
            genome['suspect'] = 1
        print("Feature Counts: ", genome['feature_counts'])
        return genome

    def _save_assembly(self, genbank_file, params):
        """Convert genbank file to fasta and sve as assembly"""
        contigs = Bio.SeqIO.parse(genbank_file, "genbank")
        print("Saving sequence as Assembly object")
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
            self.contig_seq[in_contig.id] = in_contig.seq

        Bio.SeqIO.write(out_contigs, fasta_file, "fasta")
        assembly_ref = self.aUtil.save_assembly_from_fasta(
            {'file': {'path': fasta_file},
             'workspace_name': params['workspace_name'],
             'assembly_name': assembly_id,
             'contig_info': extra_info})
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
                tags = ['locus_tag', 'gene', 'kbase_id']
            while tags and not _id:
                _id = feat.qualifiers.get(tags.pop(0), [""])[0]
            return _id

        def _location(feat):
            """Convert to KBase style location objects"""
            strand_trans = ("", "+", "-")
            loc = []
            for part in feat.location.parts:
                if part.strand >= 0:
                    begin = int(part.start) + 1
                else:
                    begin = int(part.end)
                loc.append((
                        record.id,
                        begin,
                        strand_trans[part.strand],
                        len(part)))
            return loc

        def _warn(message):
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

            msg = "The feature coordinates order are suspect and the " \
                   "feature is not flagged as being trans-spliced"
            _warn(msg)
            self.genome_warnings.append("SUSPECT {}:".format(out_feat['id']) +
                                        msg)
            self.genome_suspect = True

        genes, cdss, mrnas, noncoding = OrderedDict(), OrderedDict(), OrderedDict(), []
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
                "md5": hashlib.md5(str(feat_seq)).hexdigest(),
            }
            if not _id:
                out_feat['id'] = in_feature.type

            # validate input feature
            # note that end is the larger number regardless of strand
            if int(in_feature.location.end) > len(record):
                self.genome_warnings.append('SUSPECT {}: Feature has '
                    'invalid coordinates off of the end of the contig and was '
                    'not included'.format(out_feat['id']))
                self.genome_suspect = True
                continue

            for piece in in_feature.location.parts:
                if not isinstance(piece.start, ExactPosition) \
                        or not isinstance(piece.end, ExactPosition):
                    _warn("The coordinates supplied for this feature are "
                          "non-exact. DNA or protein translations are "
                          "approximate.")

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
                out_feat['inference_data'] = self._inferences(in_feature)

            _check_suspect_location(genes.get(_id))

            # add type specific features
            if in_feature.type == 'CDS':
                self.process_cds(
                    _id, feat_seq, genes, in_feature, mrnas, out_feat, cdss)

            elif in_feature.type == 'gene':
                self.process_gene(_id, genes, out_feat)

            elif in_feature.type == 'mRNA':
                self.process_mrna(_id, genes, out_feat, mrnas)

            else:
                noncoding.append(self.process_noncodeing(_id, genes,
                                                         in_feature, out_feat))

        # sort genes into their final arrays
        coding = []
        for g in genes.values():
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
                coding.append(g)
                self.feature_counts["protein_encoding_gene"] += 1
            # genes that have no valid CDS after location tests are excluded
            elif warnings['child_cds_failed'] in g.get('warnings', []):
                self.genome_warnings.append(warnings['gene_excluded'].format(
                    g['id']))
                self.genome_suspect = True
                continue
            else:
                del g['mrnas'], g['cdss']
                noncoding.append(g)
                self.feature_counts["non-protein_encoding_gene"] += 1

        return {'features': coding, 'non_coding_features': noncoding,
                'cdss': cdss.values(), 'mrnas': mrnas.values()}

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

    def _get_ontology_db_xrefs(self, feature):
        """Splits the ontology info from the other db_xrefs"""
        ontology = defaultdict(dict)
        db_xref = []
        for key in ("GO_process", "GO_function", "GO_component"):
            for term in feature.qualifiers.get(key, []):
                sp = term.split(" - ")
                ontology['GO'][sp[0]] = [1]
                self.ontologies_present['GO'][sp[0]] = sp[1]
        for ref in feature.qualifiers.get('db_xref', []):
            if ref.startswith('GO:'):
                ontology['GO'][ref] = [0]
                self.ontologies_present['GO'][ref] = self.go_mapping.get(ref, '')
            elif ref.startswith('PO:'):
                ontology['PO'][ref] = [1]
                self.ontologies_present['PO'][ref] = self.po_mapping.get(ref, '')
            else:
                db_xref.append(tuple(ref.split(":")))
        # TODO: Support other ontologies
        return dict(ontology), db_xref

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
                result['functions'].extend(val_list[0].split('; '))
            if key == 'product':
                result['functions'].append("product:" + val_list[0])

        return result

    @staticmethod
    def _inferences(feat):
        """Whoever designed the genbank delimitation is an idiot: starts and
        ends with a optional values and uses a delimiter ":" that is
        used to divide it's DBs in the evidence. Anyway, this sorts that"""
        result = []
        for inf in feat.qualifiers['inference']:
            try:
                sp_inf = inf.split(":")
                if sp_inf[0] in ('COORDINATES', 'DESCRIPTION', 'EXISTENCE'):
                    inference = {'category': sp_inf.pop(0)}
                else:
                    inference = {'category': ''}
                inference['type'] = sp_inf[0]
                inference['evidence'] = ":".join(sp_inf[1:])
                result.append(inference)
            except IndexError('Unparseable inference string: ' + inf):
                continue
        return result

    @staticmethod
    def process_gene(_id, genes, out_feat):
        out_feat.update({
            "id": _id,
            "type": 'gene',
            "mrnas": [],
            'cdss': [],
        })
        genes[_id] = out_feat

    def process_noncodeing(self, _id, genes, in_feature, out_feat):
        out_feat["type"] = in_feature.type
        # add increment number of each type
        if _id in genes:
            if not is_parent(genes[_id], out_feat):
                out_feat['warnings'] = out_feat.get('warnings', []) + [
                    "Feature order suggests that {} is the parent gene, but it"
                    " fails location validation".format(_id)]
                self.bad_parent_loc += 1
            else:
                if 'children' not in genes[_id]:
                    genes[_id]['children'] = []
                out_feat['id'] += "_" + str(len(genes[_id]['children']) + 1)
                genes[_id]['children'].append(out_feat['id'])
                out_feat['parent_gene'] = _id
        else:
            self.orphan_types[in_feature.type] += 1
            out_feat['id'] += "_" + str(self.orphan_types[in_feature.type])
        return out_feat

    def process_cds(self, _id, feat_seq, genes, in_feature, mrnas, out_feat, cdss):
        prot_seq = in_feature.qualifiers.get("translation", [""])[0]
        out_feat.update({
            "protein_translation": prot_seq,
            "protein_md5": hashlib.md5(prot_seq).hexdigest(),
            "protein_translation_length": len(prot_seq),
        })
        if not prot_seq:
            try:
                out_feat['protein_translation'] = Seq.translate(
                        feat_seq, self.code_table, cds=True).strip("*")
                out_feat['warnings'] = out_feat.get('warnings', []) + [
                    "This CDS did not have a supplied translation. The "
                    "translation is derived directly from DNA sequence."]
            except TranslationError as e:
                out_feat['warnings'] = out_feat.get('warnings', []) + [
                    "Protein translation not supplied. Unable to generate "
                    "protein sequence:" + str(e)]

        # allow a little slack to account for frameshift and stop codon
        if prot_seq and abs(len(prot_seq) * 3 - len(feat_seq)) > 4:
            out_feat['warnings'] = out_feat.get('warnings', []) + [
                "This CDS has a length of {} which is not consistent with the "
                "length of the translation included ({} amino acids)".format(
                    len(feat_seq), len(prot_seq))]
            self.genome_warnings.append('SUSPECT {}: This CDS has a length of '
                    '{} which is not consistent with the length of the '
                    'translation included ({} amino acids)'.format(
                out_feat['id'], len(feat_seq), len(prot_seq)))
            self.genome_suspect = True

        try:
            if prot_seq and prot_seq != Seq.translate(
                    feat_seq, self.code_table, cds=True).strip("*"):
                out_feat['warnings'] = out_feat.get('warnings', []) + [
                    "The annotated protein translation is not "
                    "consistent with the recorded DNA sequence"]
                self.cds_seq_not_matching += 1

        except TranslationError as e:
            out_feat['warnings'] = out_feat.get('warnings', []) + [
                "Unable to verify protein sequence:" + str(e)]

        if _id not in genes and self.generate_parents:
            self.process_gene(_id, genes, copy.copy(out_feat))

        if _id in genes:
            out_feat['id'] += "_" + str(len(genes[_id]['cdss']) + 1)
            if not is_parent(genes[_id], out_feat):
                genes[_id]['warnings'] = genes[_id].get('warnings', []) + [
                    warnings['child_cds_failed']]
                self.genome_warnings.append(
                    warnings['cds_excluded'].format(_id))
                self.genome_suspect = True
                self.bad_parent_loc += 1
                return
            else:
                genes[_id]['cdss'].append(out_feat['id'])
                propagate_cds_props_to_gene(out_feat, genes[_id])
                out_feat['parent_gene'] = _id
        else:
            out_feat['parent_gene'] = ""

        mrna_id = out_feat["id"].replace('CDS', 'mRNA')
        if mrna_id in mrnas:
            if not is_parent(mrnas[mrna_id], out_feat):
                out_feat['warnings'] = out_feat.get('warnings', []) + [
                    warnings['cds_mrna_cds'].format(mrna_id)]
                mrnas[mrna_id]['warnings'] = mrnas[mrna_id].get(
                    'warnings', []) + [warnings['cds_mrna_mrna']]
                self.bad_parent_loc += 1
            else:
                out_feat['parent_mrna'] = mrna_id
                mrnas[mrna_id]['cds'] = out_feat['id']
        cdss[out_feat['id']] = out_feat

    def process_mrna(self, _id, genes, out_feat, mrnas):
        if _id not in genes and self.generate_parents:
                self.process_gene(_id, genes, copy.copy(out_feat))

        if _id in genes:
            out_feat['id'] += "_" + str(len(genes[_id]['mrnas']) + 1)
            if not is_parent(genes[_id], out_feat):
                genes[_id]['warnings'] = genes[_id].get('warnings', []) + [
                    warnings['child_mrna_failed']]
                self.genome_warnings.append(
                    warnings['mrna_excluded'].format(_id))
                self.genome_suspect = True
                self.bad_parent_loc += 1
                return
            else:
                genes[_id]['mrnas'].append(out_feat['id'])
                out_feat['parent_gene'] = _id
        else:
            out_feat['warnings'] = out_feat.get('warnings', []) + [
                'Unable to find parent mrna for ' + str(out_feat)]
            out_feat['parent_gene'] = ""
        mrnas[out_feat['id']] = out_feat


warnings = {
    "cds_excluded": "SUSPECT CDS {} was excluded because the associated CDS "
                    "failed coordinates validation",

    "cds_mrna_cds": "Feature order suggests that {} is the parent mRNA, but it"
                    " fails location validation",
    "cds_mrna_mrna": "Potential child CDS relationship failed due to location "
                    "validation.",

    "child_cds_failed": "The child CDS failed location validation. That CDS "
                        "has been excluded.",
    "child_mrna_failed": "The child mRNA failed location validation. That mRNA"
                         " has been excluded.",
    "gene_excluded": "SUSPECT gene {} was excluded because the associated CDS "
                     "failed coordinates validation",
    "mrna_excluded": "SUSPECT mRNA {} was excluded because the associated mRNA "
                    "failed coordinates validation"
}