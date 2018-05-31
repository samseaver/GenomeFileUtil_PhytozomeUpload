"""
GenomeAnnotation to GenBank file conversion.
"""

# Stdlib
from collections import defaultdict
import time
import logging

from Bio import SeqIO, SeqFeature, Alphabet

# Local
from DataFileUtil.DataFileUtilClient import DataFileUtil
from AssemblyUtil.AssemblyUtilClient import AssemblyUtil
STD_PREFIX = " " * 21
CONTIG_ID_FIELD_LENGTH = 16


class GenomeToGenbank(object):

    def __init__(self, sdk_config):
        self.cfg = sdk_config
        self.dfu = DataFileUtil(self.cfg.callbackURL)

    def validate_params(self, params):
        if 'genome_ref' not in params:
            raise ValueError('required "genome_ref" field was not defined')

    def export(self, ctx, params):
        # 1) validate parameters and extract defaults
        self.validate_params(params)

        # 2) get genome info
        genome_data = self.dfu.get_objects({
            'object_refs': [params['genome_ref']]
        })['data'][0]
        info = genome_data['info']
        data = genome_data['data']

        # 3) make sure the type is valid
        if info[2].split(".")[1].split('-')[0] != 'Genome':
            raise ValueError('Object is not a Genome, it is a:' + str(info[2]))

        # 4) build the genbank file and return it
        print('not cached, building file...')
        result = self.build_genbank_file(data,
                                         "KBase_derived_" + info[1] + ".gbff")
        if result is None:
            raise ValueError('Unable to generate file.  Something went wrong')
        result['from_cache'] = 0
        return result

    def export_original_genbank(self, ctx, params):
        # 1) validate parameters and extract defaults
        self.validate_params(params)

        # 2) get genome genbank handle reference
        genome_data = self.dfu.get_objects({
            'object_refs': [params['genome_ref']]
        })['data'][0]
        info = genome_data['info']
        data = genome_data['data']

        # 3) make sure the type is valid
        if info[2].split(".")[1].split('-')[0] != 'Genome':
            raise ValueError('Object is not a Genome, it is a:' + str(info[2]))

        # 4) if the genbank handle is there, get it and return
        print('checking if genbank file is cached...')
        result = self.get_genbank_handle(data)
        return result

    def get_genbank_handle(self, data):
        if 'genbank_handle_ref' not in data:
            return None
        if data['genbank_handle_ref'] is None:
            return None

        print('pulling cached genbank file from Shock: ' +
              str(data['genbank_handle_ref']))
        file = self.dfu.shock_to_file({
                            'handle_id': data['genbank_handle_ref'],
                            'file_path': self.cfg.sharedFolder,
                            'unpack': 'unpack'
                        })
        return {
            'genbank_file': {
                'file_path': file['file_path']
            }
        }

    def build_genbank_file(self, genome_data, output_filename):
        g = GenomeFile(self.cfg, genome_data)
        file_path = self.cfg.sharedFolder + "/" + output_filename
        g.write_genbank_file(file_path)

        return {
            'genbank_file': {
                'file_path': file_path
            }
        }


class GenomeFile:
    def __init__(self, cfg, genome_object):
        self.cfg = cfg
        self.genome_object = genome_object
        self.seq_records = []
        self.features_by_contig = defaultdict(list)
        # make special dict for all mrna & cds, they will be added when their
        # parent gene is added
        self.child_dict = {}
        for x in genome_object.get('mrnas', []):
            x['type'] = 'mRNA'
            self.child_dict[x['id']] = x

        for x in genome_object.get('cdss', []):
            x['type'] = 'CDS'
            self.child_dict[x['id']] = x

        # sort other features into a dict by contig
        for feat in genome_object['features']:
            if 'type' not in feat:
                feat['type'] = 'gene'
            self.features_by_contig[feat['location'][0][0]].append(feat)
        for feat in genome_object.get('non_coding_features', []):
            self.features_by_contig[feat['location'][0][0]].append(feat)

        assembly_file_path, circ_contigs = self._get_assembly(genome_object)
        for contig in SeqIO.parse(open(assembly_file_path), 'fasta', Alphabet.generic_dna):
            if contig.id in circ_contigs:
                contig.annotations['topology'] = "circular"
            self._parse_contig(contig)

    def _get_assembly(self, genome):
        if 'assembly_ref' in genome:
            assembly_ref = genome['assembly_ref']
        else:
            assembly_ref = genome['contigset_ref']
        print('Assembly reference = ' + assembly_ref)
        print('Downloading assembly')
        dfu = DataFileUtil(self.cfg.callbackURL)
        assembly_data = dfu.get_objects({
            'object_refs': [assembly_ref]
        })['data'][0]['data']
        if isinstance(assembly_data['contigs'], dict):  # is an assembly
            circular_contigs = set([x['contig_id'] for x in assembly_data['contigs'].values()
                                    if x.get('is_circ')])
        else:  # is a contig set
            circular_contigs = set([x['id'] for x in assembly_data['contigs']
                                    if x.get('replicon_geometry') == 'circular'])
        au = AssemblyUtil(self.cfg.callbackURL)
        assembly_file_path = au.get_assembly_as_fasta(
            {'ref': assembly_ref}
        )['path']
        return assembly_file_path, circular_contigs

    def _parse_contig(self, raw_contig):
        def feature_sort(feat):
            order = ('gene', 'mRNA', 'CDS')
            if feat['type'] not in order:
                priority = len(order)
            else:
                priority = order.index(feat['type'])
            start = min(x[1] for x in feat['location'])
            return start, priority

        go = self.genome_object  # I'm lazy
        raw_contig.dbxrefs = self.genome_object.get('aliases', [])
        raw_contig.annotations.update({
            "comment": go.get('notes', ""),
            "source": "KBase_" + go.get('source', ""),
            "taxonomy": go.get('taxonomy', "").split("; "),
            "organism": go.get('scientific_name', ""),
            "date": time.strftime("%d-%b-%Y",
                                  time.localtime(time.time())).upper()
        })
        if not self.seq_records:  # Only on the first contig
            raw_contig.annotations['references'] = self._format_publications()
            print("Added {} references".format(
                len(raw_contig.annotations['references'])))
            if 'notes' in go:
                raw_contig.annotations['comment'] = go['notes']

        if raw_contig.id in self.features_by_contig:
            # sort all features except for cdss and mrnas
            self.features_by_contig[raw_contig.id].sort(key=feature_sort)
            for feat in self.features_by_contig[raw_contig.id]:
                raw_contig.features.append(self._format_feature(feat))
                # process child mrnas & cdss if present
                raw_contig.features.extend([self._format_feature(
                    self.child_dict[_id]) for _id in feat.get('mrnas', [])])
                raw_contig.features.extend([self._format_feature(
                    self.child_dict[_id]) for _id in feat.get('cdss', [])])
        if len(raw_contig.name) > CONTIG_ID_FIELD_LENGTH:
            logging.info("Trimming contig name {} to the last {} characters"
                         .format(raw_contig.name, CONTIG_ID_FIELD_LENGTH))
            raw_contig.name = raw_contig.name[-CONTIG_ID_FIELD_LENGTH:]

        self.seq_records.append(raw_contig)

    def _format_publications(self):
        references = []
        for pub in self.genome_object.get('publications', []):
            if len(pub) != 7:
                print('Skipping unparseable publication {}'.format(pub))
            ref = SeqFeature.Reference()
            if pub[0]:
                ref.pubmed_id = str(pub[0])
            ref.title = pub[2]
            ref.authors = pub[5]
            ref.journal = pub[6]
            references.append(ref)
        return references

    def _format_feature(self, in_feature):
        def _trans_loc(loc):
            if len(loc[0]) > CONTIG_ID_FIELD_LENGTH:
                loc[0] = loc[0][-CONTIG_ID_FIELD_LENGTH:]
            if loc[2] == "-":
                return SeqFeature.FeatureLocation(loc[1]-loc[3], loc[1], -1, loc[0])
            else:
                return SeqFeature.FeatureLocation(loc[1]-1, loc[1]+loc[3]-1, 1, loc[0])

        # we have to do it this way to correctly make a "CompoundLocation"
        location = _trans_loc(in_feature['location'].pop(0))
        while in_feature['location']:
            location += _trans_loc(in_feature['location'].pop(0))
        out_feature = SeqFeature.SeqFeature(location, in_feature['type'])

        if in_feature.get('functional_descriptions'):
            out_feature.qualifiers['function'] = "; ".join(
                    in_feature['functional_descriptions'])
        if in_feature.get('functions'):
            out_feature.qualifiers['product'] = "; ".join(in_feature['functions'])
        if 'function' in in_feature:
            out_feature.qualifiers['product'] = in_feature['function']

        if in_feature.get('note', False):
            out_feature.qualifiers['note'] = in_feature['note']
        if in_feature.get('protein_translation', False):
            out_feature.qualifiers['translation'] = in_feature['protein_translation']
        if in_feature.get('db_xrefs', False):
            out_feature.qualifiers['db_xref'] = ["{}:{}".format(*x) for x in
                                                 in_feature['db_xrefs']]
        if in_feature.get('ontology_terms', False):
            if 'db_xref' not in out_feature.qualifiers:
                out_feature.qualifiers['db_xref'] = []
            for ont, terms in in_feature['ontology_terms'].items():
                out_feature.qualifiers['db_xref'].extend([t for t in terms])

        for alias in in_feature.get('aliases', []):
            if len(alias) == 2:
                if not alias[0] in out_feature.qualifiers:
                    out_feature.qualifiers[alias[0]] = []
                out_feature.qualifiers[alias[0]].append(alias[1])
            else:  # back compatibility
                if 'db_xref' not in out_feature.qualifiers:
                    out_feature.qualifiers['db_xref'] = []
                out_feature.qualifiers['db_xref'].append(alias)

        for flag in in_feature.get('flags', []):
            out_feature.qualifiers[flag] = None

        if 'inference_data' in in_feature:
            out_feature.qualifiers['inference'] = [
                ":".join([x[y] for y in ('category', 'type', 'evidence') if x[y]])
                for x in in_feature['inference_data']]

        if in_feature.get('warnings', False):
            out_feature.qualifiers['note'] = out_feature.qualifiers.get(
                'note', "") + "Warnings: " + ",".join(in_feature['warnings'])

        return out_feature

    def write_genbank_file(self, file_path):
        if not self.seq_records:
            raise ValueError("No sequence data to write!")
        SeqIO.write(self.seq_records, open(file_path, 'w'), 'genbank')
