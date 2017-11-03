"""
GenomeAnnotation to GenBank file conversion.
"""

# Stdlib
from collections import defaultdict
import time

from Bio import SeqIO, SeqFeature, Alphabet

# Local
from DataFileUtil.DataFileUtilClient import DataFileUtil
from AssemblyUtil.AssemblyUtilClient import AssemblyUtil
from GenomeAnnotationAPI.GenomeAnnotationAPIClient import GenomeAnnotationAPI

STD_PREFIX = " " * 21


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

        # 2) get genome genbank handle reference
        getGenomeOptions = {
            'genomes': [{
                'ref': params['genome_ref']
            }],
            'included_fields': ['genbank_handle_ref'],
            'ignore_errors': 0  # if we can't find the genome, throw an error
        }
        if 'ref_path_to_genome' in params:
            getGenomeOptions['genomes'][0]['ref_path_to_genome'] = params['ref_path_to_genome']

        api = GenomeAnnotationAPI(self.cfg.callbackURL)
        genome_data = api.get_genome_v1(getGenomeOptions)['genomes'][0]
        info = genome_data['info']
        data = genome_data['data']

        # 3) make sure the type is valid
        if info[2].split('-')[0] != 'KBaseGenomes.Genome':
            raise ValueError('Object is not a Genome, it is a:' + str(info[2]))

        # 4) build the genbank file and return it
        print('not cached, building file...')
        result = self.build_genbank_file(params['genome_ref'],
                                         "KBase_derived_" + info[1] + ".gbff")
        if result is None:
            raise ValueError('Unable to generate file.  Something went wrong')
        result['from_cache'] = 0
        return result

    def export_original_genbank(self, ctx, params):
        # 1) validate parameters and extract defaults
        self.validate_params(params)

        # 2) get genome genbank handle reference
        getGenomeOptions = {
            'genomes':[{
                'ref': params['genome_ref']
            }],
            'included_fields':['genbank_handle_ref'],
            'ignore_errors':0 # if we can't find the genome, throw an error
        }
        if 'ref_path_to_genome' in params:
            getGenomeOptions['genomes'][0]['ref_path_to_genome'] = params['ref_path_to_genome']

        api = GenomeAnnotationAPI(self.cfg.callbackURL)
        genome_data = api.get_genome_v1(getGenomeOptions)['genomes'][0]
        info = genome_data['info']
        data = genome_data['data']

        # 3) make sure the type is valid
        if info[2].split('-')[0] != 'KBaseGenomes.Genome':
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

    def build_genbank_file(self, genome_ref, output_filename):
        genome = self.dfu.get_objects({
            'object_refs': [genome_ref]
        })['data'][0]
        g = GenomeFile(self.cfg, genome['data'])
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
        feat_arrays = (('features', 'gene'), ('cdss', 'CDS'),
                       ('mrnas', 'mRNA'), ('non_coding_features', ''))
        # sort every feature in the feat_arrays into a dict by contig
        for key, type in feat_arrays:
            if key not in genome_object:
                continue
            for feat in genome_object[key]:
                if type:
                    feat['type'] = type
                self.features_by_contig[feat['location'][0][0]].append(feat)

        assembly_file_path = self._get_assembly(genome_object['assembly_ref'])
        for contig in SeqIO.parse(open(assembly_file_path), 'fasta',
                                  Alphabet.generic_dna):
            self._parse_contig(contig)

    def _get_assembly(self, assembly_ref):
        print('Assembly reference = ' + assembly_ref)
        print('Downloading assembly')
        au = AssemblyUtil(self.cfg.callbackURL)
        assembly_file_path = au.get_assembly_as_fasta(
            {'ref': assembly_ref}
        )['path']
        return assembly_file_path

    def _parse_contig(self, raw_contig):
        go = self.genome_object  # I'm lazy
        raw_contig.dbxrefs = self.genome_object.get('aliases', [])
        raw_contig.annotations = {
            "comment": go.get('notes', ""),
            "source": "KBase_" + go.get('source', ""),
            "taxonomy": go.get('taxonomy', "").split("; "),
            "organism": go.get('scientific_name', ""),
            "date": time.strftime("%d-%b-%Y",
                                  time.localtime(time.time())).upper()
        }
        if not self.seq_records:  # Only on the first contig
            raw_contig.annotations['references'] = self._format_publications()
            print("Added {} references".format(
                len(raw_contig.annotations['references'])))

        if raw_contig.id in self.features_by_contig:
            for feat in self.features_by_contig[raw_contig.id]:
                raw_contig.features.append(self._format_feature(feat))
            raw_contig.features.sort()

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

    @staticmethod
    def _format_feature(in_feature):
        def _trans_loc(loc):
            if loc[2] == "-":
                return SeqFeature.FeatureLocation(loc[1]-1, loc[1]-loc[3]-1, -1)
            else:
                return SeqFeature.FeatureLocation(loc[1]-1, loc[1]+loc[3]-1, 1)

        location = _trans_loc(in_feature['location'].pop())
        while in_feature['location']:
            location += _trans_loc(in_feature['location'].pop())

        out_feature = SeqFeature.SeqFeature(location, in_feature['type'])
        if 'function' in in_feature and in_feature['function']:
            out_feature.qualifiers['product'] = in_feature['function']
        if 'aliases' in in_feature and in_feature['aliases']:
            out_feature.qualifiers['db_xref'] = in_feature['aliases']
        if 'note' in in_feature and in_feature['note']:
            out_feature.qualifiers['note'] = in_feature['note']
        if 'protein_translation' in in_feature and in_feature['protein_translation']:
            out_feature.qualifiers['translation'] = in_feature['protein_translation']
        # TODO: Ontologies
        return out_feature

    def write_genbank_file(self, file_path):
        if not self.seq_records:
            raise ValueError("No sequence data to write!")
        self._sort_features()
        SeqIO.write(self.seq_records, open(file_path, 'w'), 'genbank')

    def _sort_features(self):
        def feature_sort(feat):
            order = ('gene', 'mRNA', 'CDS')
            if feat.type not in order:
                priority = len(order)
            else:
                priority = order.index(feat.type)
            return feat.location.start, priority
        for contig in self.seq_records:
            contig.features.sort(key=feature_sort)
