
import os
import time

from biokbase.workspace.client import Workspace
from DataFileUtil.DataFileUtilClient import DataFileUtil

from GenomeAnnotationAPI.GenomeAnnotationAPIClient import GenomeAnnotationAPI


class GenomeToGFF:
    '''
    typedef structure {
        string genome_ref;
        list <string> ref_path_to_genome;
        int is_gtf;
    } GenomeToGFFParams;

    /* from_cache is 1 if the file already exists and was just returned, 0 if
    the file was generated during this call. */
    typedef structure {
        File file_path;
        boolean from_cache;
    } GenomeToGFFResult;

    funcdef genome_to_gff(GenomeToGFFParams params)
                returns (GenomeToGFFResult result) authentication required;
    '''

    def __init__(self, sdk_config):
        self.cfg = sdk_config;


    def export(self, ctx, params):
        # 1) validate parameters and extract defaults
        self.validate_params(params)

        # 2) get genome gff handle reference
        getGenomeOptions = {
            'genomes':[{
                'ref':params['genome_ref']
            }],
            'included_fields':['gff_handle_ref'],
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

        is_gtf = params.get('is_gtf', 0)

        target_dir = params.get('target_dir')
        if not target_dir:
            target_dir = os.path.join(self.cfg.sharedFolder, "gff_" + str(int(time.time() * 1000)))
        if not os.path.exists(target_dir):
            os.makedirs(target_dir)

        # 4) if the GFF handle is there, get it and return
        if is_gtf != 1:
            print('checking if GFF file is cached...')
            result = self.get_gff_handle(data, target_dir)
            if result is not None:
                result['from_cache'] = 1
                return result
            print('not cached, building file...')

        # 5) otherwise, build the GFF file and return it
        result = self.build_gff_file(getGenomeOptions, target_dir, info[1], is_gtf == 1)
        if result is None:
            raise ValueError('Unable to generate file.  Something went wrong')
        result['from_cache'] = 0
        return result


    def get_gff_handle(self, data, output_dir):

        if 'gff_handle_ref' not in data:
            return None
        if data['gff_handle_ref'] is None:
            return None

        print('pulling cached GFF file from Shock: '+str(data['gff_handle_ref']))
        dfu = DataFileUtil(self.cfg.callbackURL)
        file_ret = dfu.shock_to_file({'handle_id':data['gff_handle_ref'],
                                      'file_path':output_dir,
                                      'unpack': 'unpack'})
        return {
            'file_path': file_ret['file_path']
        }


    ###  see logic from: https://github.com/kbase/KBaseRNASeq/blob/e2d69e4137903c68b5a1fedeab57f7900aad7253/lib/biokbase/RNASeq/KBaseRNASeqImpl.py#L443-L562
    def build_gff_file(self, getGenomeOptions, output_dir, output_filename, is_gtf):

        # first get subdata needed; forget about the metadata
        #getGenomeOptions['included_fields'] = []
        #getGenomeOptions['included_feature_fields'] = ['id', 'type', 'location']
        getGenomeOptions['no_metadata'] = 1
        if 'included_fields' in getGenomeOptions:
            del getGenomeOptions['included_fields']
        if 'included_feature_fields' in getGenomeOptions:
            del getGenomeOptions['included_feature_fields']

        api = GenomeAnnotationAPI(self.cfg.callbackURL)
        genome_data = api.get_genome_v1(getGenomeOptions)['genomes'][0]['data']

        # create the file
        try:
            file_ext = ".gtf" if is_gtf else ".gff"
            out_file_path = os.path.join(output_dir, output_filename + file_ext)
            print('Creating file: '+ str(out_file_path))
            output = open(out_file_path,'w')
            features = []
            if 'features' in genome_data:
                for f in genome_data['features']:
                    features.append({'id': f['id'], 'type': f['type'], 'location': f['location']})
            if 'cdss' in genome_data:
                for f in genome_data['cdss']:
                    features.append({'id': f['id'], 'type': 'CDS', 'location': f['location'], 
                                    'parent_gene': f['parent_gene'], 
                                    'parent_mrna': f['parent_mrna']})
            if 'mrnas' in genome_data:
                for f in genome_data['mrnas']:
                    features.append({'id': f['id'], 'type': 'mRNA', 'location': f['location'],
                                     'parent_gene': f['parent_gene']})
            mrna_map = {}   ## mrna_id -> <mRNA>
            gene_map = {}   ## gene_id -> <gene>
                            ## gene is {'id': <>, 'location': [[contig,start,strand,len], ...], 
                            ##          'mrna_cds_pairs': [[<mRNA>, <CDS>], ...]}
            #gene_id_generation = 1
            #mrna_id_generation = 1
            for f in features:
                if f['type'] == 'mRNA':
                    mrna_map[f['id']] = f
                elif f['type'] != 'CDS':
                    gene_map[f['id']] = f
            ## Now let's go over CDSs
            for f in features:
                if f['type'] == 'CDS':
                    gene_id = f.get('parent_gene')
                    gene = None
                    if gene_id:
                        gene = gene_map.get(gene_id)
                    rename_cds = False
                    if gene is None:
                        if gene_id is None:
                            gene_id = f['id']  #'gene_' + str(gene_id_generation)
                            #gene_id_generation += 1
                            rename_cds = True
                        gene = {'id': gene_id, 'location': self.get_common_location(f['location'])}
                        gene_map[gene_id] = gene
                    mrna_id = f.get('parent_mrna')
                    mrna = None
                    if mrna_id:
                        mrna = mrna_map.get(mrna_id)
                    if mrna is None:
                        if mrna_id is None:
                            mrna_id = f['id'] + '_mRNA'  # 'mRNA_' + str(mrna_id_generation)
                            #mrna_id_generation += 1
                        mrna = {'id': mrna_id, 'location': f['location']}
                        mrna_map[mrna_id] = mrna
                    if rename_cds:
                        f['id'] = f['id'] + '_CDS'
                    mrna_cds_pairs = gene.get('mrna_cds_pairs')
                    if mrna_cds_pairs is None:
                        mrna_cds_pairs = []
                        gene['mrna_cds_pairs'] = mrna_cds_pairs
                    mrna_cds_pairs.append([mrna, f])
            ## Let's sort genes by contigs
            contigs = []  ## contig is {'genes': []}
            contig_map = {}
            for gene_id in gene_map:
                gene = gene_map[gene_id]
                gene['start'] = self.get_start(gene['location'][0])
                contig_id = gene['location'][0][0]
                contig = contig_map.get(contig_id)
                if contig is None:
                    contig = {'id': contig_id, 'genes': []}
                    contig_map[contig_id] = contig
                    contigs.append(contig)
                contig['genes'].append(gene)

            for contig in contigs:
                contig['genes'].sort(key=lambda gene: gene['start'])

            # write the file
            exon_id_generation = 1
            for contig in contigs:
                contig_id = contig['id']
                for gene in contig['genes']:
                    gene_id = gene['id']
                    strand = gene['location'][0][2]
                    if not is_gtf:
                        self.write_gff_line(output, contig_id, 'gene', gene['start'],
                                            self.get_end(gene['location'][0]), strand, '.', 
                                            gene_id, None)
                    if 'mrna_cds_pairs' not in gene:
                        continue
                    for [mrna, cds] in gene['mrna_cds_pairs']:
                        mrna_id = mrna['id']
                        mrna_loc = self.get_common_location(mrna['location'])[0]
                        if not is_gtf:
                            self.write_gff_line(output, contig_id, 'mRNA', 
                                                self.get_start(mrna_loc), self.get_end(mrna_loc),
                                                strand, '.', mrna_id, gene_id)
                        mrna_exons = self.get_location_as_sorted_exons(mrna['location'], strand)
                        for exon in mrna_exons:
                            exon_id = 'exon_' + str(exon_id_generation)
                            exon_id_generation += 1
                            if is_gtf:
                                self.write_gtf_line(output, contig_id, 'exon', exon['start'],
                                                    exon['end'], strand, '.', gene_id, mrna_id)
                            else:
                                self.write_gff_line(output, contig_id, 'exon', exon['start'],
                                                    exon['end'], strand, '.', exon_id, mrna_id)
                        cds_exons = self.get_location_as_sorted_exons(cds['location'], strand)
                        cds_id = cds['id']
                        frame = 0
                        for exon in cds_exons:
                            f_start = exon['start']
                            f_end = exon['end']
                            f_length = f_end - f_start + 1
                            if is_gtf:
                                self.write_gtf_line(output, contig_id, 'CDS', f_start, 
                                                    f_end, strand, frame, gene_id, mrna_id)
                            else:
                                self.write_gff_line(output, contig_id, 'CDS', f_start, 
                                                    f_end, strand, frame, cds_id, mrna_id)
                            frame = (3 - ((f_length - frame) % 3)) % 3

        except Exception,e:
            raise ValueError("Failed to create file: {0}".format(e)) 
        finally:
            output.close()

        return {
            'file_path': str(out_file_path)
        }


    def get_location_as_sorted_exons(self, location_array, strand):
        ret = []
        for loc in location_array:
            ret.append({'start': self.get_start(loc), 'end': self.get_end(loc)})
        ret.sort(key = lambda item: item['start'], reverse = (strand != '+'))
        return ret


    def write_gff_line(self, output, contig_id, f_type, item_start, item_end, f_strand, frame,
                       item_id, parent_id):
        parent_suffix = (";Parent=" + parent_id) if parent_id else ""
        output.write(contig_id + "\tKBase\t" + f_type + "\t" + str(item_start) + "\t" + 
                     str(item_end) + "\t.\t" + f_strand + "\t"+ str(frame) + "\tID=" + item_id +
                     parent_suffix + "\n")


    def write_gtf_line(self, output, contig_id, f_type, item_start, item_end, f_strand, frame,
                       gene_id, trans_id):
        gene_id = gene_id if gene_id else ""
        trans_id = trans_id if trans_id else ""
        output.write(contig_id + "\tKBase\t" + f_type + "\t" + str(item_start) + "\t" + 
                     str(item_end) + "\t.\t" + f_strand + "\t"+ str(frame) + "\t" + 
                     "gene_id \"" + gene_id + "\"; transcript_id \"" + trans_id + "\";\n")


    def get_start(self, loc):
        start = loc[1]
        strand = loc[2]
        leng = loc[3]
        if strand == '+':
            return start
        if strand == '-':
            return start - ( leng - 1 )
        return 0


    # copied from RNASeq script_utils
    def get_end(self, loc):
        start = loc[1]
        strand = loc[2]
        leng = loc[3]
        if strand == '+':
            return start + ( leng - 1 )
        if strand == '-':
            return start
        return 0


    def get_common_location(self, location_array):
        contig = location_array[0][0]
        strand = location_array[0][2]
        min_pos = min([self.get_start(loc) for loc in location_array])
        max_pos = max([self.get_end(loc) for loc in location_array])
        common_length = max_pos - min_pos + 1
        common_start = min_pos if strand == '+' else max_pos
        return [[contig, common_start, strand, common_length]]


    def validate_params(self, params):
        if 'genome_ref' not in params:
            raise ValueError('required "genome_ref" field was not defined')


