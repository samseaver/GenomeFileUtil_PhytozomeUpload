
import os
import sys
import shutil
import traceback

from pprint import pprint, pformat

from biokbase.workspace.client import Workspace
from DataFileUtil.DataFileUtilClient import DataFileUtil

from GenomeAnnotationAPI.GenomeAnnotationAPIClient import GenomeAnnotationAPI


class GenomeToGFF:
    '''
    typedef structure {
        string genome_ref;
        list <string> ref_path_to_genome;
    } GenomeToGTFParams;

    /* from_cache is 1 if the file already exists and was just returned, 0 if
    the file was generated during this call. */
    typedef structure {
        File gtf_file;
        boolean from_cache;
    } GenomeToGTFResult;

    funcdef genome_to_gtf(GenomeToGTFParams params)
                returns (GenomeToGTFResult result) authentication required;
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

        # 4) if the GFF handle is there, get it and return
        print('checking if GFF file is cached...')
        result = self.get_gff_handle(data)
        if result is not None:
            result['from_cache'] = 1
            return result

        # 5) otherwise, build the GFF file and return it
        print('not cached, building file...')
        result = self.build_gff_file(getGenomeOptions, info[1])
        if result is None:
            raise ValueError('Unable to generate file.  Something went wrong')
        result['from_cache'] = 0
        return result


    def get_gff_handle(self, data):

        if 'gff_handle_ref' not in data:
            return None
        if data['gff_handle_ref'] is None:
            return None

        print('pulling cached GFF file from Shock: '+str(data['gff_handle_ref']))
        dfu = DataFileUtil(self.cfg.callbackURL)
        file = dfu.shock_to_file({
                            'shock_id':data['gff_handle_ref'],
                            'file_path':self.cfg.sharedFolder,
                            'unpack': 'unpack'
                        })
        return {
            'file_path': file['file_path']
        }


    ###  see logic from: https://github.com/kbase/KBaseRNASeq/blob/e2d69e4137903c68b5a1fedeab57f7900aad7253/lib/biokbase/RNASeq/KBaseRNASeqImpl.py#L443-L562
    def build_gff_file(self, getGenomeOptions, output_filename):

        # first get subdata needed; forget about the metadata
        getGenomeOptions['included_fields'] = []
        getGenomeOptions['included_feature_fields'] = ['id', 'type', 'location']
        getGenomeOptions['no_metadata'] = 1

        api = GenomeAnnotationAPI(self.cfg.callbackURL)
        genome_data = api.get_genome_v1(getGenomeOptions)['genomes'][0]['data']

        # create the file
        try:         
            out_file_path = os.path.join(self.cfg.sharedFolder, output_filename+'.gtf')
            print('Creating file: '+ str(out_file_path))
            output = open(out_file_path,'w')

            # write the file
            if 'features' in genome_data:
                for f in genome_data['features']:
                    if "type" in f and  f['type'] == 'CDS': f_type = f['type']
                    if "id" in f: f_id = f['id']
                    if "location" in f:
                        for contig_id,f_start,f_strand,f_len  in f['location']:
                            f_end = self.get_end(int(f_start),int(f_len),f_strand)
                            output.write(contig_id + "\tKBase\t" + f_type + "\t" + str(f_start) + "\t" + str(f_end) + "\t.\t" + f_strand + "\t"+ str(0) + "\ttranscript_id " + f_id + "; gene_id " + f_id + ";\n")

        except Exception,e:
            raise ValueError("Failed to create file: {0}".format(e)) 
        finally:
            output.close()

        return {
            'file_path': str(out_file_path)
        }

    # copied from RNASeq script_utils
    def get_end(self, start, leng, strand):
        stop = 0
        if strand == '+':
            stop = start + ( leng - 1 )
        if strand == '-':
            stop = start - ( leng + 1)
        return stop


    def validate_params(self, params):
        if 'genome_ref' not in params:
            raise ValueError('required "genome_ref" field was not defined')


