"""
Genome CDS Protein Sequence to Fasta file conversion.
"""

import time
from collections import defaultdict

from installed_clients.DataFileUtilClient import DataFileUtil


def log(message, prefix_newline=False):
    time_str = time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(time.time()))
    print(('\n' if prefix_newline else '') + time_str + ': ' + message)


class GenomeFeaturesProteinToFasta(object):

    def __init__(self, sdk_config):
        self.cfg = sdk_config
        self.dfu = DataFileUtil(self.cfg.callbackURL)


    def validate_params(self, params):
        if 'genome_ref' not in params:
            raise ValueError('required "genome_ref" field was not defined')
        
    def export(self, ctx, params):
        # 1) validate parameters and extract defaults
        self.validate_params(params)

        include_functions = False
        if 'include_functions' in params:
            include_functions = True
        include_aliases = False
        if 'include_aliases' in params:
            include_aliases = True

        # 2) get genome info
        genome_data = self.dfu.get_objects({
            'object_refs': [params['genome_ref']]
        })['data'][0]
        info = genome_data['info']
        data = genome_data['data']


        # 3) make sure the type is valid
        if info[2].split(".")[1].split('-')[0] != 'Genome':
            raise ValueError('Object is not a Genome, it is a:' + str(info[2]))

        # 4) make sure has CDSs
        if 'cdss' not in data:
            raise ValueError('This genome does not have a CDS section to it')
        elif len(data['cdss']) == 0:
            raise ValueError('This genome does not have any CDS features')
        else:
            cds_data = data['cdss']

        # 5) build the fasta file and return it
        log('not cached, building file...')
        result = self.build_protein_fasta_file(cds_data, info[1] + "_protein.faa",
                                         include_functions, include_aliases)
        if result is None:
            raise ValueError('Unable to generate file.  Something went wrong')
        result['from_cache'] = 0
        return result

    def build_protein_fasta_file(self, cds_data, output_filename, include_functions, include_aliases):
        file_path = self.cfg.sharedFolder + "/" + output_filename
        f = open("file_path", "w")
        for cds in cds_data:
            if len(cds["protein_translation"]) > 0:
                #Build up header line
                header_line = ">{}".format(cds["id"])
                if include_functions:
                    if "functions" in cds and len(cds["functions"]) > 0:
                        header_line += " Functions:{}".format(','.join(cds["functions"]))
                    if "functional_descriptions" in cds and len(cds["functional_descriptions"]) > 0:
                        header_line += " Functional_descriptions:{}".format(','.join(cds["functional_descriptions"]))
                if include_aliases:
                    if "aliases" in cds and len(cds["aliases"]) > 0:
                        alias_set = set()
                        for alias in cds["aliases"]:
                            alias_set.add[alias[1]]
                        header_line += " Aliases:{}".format(','.join(alias_set))
                    if "db_xrefs" in cds and len(cds["db_xrefs"]) > 0:
                        db_xref_info = "db_xrefs"
                        for db_xref in cds["db_xrefs"]:
                            db_xref_info += "{}:{},".format(db_xref[0],db_xref[1])
                        header_line += db_xref_info[:-1]
                f.write(header_line)
                f.write(cds["protein_translation"])
        f.close
        return {
            'protein_fasta_file': {
                'file_path': file_path
            }
        }