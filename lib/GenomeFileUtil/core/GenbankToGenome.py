
import os
import sys
import shutil
import traceback
import uuid
import urllib2

from urlparse import urlparse
from pprint import pprint, pformat

from biokbase.workspace.client import Workspace
from DataFileUtil.DataFileUtilClient import DataFileUtil

from GenomeFileUtil.core.GenbankUploaderScript import upload_genome

class GenbankToGenome:

    def __init__(self, sdk_config):
        self.cfg = sdk_config;

    def import_file(self, ctx, params):

        # 1) validate parameters and extract defaults
        self.validate_params(params)

        # 2) construct the input directory staging area
        input_directory = self.stage_input(params)


        # 3) extract out the parameters
        workspace_name = params['workspace_name']
        genome_name = params['genome_name']
        source = 'Genbank'
        if 'source' in params:
            source = source;
        taxon_wsname = 'ReferenceTaxons'
        if 'taxon_wsname' in params:
            taxon_wsname = params['taxon_wsname']

        if 'release' in params:
            release = params['release']

        if 'genetic_code' in params:
            genetic_code = params['genetic_code']

        # other options to handle
        # release
        # taxon_reference
        # exclude_feature_types
        # type

        # 4) Do the upload
        info = upload_genome(
                logger=None,
        
                shock_service_url = self.cfg.shockURL,
                handle_service_url = self.cfg.handleURL,
                workspace_service_url = self.cfg.workspaceURL,
        
                input_directory=input_directory,
        
                workspace_name   = workspace_name,
                core_genome_name = genome_name,
                source           = source,
                taxon_wsname     = taxon_wsname,
                release          = release,
                genetic_code     = genetic_code,

                provenance = ctx['provenance']
            )

        # 5) clear the temp directory
        shutil.rmtree(input_directory)

        # 6) return the result
        details = {
            'genome_ref': str(info[6]) + '/' + str(info[0]) + '/' + str(info[4]),
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
            raise ValueError('required "file" field must include one source: path | shock_id | ftp_url')
        if n_valid_fields > 1:
            raise ValueError('required "file" field has too many sources specified: ' + str(file.keys()))



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
            print('Downloading file from SHOCK node: ' + str(self.cfg.shockURL) + ' - ' + str(file['shock_id']))
            sys.stdout.flush()
            dfUtil = DataFileUtil(self.cfg.callbackURL)
            file_name = dfUtil.shock_to_file({
                                    'file_path': input_directory,
                                    'shock_id': file['shock_id']
                                })['node_file_name']
            genbank_file_path = os.path.join(input_directory, file_name)

        if 'ftp_url' in file and file['ftp_url'] is not None:
            # Note that the Transform originally had a script_utils.download_from_urls method
            # that, if the url is a folder, pulls all subfiles.  That code recently broke when
            # fetching from NCBI (not clear if it is our issue or NCBI), but for now just
            # support the most common case- an FTP to a single file.
            print('Downloading file from: ' + str(file['ftp_url']))
            sys.stdout.flush()

            url = urlparse(file['ftp_url'])
            if url.scheme != 'ftp' and url.scheme != 'http':
                raise ValueError('Only FTP/HTTP servers are supported')
            file_name = 'genome.gbk'
            if url.path != '':
                file_name = url.path.split('/')[-1]

            req = urllib2.Request(file['ftp_url'])
            response = urllib2.urlopen(req)
            file_data = response.read()

            genbank_file_path = os.path.join(input_directory, file_name)
            with open(genbank_file_path, "w") as genbank_file:
                genbank_file.write(file_data)

        # extract the file if it is compressed
        if genbank_file_path is not None:
            print("staged input file =" + genbank_file_path)
            sys.stdout.flush()
            dfUtil = DataFileUtil(self.cfg.callbackURL)
            dfUtil.unpack_file({ 'file_path': genbank_file_path })

        else:
            raise ValueError('No valid files could be extracted based on the input')

        return input_directory




