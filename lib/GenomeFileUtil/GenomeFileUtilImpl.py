# -*- coding: utf-8 -*-
#BEGIN_HEADER

import os
import sys
import shutil
import traceback
import uuid
from pprint import pprint, pformat

from GenomeFileUtil.core.GenbankToGenome import GenbankToGenome
from GenomeFileUtil.core.GenomeToGFF import GenomeToGFF

# Used to store and pass around configuration URLs more easily
class SDKConfig:
    def __init__(self, config):
        self.workspaceURL = config['workspace-url']
        self.shockURL = config['shock-url']
        self.handleURL = config['handle-service-url']
        self.callbackURL = os.environ['SDK_CALLBACK_URL']
        self.sharedFolder = config['scratch']

#END_HEADER


class GenomeFileUtil:
    '''
    Module Name:
    GenomeFileUtil

    Module Description:
    
    '''

    ######## WARNING FOR GEVENT USERS #######
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    #########################################
    VERSION = "0.0.1"
    GIT_URL = "git@github.com:kbaseapps/GenomeFileUtil.git"
    GIT_COMMIT_HASH = "431de8850afbb007b83e0c9fb197c103a974cb00"
    
    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.cfg = SDKConfig(config)
        #END_CONSTRUCTOR
        pass
    

    def genbank_to_genome(self, ctx, params):
        """
        :param params: instance of type "GenbankToGenomeParams" -> structure:
           parameter "file" of type "File" -> structure: parameter "path" of
           String, parameter "shock_id" of String, parameter "ftp_url" of
           String, parameter "genome_name" of String, parameter
           "workspace_name" of String, parameter "source" of String,
           parameter "taxon_wsname" of String
        :returns: instance of type "GenomeSaveResult" -> structure: parameter
           "ref" of String
        """
        # ctx is the context object
        # return variables are: result
        #BEGIN genbank_to_genome
        print('genbank_to_genome -- paramaters = ')
        pprint(params)

        importer = GenbankToGenome(self.cfg)
        result = importer.import_file(ctx, params)

        print('import complete -- result = ')
        pprint(result)
        #END genbank_to_genome

        # At some point might do deeper type checking...
        if not isinstance(result, dict):
            raise ValueError('Method genbank_to_genome return value ' +
                             'result is not type dict as required.')
        # return the results
        return [result]

    def genome_to_gff(self, ctx, params):
        """
        :param params: instance of type "GenomeToGFFParams" -> structure:
           parameter "genome_ref" of String, parameter "ref_path_to_genome"
           of list of String
        :returns: instance of type "GenomeToGFFResult" (from_cache is 1 if
           the file already exists and was just returned, 0 if the file was
           generated during this call.) -> structure: parameter "gff_file" of
           type "File" -> structure: parameter "path" of String, parameter
           "shock_id" of String, parameter "ftp_url" of String, parameter
           "from_cache" of type "boolean" (A boolean - 0 for false, 1 for
           true. @range (0, 1))
        """
        # ctx is the context object
        # return variables are: result
        #BEGIN genome_to_gff
        print('genome_to_gff -- paramaters = ')
        pprint(params)

        exporter = GenomeToGFF(self.cfg)
        result = exporter.export(ctx, params)

        print('export complete -- result = ')
        pprint(result)
        #END genome_to_gff

        # At some point might do deeper type checking...
        if not isinstance(result, dict):
            raise ValueError('Method genome_to_gff return value ' +
                             'result is not type dict as required.')
        # return the results
        return [result]

    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK", 'message': "", 'version': self.VERSION, 
                     'git_url': self.GIT_URL, 'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
