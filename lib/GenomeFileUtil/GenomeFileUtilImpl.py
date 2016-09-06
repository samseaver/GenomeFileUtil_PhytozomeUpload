# -*- coding: utf-8 -*-
#BEGIN_HEADER

import os
import sys
import shutil
import traceback
import uuid
from pprint import pprint, pformat

from GenomeFileUtil.GenbankToGenome import GenbankToGenome


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
    GIT_URL = ""
    GIT_COMMIT_HASH = ""
    
    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.workspaceURL = config['workspace-url']
        self.shockURL = config['shock-url']
        self.handleURL = config['handle-service-url']
        self.callbackURL = os.environ['SDK_CALLBACK_URL']
        self.sharedFolder = config['scratch']
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

        importer = GenbankToGenome(
                        self.workspaceURL,
                        self.shockURL,
                        self.handleURL,
                        self.callbackURL,
                        self.sharedFolder)

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

    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK", 'message': "", 'version': self.VERSION, 
                     'git_url': self.GIT_URL, 'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
