/*

*/
module GenomeFileUtil {

    /* A boolean - 0 for false, 1 for true.
       @range (0, 1)
    */
    typedef int boolean;

    typedef structure {
        string path;
        string shock_id;
        string ftp_url;
    } File;

    /* */
    typedef structure {
        File file;

        string genome_name;
        string workspace_name;

        string source;
        string taxon_wsname;

    } GenbankToGenomeParams;

    typedef structure {
        string ref;
    } GenomeSaveResult;

    funcdef genbank_to_genome(GenbankToGenomeParams params)
                returns (GenomeSaveResult result) authentication required;
};
