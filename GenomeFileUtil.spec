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
        string genome_ref;
    } GenomeSaveResult;

    funcdef genbank_to_genome(GenbankToGenomeParams params)
                returns (GenomeSaveResult result) authentication required;

    typedef structure {
        string genome_ref;
        list <string> ref_path_to_genome;
    } GenomeToGFFParams;

    /* from_cache is 1 if the file already exists and was just returned, 0 if
    the file was generated during this call. */
    typedef structure {
        File gff_file;
        boolean from_cache;
    } GenomeToGFFResult;

    funcdef genome_to_gff(GenomeToGFFParams params)
                returns (GenomeToGFFResult result) authentication required;

    typedef structure {
        string genome_ref;
        list <string> ref_path_to_genome;
    } GenomeToGenbankParams;

    /* from_cache is 1 if the file already exists and was just returned, 0 if
    the file was generated during this call. */
    typedef structure {
        File gff_file;
        boolean from_cache;
    } GenomeToGenbankResult;

    funcdef genome_to_genbank(GenomeToGenbankParams params)
                returns (GenomeToGenbankResult result) authentication required;
};
