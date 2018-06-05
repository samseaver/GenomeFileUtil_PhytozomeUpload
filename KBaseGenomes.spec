/*
@author chenry,kkeller
*/
module KBaseGenomes {
    /*
    	Reference to a ContigSet object containing the contigs for this genome in the workspace
		@id ws KBaseGenomes.ContigSet
	*/
    typedef string ContigSet_ref;
    /*
		Reference to a ProteinSet object containing the proteins for this genome in the workspace
		@id ws KBaseGenomes.ProteinSet
	*/
    typedef string ProteinSet_ref;
    /*
		Reference to a TranscriptSet object containing the transcripts for this genome in the workspace
		@id ws KBaseGenomes.TranscriptSet
	*/
    typedef string TranscriptSet_ref;
    /*
		Reference to a Feature object of a genome in the workspace
		@id subws KBaseGenomes.Genome.features.[*].id
	*/
    typedef string Feature_ref;
    /*
		Reference to a Genome object in the workspace
		@id ws KBaseGenomes.Genome KBaseGenomeAnnotations.GenomeAnnotation
	*/
    typedef string Genome_ref;
    /*
		Reference to an Assembly object in the workspace
		@id ws KBaseGenomeAnnotations.Assembly
	*/
    typedef string Assembly_ref;
    /*
		Reference to a Pangenome object in the workspace
		@id ws KBaseGenomes.Pangenome
	*/
    typedef string Pangenome_ref;
    /*
		Reference to a Proteome Comparison object in the workspace
		@id ws GenomeComparison.ProteomeComparison
	*/
    typedef string Protcomp_ref;
    /*
		Reference to a source_id
		@id external
	*/
    typedef string source_id;
    /*
		KBase genome ID
		@id kb
	*/
    typedef string Genome_id;
    /*
		KBase Reaction ID
		@id external
	*/
	typedef string Reaction_id;
    /*
		KBase Feature ID
		@id external
	*/
	typedef string Feature_id;
    /*
		KBase ProteinSet ID
		@id kb
	*/
    typedef string ProteinSet_id;
	/*
		ProbabilisticAnnotation ID
		@id kb
	*/
	typedef string ProbabilisticAnnotation_id;
	/*
		Genome protein ID
		@id external
	*/
    typedef string Protein_id;
    /*
		Reference to an individual contig in a ContigSet object
		@id subws KBase.ContigSet.contigs.[*].id
	*/
    typedef string Contig_ref;
    /*
		ContigSet contig ID
		@id external
	*/
    typedef string Contig_id;
    /*
		KBase contig set ID
		@id kb
	*/
    typedef string ContigSet_id;

    /*
		Reference to a reads file in shock
		@id shock
	*/
    typedef string Reads_ref;
    /*
		Reference to a fasta file in shock
		@id shock
	*/
    typedef string Fasta_ref;

    typedef string Feature_type;
    typedef int Bool;

    /* Type spec for a "Contig" subobject in the "ContigSet" object

		Contig_id id - ID of contig in contigset
		string md5 - unique hash of contig sequence
		string sequence - sequence of the contig
		string description - Description of the contig (e.g. everything after the ID in a FASTA file)

		@optional length md5 genetic_code cell_compartment replicon_geometry replicon_type name description complete
	*/
	typedef structure {
		Contig_id id;
		int length;
		string md5;
		string sequence;/*using "sequence" instead of "dna"*/
		int genetic_code;
		string cell_compartment;
		string replicon_type;
		/* circular / linear */
		string replicon_geometry;
		string name;
		string description;
		Bool complete;
    } Contig;

    /* Type spec for the "ContigSet" object

		contigset_id id - unique kbase ID of the contig set
		string name - name of the contig set
		string type - type of the contig set (values are: Genome,Transcripts,Environment,Collection)
		source_id source_id - source ID of the contig set
		string source - source of the contig set
		list<Contig> contigs - list of contigs in the contig set
		reads_ref reads_ref - reference to the shocknode with the rawreads from which contigs were assembled
		fasta_ref fasta_ref - reference to fasta file from which contig set were read

		@optional name type reads_ref fasta_ref
    	@metadata ws type as Type
		@metadata ws source_id as Source ID
		@metadata ws source as Source
		@metadata ws name as Name
		@metadata ws length(contigs) as Number contigs
	*/
	typedef structure {
		ContigSet_id id;
		string name;
		string md5;
		source_id source_id;
		string source;
		string type;
		Reads_ref reads_ref;
		Fasta_ref fasta_ref;
		list<Contig> contigs;
    } ContigSet;

    /*
    Structure for a publication
    (float pubmedid
    string source (ex. Pubmed)
    string title
    string web address
    string  publication year
    string authors
    string journal)
    */
    typedef tuple<float pubmedid,string source,string title, string url,string year,string authors, string journal> publication;

    /*
    KBase CDS ID
    @id external
    */
    typedef string cds_id;

    /*
    KBase mRNA ID
    @id external
    */
    typedef string mrna_id;

        /*
        category;#Maybe a controlled vocabulary
    type;#Maybe a controlled vocabulary
    */
    typedef structure {
      string category;
      string type;
      string evidence;
    } InferenceInfo;
    /*
    Structure for a single CDS encoding “gene” of a genome
    ONLY PUT GENES THAT HAVE A CORRESPONDING CDS IN THIS ARRAY

        NOTE: Sequence is optional. Ideally we can keep it in here, but
    Recognize due to space constraints another solution may be needed.
        We may want to add additional fields for other CDM functions
        (e.g., atomic regulons, coexpressed fids, co_occurring fids,...)

        protein_translation_length and protein_translation are
        for longest coded protein (representative protein for splice variants)

        NOTE: New Aliases field definitely breaks compatibility.
              As Does Function.

        flags are flag fields in GenBank format. This will be a controlled vocabulary.
        Initially Acceptable values are pseudo, ribosomal_slippage, and trans_splicing

        Md5 is the md5 of dna_sequence.

        @optional functions ontology_terms note protein_translation mrnas flags warnings
        @optional inference_data dna_sequence aliases db_xrefs children functional_descriptions
    */
    typedef structure {
      Feature_id id;
      list<tuple<Contig_id,int,string,int>> location;
      list<string> functions;
      list<string> functional_descriptions;
      mapping<string ontology_namespace,mapping<string ontology_id,list<int> evidence_events>> ontology_terms;
      string note;
      string md5;
      string protein_translation;
      int protein_translation_length;
      list<string> cdss;
      list<string> mrnas;
      list<string> children;
      list<string> flags;
      list<string> warnings;
      list <InferenceInfo> inference_data;
      string dna_sequence;
      int dna_sequence_length;
      list<tuple<string fieldname,string alias>> aliases;
      list<tuple<string db_source,string db_identifier>> db_xrefs;
    } Feature;

        /*
          Structure for a single feature that is NOT one of the following:
          Protein encoding gene (gene that has a corresponding CDS)
          mRNA
          CDS

          Note pseudo-genes and Non protein encoding genes are put into this

      flags are flag fields in GenBank format. This will be a controlled vocabulary.
        Initially Acceptable values are pseudo, ribosomal_slippage, and trans_splicing
        Md5 is the md5 of dna_sequence.

        @optional functions ontology_terms note flags warnings functional_descriptions
        @optional inference_data dna_sequence aliases db_xrefs children parent_gene
    */
    typedef structure {
      Feature_id id;
      list<tuple<Contig_id,int,string,int>> location;
      string type;
      list<string> functions;
      list<string> functional_descriptions;
      mapping<string ontology_namespace,mapping<string ontology_id,list<int> evidence_event>> ontology_terms;
      string note;
      string md5;
      string parent_gene;
      list<string> children;
      list<string> flags;
      list<string> warnings;
      list <InferenceInfo> inference_data;
      string dna_sequence;
      int dna_sequence_length;
      list<tuple<string fieldname,string alias>> aliases;
      list<tuple<string db_source,string db_identifier>> db_xrefs;
    } NonCodingFeature;


        /*
          Structure for a single feature CDS

      flags are flag fields in GenBank format. This will be a controlled vocabulary.
        Initially Acceptable values are pseudo, ribosomal_slippage, and trans_splicing
        Md5 is the md5 of dna_sequence.

        @optional parent_mrna functions ontology_terms note flags warnings
        @optional inference_data dna_sequence aliases db_xrefs functional_descriptions
    */
    typedef structure {
      cds_id id;
      list<tuple<Contig_id,int,string,int>> location;
      string md5;
      string protein_md5;
      Feature_id parent_gene;
      mrna_id parent_mrna;
      string note;
      list<string> functions;
      list<string> functional_descriptions;
      mapping<string ontology_namespace,mapping<string ontology_id,list<int> evidence_events>> ontology_terms;
             list<string> flags;
      list<string> warnings;
      list <InferenceInfo> inference_data;
      string protein_translation;
      int protein_translation_length;
      list<tuple<string fieldname,string alias>> aliases;
      list<tuple<string db_source,string db_identifier>> db_xrefs;
      string dna_sequence;
      int dna_sequence_length;
    } CDS;


        /*
          Structure for a single feature mRNA

      flags are flag fields in GenBank format. This will be a controlled vocabulary.
        Initially Acceptable values are pseudo, ribosomal_slippage, and trans_splicing
        Md5 is the md5 of dna_sequence.

        @optional cds functions ontology_terms note flags warnings
        @optional inference_data dna_sequence aliases db_xrefs functional_descriptions
    */
    typedef structure {
      mrna_id id;
      list<tuple<Contig_id,int,string,int>> location;
      string md5;
      Feature_id parent_gene;
      cds_id cds;
      string dna_sequence;
      int dna_sequence_length;
      string note;
      list<string> functions;
      list<string> functional_descriptions;
      mapping<string ontology_namespace,mapping<string ontology_id,list<int> evidence_events>> ontology_terms;
             list<string> flags;
      list<string> warnings;
      list <InferenceInfo> inference_data;
      list<tuple<string fieldname,string alias>> aliases;
      list<tuple<string db_source,string db_identifier>> db_xrefs;
    } mRNA;

    /*
    Reference to a taxon object
        @id ws KBaseGenomeAnnotations.Taxon
    */
    typedef string Taxon_ref;

    /*
    Reference to a handle to the Genbank file on shock
        @id handle
    */
    typedef string genbank_handle_ref;

    /*
    Reference to a handle to the GFF file on shock
        @id handle
    */
    typedef string gff_handle_ref;

    /*
    Reference to a ontology object
        @id ws KBaseOntology.OntologyDictionary
    */
    typedef string Ontology_ref;

    /*
    Reference to a report object
        @id ws KBaseReport.Report
    */
    typedef string Method_report_ref;

    /*
        @optional ontology_ref method_version eco
    */
    typedef structure {
        string id;
        Ontology_ref ontology_ref;
        string method;
        string method_version;
        string timestamp;
        string eco;
    } Ontology_event;

    /*
    Score_interpretation : fraction_complete - controlled vocabulary managed by API
        @optional method_report_ref method_version
    */
    typedef structure {
        string method;
    Method_report_ref method_report_ref;
        string method_version;
        string score;
        string score_interpretation;
    string timestamp;
    } GenomeQualityScore;


    /*
    Genome object holds much of the data relevant for a genome in KBase
        Genome publications should be papers about the genome
    Should the Genome object contain a list of contig_ids too?
    Source: allowed entries RefSeq, Ensembl, Phytozome, RAST, Prokka, User_upload
        #allowed entries RefSeq, Ensembl, Phytozome, RAST, Prokka,
    User_upload controlled vocabulary managed by API

    Domain is a controlled vocabulary
    Warnings : mostly controlled vocab but also allow for unstructured
    Genome_tiers : controlled vocabulary (based on ap input and API checked)
    Allowed values: #Representative, Reference, ExternalDB, User

    Examples Tiers:
    All phytozome - Representative and ExternalDB
    Phytozome flagship genomes - Reference, Representative and ExternalDB
    Ensembl - Representative and ExternalDB
    RefSeq Reference - Reference, Representative and ExternalDB
    RefSeq Representative - Representative and ExternalDB
    RefSeq Latest or All Assemblies folder - ExternalDB
    User Data - User tagged

    Example Sources:
    RefSeq, Ensembl, Phytozome, Microcosm, User, RAST, Prokka, (other annotators)


    @optional warnings contig_lengths contig_ids source_id taxonomy publications
    @optional ontology_events ontologies_present non_coding_features mrnas
    @optional genbank_handle_ref gff_handle_ref external_source_origination_date
    @optional release original_source_file_name notes quality_scores suspect assembly_ref


    @metadata ws gc_content as GC content
        @metadata ws taxonomy as Taxonomy
        @metadata ws md5 as MD5
        @metadata ws dna_size as Size
        @metadata ws genetic_code as Genetic code
        @metadata ws domain as Domain
        @metadata ws source_id as Source ID
        @metadata ws source as Source
        @metadata ws scientific_name as Name
        @metadata ws length(features) as Number of Protein Encoding Genes
    @metadata ws length(cdss) as Number of CDS
        @metadata ws assembly_ref as Assembly Object
    @metadata ws num_contigs as Number contigs
    */
    typedef structure {
      Genome_id id;
      string scientific_name;
      string domain;
      list<string> warnings;
      list<string> genome_tiers;
    mapping<string type, int count> feature_counts;
      int genetic_code;
      int dna_size;
      int num_contigs;
      string molecule_type;
      list<int> contig_lengths;
      list<string> contig_ids;
      string source;
      source_id source_id;
      string md5;
      string taxonomy;
      float gc_content;
      list<publication> publications;
      list<Ontology_event> ontology_events;
      mapping<string ontology_namespace, mapping<string ontology_id,string termname>> ontologies_present;
      list<Feature> features;
      list<NonCodingFeature> non_coding_features;
      list<CDS> cdss;
      list<mRNA> mrnas;
      Assembly_ref assembly_ref;
      Taxon_ref taxon_ref;
      genbank_handle_ref genbank_handle_ref;
      gff_handle_ref gff_handle_ref;
      string external_source_origination_date;
      string release;
      string original_source_file_name;
      string notes;
      list<GenomeQualityScore> quality_scores;
      Bool suspect;
    } Genome;

    /*
        Structure for a protein family

		@optional query_begin query_end subject_begin subject_end score evalue subject_description release_version
    */
    typedef structure {
		string id;
		string subject_db;
		string release_version;
		string subject_description;
		int query_begin;
		int query_end;
		int subject_begin;
		int subject_end;
		float score;
		float evalue;
    } ProteinFamily;

    typedef tuple<string comment, string annotator, float annotation_time> annotation;

	/* Type spec for the "Protein" object

		Protein_id id - unique external ID of protein
		string function - annotated function for protein
		string md5 - md5 hash of protein sequence
		string sequence - amino acid sequence of protein
		int length - length of protein
		list<ProteinFamily> protein_families - families to which the protein belongs
		list<string> aliases - aliases for the protein
		list<annotation> annotations - curator annotations on protein
		list<subsystem_data> subsystem_data;

		@optional function
    	@searchable ws_subset id md5 function length aliases
	*/
	typedef structure {
		Protein_id id;
		string function;
		string md5;
		string sequence;
		int length;
		list<ProteinFamily> protein_families;
		list<string> aliases;
		list<annotation> annotations;
    } Protein;

   /* Type spec for the "ProteinSet" object

		proteinset_id id - unique kbase ID of the protein set
		string name - name of the protein set
		string type - type of the protein set (values are: Organism,Environment,Collection)
		source_id source_id - source ID of the protein set
		string source - source of the protein set
		list<Protein> proteins - list of proteins in the protein set
		fasta_ref fasta_ref - reference to fasta file from which contig set were read

		@optional name type fasta_ref
    	@searchable ws_subset proteins.[*].(id,md5,function,length,aliases) md5 id name source_id source type
	*/
	typedef structure {
		ProteinSet_id id;
		string name;
		string md5;
		source_id source_id;
		string source;
		string type;
		Fasta_ref fasta_ref;
		list<Protein> proteins;
    } ProteinSet;

    /*
       A function_probability is a (annotation, probability) pair associated with a gene
       An annotation is a "///"-delimited list of roles that could be associated with that gene.
    */
    typedef tuple<string annotation, float probability> function_probability;

    /* Object to carry alternative functions and probabilities for genes in a genome

        probanno_id id - ID of the probabilistic annotation object
        Genome_ref genome_ref - reference to genome probabilistic annotation was built for
        mapping<Feature_id, list<function_probability>> roleset_probabilities - mapping of features to list of alternative function_probability objects
        list<Feature_id> skipped_features - list of features in genome with no probability

    	@searchable ws_subset id genome_ref skipped_features

    */
    typedef structure {
		ProbabilisticAnnotation_id id;
		Genome_ref genome_ref;
		mapping<Feature_id,list<function_probability>> roleset_probabilities;
		list<Feature_id> skipped_features;
    } ProbabilisticAnnotation;

    /* Structure for the "MetagenomeAnnotationOTUFunction" object

		list<string> reference_genes - list of genes associated with hit
		string functional_role - annotated function
		string kbid - kbase ID of OTU function in metagenome
		int abundance - number of hits with associated role and OTU
		float confidence - confidence of functional role hit
		string confidence_type - type of functional role hit

    	@searchable ws_subset id abundance confidence functional_role
	*/
    typedef structure {
		string id;
		list<string> reference_genes;
		string functional_role;
		int abundance;
		float confidence;
    } MetagenomeAnnotationOTUFunction;

    /* Structure for the "MetagenomeAnnotationOTU" object

		string name - name of metagenome OTU
		string kbid - KBase ID of OTU of metagenome object
		string source_id - ID used for OTU in metagenome source
		string source - source OTU ID
		list<MetagenomeAnnotationOTUFunction> functions - list of functions in OTU

    	@searchable ws_subset id name source_id source functions.[*].(id,abundance,confidence,functional_role)

	*/
    typedef structure {
    	float ave_confidence;
		float ave_coverage;
		string id;
		string name;
		string source_id;
		string source;
		list<MetagenomeAnnotationOTUFunction> functions;
    } MetagenomeAnnotationOTU;

    /* Structure for the "MetagenomeAnnotation" object

		string type - type of metagenome object
		string name - name of metagenome object
		string kbid - KBase ID of metagenome object
		string source_id - ID used in metagenome source
		string source - source of metagenome data
		string confidence_type - type of confidence score
		list<MetagenomeAnnotationOTU> otus - list of otus in metagenome

    	@searchable ws_subset type name id source_id source confidence_type otus.[*].(id,name,source_id,source,functions.[*].(id,abundance,confidence,functional_role))
		@metadata ws type as Type
		@metadata ws name as Name
		@metadata ws source_id as Source ID
		@metadata ws source as Source
		@metadata ws length(otus) as Number OTUs
	*/
    typedef structure {
		string type;
		string name;
		string id;
		string source_id;
		string source;
		string confidence_type;
		list<MetagenomeAnnotationOTU> otus;
    } MetagenomeAnnotation;

    /*
		Domain - a subobject holding information on a single protein domain
		string id - numerical ID assigned by KBase
		string source_id - assession ID from CDD database;
		string type - type of CDD, possible values are cd, pfam, smart, COG, PRK, CHL
		string name - name of CDD
		string description - description of CDD
    */
    typedef structure {
		string id;
		string source_id;
		string type;
		string name;
		string description;
    } Domain;

    /*
		FeatureDomain - a subobject holding information on how a domain appears in a gene
		string id - numerical ID assigned by KBase
		string source_id - assession ID from CDD database;
		string type - type of CDD, possible values are cd, pfam, smart, COG, PRK, CHL
		string name - name of CDD
		string description - description of CDD

		@optional feature_ref domains
    */
    typedef structure {
		string id;
		string feature_id;
		string feature_ref;
		string function;
		int feature_length;
		list<tuple<string domain_ref,int identity,int alignment_length,int mismatches,int gaps,float protein_start,float protein_end,float domain_start,float domain_end,float evalue,float bit_score>> domains;
    } FeatureDomainData;

    /*
    	GenomeDomainData object: this object holds all data regarding protein domains in a genome in KBase

		@optional genome_ref
    	@searchable ws_subset id genome_id scientific_name genome_ref num_domains num_features
    */
    typedef structure {
    	string id;
		Genome_id genome_id;
		string scientific_name;
		Genome_ref genome_ref;
		int num_domains;
		int num_features;

		list<Domain> domains;
		list<FeatureDomainData> featuredomains;
	} GenomeDomainData;

	/*
    	OrthologFamily object: this object holds all data for a single ortholog family in a metagenome

    	@optional type function md5 protein_translation
    */
	typedef structure {
    	string id;
		string type;
		string function;
		string md5;
		string protein_translation;
		list<tuple<string,float,string>> orthologs;
	} OrthologFamily;

	/*
    	Pangenome object: this object holds all data regarding a pangenome

    	@searchable ws_subset id name
		@metadata ws type as Type
		@metadata ws name as Name
		@metadata ws length(orthologs) as Number orthologs
		@metadata ws length(genome_refs) as Number genomes
    */
    typedef structure {
    	string id;
    	string name;
    	string type;
    	list<Genome_ref> genome_refs;
    	list<OrthologFamily> orthologs;
	} Pangenome;

	/*
    	GenomeComparisonGenome object: this object holds information about a genome in a genome comparison
    */
    typedef structure {
		string id;
		Genome_ref genome_ref;
		mapping<string genome_id,tuple<int commonfamilies,int commonfunctions> > genome_similarity;
		string name;
		string taxonomy;
		int features;
		int families;
		int functions;
    } GenomeComparisonGenome;

    /*
    	GenomeComparisonFunction object: this object holds information about a genome in a function across all genomes
    */
    typedef structure {
		int core;
		mapping<string genome_id,list<tuple<Feature_id,int famindex,float score> > > genome_features;
		string id;
		list<tuple<Reaction_id, string equation>> reactions;
		string subsystem;
		string primclass;
		string subclass;
		int number_genomes;
		float fraction_genomes;
		float fraction_consistent_families;
		string most_consistent_family;
    } GenomeComparisonFunction;

    /*
    	GenomeComparisonFamily object: this object holds information about a protein family across a set of genomes
    */
    typedef structure {
		int core;
		mapping<string genome_id,list< tuple<Feature_id,list<int> funcindecies,float score > > > genome_features;
		string id;
		string type;
		string protein_translation;
		int number_genomes;
		float fraction_genomes;
		float fraction_consistent_annotations;
		string most_consistent_role;
    } GenomeComparisonFamily;

    /*
    	GenomeComparisonData object: this object holds information about a multigenome comparison

    	@optional protcomp_ref pangenome_ref
    	@metadata ws core_functions as Core functions
		@metadata ws core_families as Core families
		@metadata ws name as Name
		@metadata ws length(genomes) as Number genomes
    */
    typedef structure {
		string id;
		string name;
		int core_functions;
		int core_families;
		Protcomp_ref protcomp_ref;
		Pangenome_ref pangenome_ref;
		list<GenomeComparisonGenome> genomes;
		list<GenomeComparisonFamily> families;
		list<GenomeComparisonFunction> functions;
    } GenomeComparison;
};