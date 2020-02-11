This directory should contain scripts and files needed to test your module's code.

NOTE: make sure that your test environment is 'ci', the tests won't function otherwise.
 
Unfortunately not all tests can be run at once as it overloads the callback server.
So the tests have been broken up into three main areas.
The first is the \*test.py files under the "GenomeFileUtil/test" directory.
Those tests will be run with doing "kb-sdk test" from the command line in the repo.
To do the other tests, first move the current \*test.py files to "GenomeFileUtil/skip" directory
There are two other major groups of tests that can be done.  
One is located under the "GenomeFileUtil/skip/supplemental_genbank_tests" directory.
The other is located under the "GenomeFileUtil/skip/supplemental_gff_tests" directory.
Simply move the files from one of these two directories and put them into the
"GenomeFileUtil/test" directory. Then you can run "kb-sdk test" on them as normal.

## State of tests

### Core tests
	32 ran - all Pass!

### GFF supplemental tests
	37 ran - 2 Fails
		Fail - test_simple_fasta_gff_to_genome_w_null_params
			AssertionError: 'Unconfirmed Organism' != 'Unconfirmed Organism: unknown_taxon'
		Fail - test_fasta_gff_to_genome_json
			AssertionError: 'Unknown' != 'Eukaryota'

### Genbank tests
	65 ran - 1 Fail, 1 Error
		Error - user slebras does not have access to object 31767/5
		Fail  - test_upload (core.genbank_upload_parameter_test.MinimalGenbankUploadTest) - AssertionError: 'Unknown' != 'Eukaryota'

### Utility tests
	3 ran - 3 Errors
		Error - ModuleNotFoundError: No module named 'mock'
		Error - FileNotFoundError: [Errno 2] No such file or directory: 'data/mouse_chr1.gbff'
		Error - ModuleNotFoundError: No module named 'mock'

### Problematic tests
	Ran 9 tests in 38.080s - 4 Fails, 1 Error

	Error - test suite for <class 'core.save_genome_test.SaveGenomeTest'>
		ValueError: Configuration in <module>/test_local/test.cfg file should include second user credentials ('test_token2')

	======================================================================
	FAIL: test_for_empty_feature_warnings (core.gff_upload_quality_test.GenomeFileUtilTest)
	----------------------------------------------------------------------
	Traceback (most recent call last):
	  File "/kb/module/test/core/gff_upload_quality_test.py", line 194, in test_for_empty_feature_warnings
	    self.assertTrue(found_warning_count > 0, "No features had warnings.")
	AssertionError: False is not true : No features had warnings.

	======================================================================
	FAIL: test_for_empty_functions (core.gff_upload_quality_test.GenomeFileUtilTest)
	----------------------------------------------------------------------
	Traceback (most recent call last):
	  File "/kb/module/test/core/gff_upload_quality_test.py", line 143, in test_for_empty_functions
	    self.assertTrue(found_function_count > 0, "No features had functions.")
	AssertionError: False is not true : No features had functions.

	======================================================================
	FAIL: test_for_gene_synonyms (core.gff_upload_quality_test.GenomeFileUtilTest)
	----------------------------------------------------------------------
	Traceback (most recent call last):
	  File "/kb/module/test/core/gff_upload_quality_test.py", line 105, in test_for_gene_synonyms
	    self.assertTrue( found_synonyms, "Expected Gene Synonyms were not found")
	AssertionError: False is not true : Expected Gene Synonyms were not found

	======================================================================
	FAIL: test_getting_all_go_ontologies (core.gff_upload_quality_test.GenomeFileUtilTest)
	----------------------------------------------------------------------
	Traceback (most recent call last):
	  File "/kb/module/test/core/gff_upload_quality_test.py", line 166, in test_getting_all_go_ontologies
	    self.assertTrue(all_ontologies_accounted_for, "Not all expected ontologies were accounted for : " + str(check_all_go_ontologies))
	AssertionError: False is not true : Not all expected ontologies were accounted for : {'GO:0005737': 0, 'GO:0016563': 0, 'GO:0016564': 0, 'GO:0006350': 0}
			Error -
