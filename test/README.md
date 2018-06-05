This directory should contain scripts and files needed to test your module's code.
 
Unfortunately not all tests can be run at once as it overloads the callback server.
So the tests have been broken up into three main areas.
The first is the *test.py files under the "GenomeFileUtil/test" directory.
Those tests will be run with doing "kb-sdk test" from the command line in the repo.
To do the other tests, first move the current *test.py files to "GenomeFileUtil/skip" directory
There are two other major groups of tests that can be done.  
One is located under the "GenomeFileUtil/skip/supplemental_genbank_tests" directory.
The other is located under the "GenomeFileUtil/skip/supplemental_gff_tests" directory.
Simply move the files from one of these two directories and put them into the
"GenomeFileUtil/test" directory. Then you can run "kb-sdk test" on them as normal.
