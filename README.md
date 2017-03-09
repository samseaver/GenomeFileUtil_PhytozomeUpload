[![Build Status](https://travis-ci.org/msneddon/GenomeFileUtil.svg?branch=master)](https://travis-ci.org/msneddon/GenomeFileUtil)

# GenomeFileUtil
---

This is the basic readme for this module. Include any usage or deployment instructions and links to other documentation here.


KNOWN ISSUES: (attention developers) as of 3/9/2017 (jkbaumohl)
1) Currently the GenBank downloader will only include feature types of "gene" and "pseudogene", as
these are the features included in the "features" array of the GenomeObject. No mRNA and CDS features
will be included.

2) Per Paramvir's request features in the original GenBank file that were features of type "gene"
and had the "pseudogene" key in them were made into type "pseudogene" in the KBase GenomeObject.
So if someone uploaded a RefSeq GenBank to KBase, then downloaded it to GenBank format.
Then they reuploaded the KBase_derived GenBank file again to KBase then the now pseudogene
features would not be reuploaded (as the uploader only uploads "gene", "mRNA" and "CDS" feature types.)

3) The features in the KBase Derived downloaded GenBank file will be in an unexpected order to the user.
Typically a GenBank contig has the features sorted by location on the contig.
How Matt did the uploader it is sorted by "+" strand, then position along the contig, then
all the "-" strand features sorted by position on the contig.