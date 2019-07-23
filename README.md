[![Build Status](https://travis-ci.org/msneddon/GenomeFileUtil.svg?branch=master)](https://travis-ci.org/msneddon/GenomeFileUtil)

# GenomeFileUtil
---

This is the basic readme for this module. Include any usage or deployment instructions and links to other documentation here.


# Testing
Important testing information. The Testing environment for the GenomeFileUtil is 'ci'


##Release Notes

### 0.8.10 (11/2018)
Use WSLargeDataIO for saving and pulling genomes
Expand and refactor import of ontology ids
Fix a few genome warnings

### 0.8.9 (06/04/2018)
Correct genbank download of old genomes with contig sets and long contig IDs
Make reference path though genome for assembly access to prevent permissions issues

### 0.8.8 (05/30/2018)
use searchapi2 in production & appdev

### 0.8.7 (05/28/2018)
Fix go term splitting for GFF uploads
Add optional "upgrade" parameter to save_genome

### 0.8.6 (05/15/2018)
Genome object refactor. Baseline for this log