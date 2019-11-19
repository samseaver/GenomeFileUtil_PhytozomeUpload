# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.9.0] - 2019-10-24

### Changed

- Use the Relation Engine API for finding and populating taxonomy data
- Use the NCBI taxonomy ID as the primary way to assign taxonomy data instead of sciname

## [0.8.10] - 2018-11

### Changed

- Use WSLargeDataIO for saving and pulling genomes
- Expand and refactor import of ontology ids

### Fixed

- Fix a few genome warnings

## [0.8.9] - 2018-06-04

### Fixed

- Correct genbank download of old genomes with contig sets and long contig IDs
- Make reference path though genome for assembly access to prevent permissions issues

## [0.8.8] - 2018-05-30

### Changed

- Use searchapi2 in production & appdev

## [0.8.7] - 2018-05-28

### Fixed

- Fix go term splitting for GFF uploads

### Added

- Add optional "upgrade" parameter to save_genome

## [0.8.6] - 2018-05-15

### Changed

- Genome object refactor. Baseline for this log
