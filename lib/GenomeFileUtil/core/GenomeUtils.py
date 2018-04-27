from itertools import izip_longest

warnings = {
    "cds_excluded": "SUSPECT: CDS from {} was excluded because the associated "
                    "CDS failed coordinates validation",
    "cds_mrna_cds": "Feature order suggests that {} is the parent mRNA, but it"
                    " fails location validation",
    "cds_mrna_mrna": "Potential child CDS relationship failed due to location "
                    "validation.",
    "child_cds_failed": "The child CDS failed location validation. That CDS "
                        "has been excluded.",
    "child_mrna_failed": "The child mRNA failed location validation. That mRNA"
                         " has been excluded.",
    "genome_excluded": "SUSPECT: gene {} had some of its child features "
                       "(CDS and/or mRNAs) excluded because of failed "
                       "coordinates validation",
    "gene_excluded": "SUSPECT: gene {} was excluded because the associated CDS "
                     "failed coordinates validation",
    "mrna_excluded": "SUSPECT: mRNA from {} was excluded because the associated "
                     "mRNA failed coordinates validation",
    "no_spoof": "Some CDS features in the file do not have a parent gene. "
                "Ensure the correct file source is selected, correct the source file "
                "or select the 'generate_missing_genes' option.",
    "spoofed_gene": "This gene was not in the source GenBank file. It was "
                    "added to be the parent of the CDS {}.",
    "spoofed_genome": "SUSPECT: This genome has {} genes that needed to be "
                      "spoofed for existing parentless CDS.",
    "not_trans_spliced": "The feature coordinates order are suspect and the "
                         "feature is not flagged as being trans-spliced",
    "genome_not_trans_spliced": "SUSPECT: This Genome has {} features with "
                                "coordinates that are out of order and are "
                                "not trans_splicing.",
    "inconsistent_CDS_length": "This CDS has a length of {} which is not "
                               "consistent with the length of the translation "
                               "included ({} amino acids).",
    "genome_inc_CDS_length": "SUSPECT: CDS {} has a length of {} which is "
                             "not consistent with the length of the "
                             "translation included ({} amino acids).",
    "inconsistent_translation": "The annotated protein translation is not "
                                "consistent with the recorded DNA sequence.",
    "genome_inc_translation": "SUSPECT: This Genome has a high proportion "
                              "({} out of {}) CDS features that do not "
                              "translate the supplied translation.",
    "no_translation_supplied": "This CDS did not have a supplied "
                               "translation. The translation is derived "
                               "directly from DNA sequence.",
    "coordinates_off_end": "SUSPECT: Feature {} has invalid coordinates off "
                           "of the end of the contig and was not included.",
    "non_exact_coordinates": "The coordinates supplied for this feature are "
                             "non-exact. DNA or protein translations are "
                             "approximate.",
    "not_multiple_of_3CDS": "Sequence length {} is not a multiple of three",
    "non_standard_start_codon": "First codon '{}' is not a start codon",
    "out_of_order": "The feature coordinates order are out of order. GFF typically does not "
                    "designate trans_splicing.",
    "both_strand_coordinates": "The feature coordinates are both strands. GFF typically does not "
                    "designate trans_splicing.",
    "premature_stop_codon": "Extra in frame stop codon found.",
    "mRNA_fail_parent_coordinate_validation": "This mRNA lists CDS {} as its "
                    "corresponding CDS, but it fails coordinate validation.",
    "CDS_fail_child_of_mRNA_coordinate_validation": "This CDS lists mRNA {} as its "
                    "corresponding mRNA, but it fails coordinate validation.",
    "CDS_fail_child_of_gene_coordinate_validation": "This CDS lists gene {} as its "
                    "corresponding gene, but it fails coordinate validation.",
    "genes_mRNA_child_fails_location_validation": "The mRNA {} lists this gene as its "
                    "corresponding parent gene, but it fails coordinate validation.",
    "genes_CDS_child_fails_location_validation": "The CDS {} lists this gene as a "
                    "corresponding parent gene, but it fails coordinate validation.",
    "mRNAs_parent_gene_fails_location_validation": "This mRNA lists gene {} as its "
                    "corresponding parent gene, but it fails coordinate validation.",
    "generic_parents_child_fails_location_validation": "This feature lists feature {} as its "
                    "corresponding child, but it fails coordinate validation.",
    "generic_childs_parent_fails_location_validation": "This feature lists feature {} as its "
                    "corresponding parent, but it fails coordinate validation.",
    "gff_odd_strand_type": "This feature had \"{}\" as the strand designation and not + or -. "
                    "The location and sequence was defaulted to the + strand.",
    "contig_length_feature": "This feature spans entire contig length."
}


def get_start(loc):
    start = loc[1]
    strand = loc[2]
    leng = loc[3]
    if strand == '+':
        return start
    if strand == '-':
        return start - (leng - 1)
    return 0


def get_end(loc):
    start = loc[1]
    strand = loc[2]
    leng = loc[3]
    if strand == '+':
        return start + (leng - 1)
    if strand == '-':
        return start
    return 0


def get_bio_end(loc):
    if loc[2] == "+":
        return loc[1] + loc[3]
    else:
        return loc[1] - loc[3]

def is_parent(feat1, feat2):
    """Check if all locations in feat2 fall within a location in
    feat1"""

    def _contains(loc1, loc2):
        if loc1[2] != loc2[2]:  # different strands
            return False
        elif loc1[2] == "+":
            return loc2[1] >= loc1[1] and (
                loc2[1] + loc2[3] <= loc1[1] + loc1[3])
        else:
            return loc2[1] <= loc1[1] and (
                loc2[1] - loc2[3] >= loc1[1] - loc1[3])

    if feat1.get('type') == 'gene':
        for l2 in feat2['location']:
            if not any(_contains(l1, l2) for l1 in feat1['location']):
                return False
    else:
        # for a mrna, the first and last part match loosely (like a gene) but
        # the internal coordinates must be an exact match
        if not any(_contains(l1, feat2['location'][0]) for l1 in feat1['location']):
            return False

        if not any(_contains(l1, feat2['location'][-1]) for l1 in feat1['location']):
            return False

        if len(feat2['location']) > 1:
            if get_bio_end(feat2['location'][0]) != get_bio_end(feat1['location'][0]):
                return False

            if feat2['location'][-1][1] != feat1['location'][-1][1]:
                return False

            for l1, l2 in izip_longest(feat1['location'][1:-1],
                                       feat2['location'][1:-1]):
                if l1 != l2:
                    return False
    return True


def parse_inferences(inferences):
    """Whoever designed the genbank delimitation is an idiot: starts and
    ends with a optional values and uses a delimiter ":" that is
    used to divide it's DBs in the evidence. Anyway, this sorts that"""
    result = []
    for inf in inferences:
        try:
            sp_inf = inf.split(":")
            if sp_inf[0] in ('COORDINATES', 'DESCRIPTION', 'EXISTENCE'):
                inference = {'category': sp_inf.pop(0)}
            else:
                inference = {'category': ''}
            inference['type'] = sp_inf[0]
            inference['evidence'] = ":".join(sp_inf[1:])
            result.append(inference)
        except IndexError('Unparseable inference string: ' + inf):
            continue
    return result


def propagate_cds_props_to_gene(cds, gene):
    # Put longest protein_translation to gene
    if "protein_translation" not in gene or (
                len(gene["protein_translation"]) <
                len(cds["protein_translation"])):
        gene["protein_translation"] = cds["protein_translation"]
        gene["protein_translation_length"] = len(
            cds["protein_translation"])
    # Merge cds list attributes with gene
    for key in ('functions', 'aliases', 'db_xref'):
        if cds.get(key, []):
            gene[key] = cds.get(key, []) + gene.get(key, [])
    # Merge cds["ontology_terms"] -> gene["ontology_terms"]
    terms2 = cds.get("ontology_terms")
    if terms2 is not None:
        terms = gene.get("ontology_terms")
        if terms is None:
            gene["ontology_terms"] = terms2
        else:
            for source in terms2:
                if source in terms:
                    terms[source].update(terms2[source])
                else:
                    terms[source] = terms2[source]

def check_full_contig_length_or_multi_strand_feature(feature, is_transpliced, contig_length, skip_types):
    ''' 
    Tests for full contig length features and if on both strands.
    '''
    feature_min_location = None
    feature_max_location = None
    strand_set = set()
    contig_id = feature["location"][0][0]
    for location in feature["location"]:
        if location[0] != contig_id:
            return feature
        location_min = get_start(location)
        location_max = get_end(location)
        strand_set.add(location[2])                 
        if feature_min_location is None or feature_min_location > location_min:
            feature_min_location = location_min
        if feature_max_location is None or feature_max_location < location_max:
            feature_max_location = location_max
    if feature_min_location == 1 \
        and feature_max_location == contig_length \
        and feature['type'] not in skip_types: 
        feature["warnings"] = feature.get('warnings', []) + [warnings["contig_length_feature"]]  
    if len(strand_set) > 1 and not is_transpliced:
        feature["warnings"] = feature.get('warnings', []) + [warnings["both_strand_coordinates"]]
    return feature 
