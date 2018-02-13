from itertools import izip_longest

warnings = {
    "cds_excluded": "SUSPECT CDS from {} was excluded because the associated "
                    "CDS failed coordinates validation",
    "cds_mrna_cds": "Feature order suggests that {} is the parent mRNA, but it"
                    " fails location validation",
    "cds_mrna_mrna": "Potential child CDS relationship failed due to location "
                    "validation.",
    "child_cds_failed": "The child CDS failed location validation. That CDS "
                        "has been excluded.",
    "child_mrna_failed": "The child mRNA failed location validation. That mRNA"
                         " has been excluded.",
    "genome_excluded": "SUSPECT gene {} had some of its child features "
                       "(CDS and/or mRNAs) excuded because of failed "
                       "coordinates validation",
    "gene_excluded": "SUSPECT gene {} was excluded because the associated CDS "
                     "failed coordinates validation",
    "mrna_excluded": "SUSPECT mRNA from {} was excluded because the associated "
                     "mRNA failed coordinates validation",
    "no_spoof": "Some CDS features in the file do not have a parent gene. "
                "Either fix the source file or select the "
                "'generate_missing_genes' option.",
    "spoofed_gene": "This gene was not in the source GenBank file. It was "
                    "added to be the parent of the CDS {}.",
    "spoofed_genome": "SUSPECT this genome has {} genes that needed to be "
                      "spoofed for existing parentless CDS."
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