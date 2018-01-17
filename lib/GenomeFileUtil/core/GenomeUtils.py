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
    for l2 in feat2['location']:
        if not any(_contains(l1, l2) for l1 in feat1['location']):
            print("No parent was found for location {}".format(l2))
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