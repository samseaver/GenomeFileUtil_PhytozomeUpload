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
