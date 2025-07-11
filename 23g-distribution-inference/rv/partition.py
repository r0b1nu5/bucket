def get_partition(list_meter: list[str], list_pe: list[list[str]], list_fs: list[list[str]]) -> list[list[str]]:
    """
    Constructs a partition of `list_meter` following the rules in Park et al. (2023), i.e., any two nodes in a given cluster are either a parent-child pair, where the child is leaf node, or  a siblings pair, where both nodes are leaves. 

    _INPUT_:
        list_meter (list[str]): List all meters' names.
        list_pe (list[list[str]]): List of all identified parent-child pairs (parent first).
        list_fs (list[list[str]]): List of all identified siblings pairs (no particular order).

    _OUTPUT_:
        partition (list[list[str]]): List of all created partition. In each cluster, the first node in the list is the parent that will be kept for the next level iteration. This parent is either the meter that was identified in a parent-child relation, or a newly created parent when a siblings relationship was discovered with a parent.
    """

    # Initialize the variables
    partition = []                                      # List of clusters we are creating.
    remaining = list_meter.copy()                       # List of nodes that have not yet been attributed to a cluster.
    multi_p = False
    
    # Dealing with parent-child relationships ###############################################################
    list_e, list_p = get_pe(list_pe)                    # Extracting the list of all children and the list of all parents (which can overlap).
    list_ef = list(set(list_e) - set(list_p))           # Extract the list of children that are not parents.

    for e in list_ef:                                   # Loop over all children that are not parents.
        p = [pe[0] for pe in list_pe if pe[1] == e]         # Find the list of parents for e.
        if len(p) > 1:
            multi_p = True                                  # If multiple parents detected, then raise a flag.
        if p[0] not in remaining:                           
            for part in partition:
                if p[0] in part:
                    part.append(e)                          # If the parent of e has already been treated, then add e to its cluster...
                    remaining.remove(e)                     # ... and remove it from remaining.
        else:
            partition.append([p[0], e])                     # Else create a new partition for e and its parent...
            for item in [p[0], e]:
                remaining.remove(item)                      # ... and remove them from remaining.

    # Dealing with siblings relationships ###############################################################
    for fs in list_fs:                                  # Loop over all siblings relationships.
        if fs[0] in list_p or fs[1] in list_p:              # If either of the siblings is a parent, then do nothing.
            continue
        if fs[0] in list_e and fs[1] in list_e:             # If both are children, then do nothing, because they already belong to a cluster.
            continue
        if fs[0] not in remaining:
            for part in partition:
                if fs[0] in part:
                    part.append(fs[1])                      # If fs[0] has been treated but not fs[1], then add fs[1] to the same cluster as fs[0].
                    if fs[1] in remaining:
                        remaining.remove(fs[1])
        elif fs[1] not in remaining:
            for part in partition:
                if fs[1] in part:
                    part.append(fs[0])                      # If fs[1] has been treated but not fs[0], then add fs[0] to the same cluster as fs[1].
                    if fs[0] in remaining:
                        remaining.remove(fs[0])
        else:
            partition.append(["h", fs[0], fs[1]])           # If neither fs[0] nor fs[1] have been treated, the we need to add their parent "h" and they all belong to the same cluster.
            for item in fs:
                if item in remaining:
                    remaining.remove(item)

    for r in remaining:
        partition.append([r])                           # All the remaining nodes belong to their own cluster.

    if multi_p:
        print("WARNING: Some children had multiple parents!")

    return partition




def get_pe(list_pe: list[list[str]]) -> tuple[list[str], list[str]]:
    list_e = []
    list_p = []
    for pe in list_pe:
        list_p.append(pe[0])
        list_e.append(pe[1])
    return list_e, list_p

