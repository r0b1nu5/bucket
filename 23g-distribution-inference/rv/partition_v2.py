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
    
    list_e, list_p = get_pe(list_pe)                    # Extracting the list of all children and the list of all parents (which can overlap).
    list_ef = list(set(list_e) - set(list_p))           # Extract the list of children that are not parents.
    list_f = list(set([x[0] for x in list_fs]) | set([x[1] for x in list_fs]))   # Extract the list of siblings.
    list_enf = list(set(list_ef) - set(list_f))         # Extract the list of children that are not in a siblings relationship.

    # Dealing with siblings relationships ###############################################################
    # The goal is first to create groups of siblings possibly more that groups of two).
    for fs in list_fs:                                  # Loop over all siblings pairs.
        if fs[0] in list_p or fs[1] in list_p:          # If either of the siblings is a parent, then do nothing.
            continue
        
        if fs[0] in remaining and fs[1] in remaining: 
            partition.append(fs)                        # If none of the siblings have been treated before, then create a new cluster for them.
            for item in fs:
                remaining.remove(item)                  
        elif fs[0] in remaining:                        # If one of the siblings already belongs to a group, add the other one to the same group.
            for part in partition:
                if fs[1] in part:
                    partt = list(set(part) | set([fs[0]]))
                    partition.remove(part)
                    partition.append(partt)
                    remaining.remove(fs[0])
        elif fs[1] in remaining:                        # Idem
            for part in partition:
                if fs[0] in part:
                    partt = list(set(part) | set([fs[1]]))
                    partition.remove(part)
                    partition.append(partt)
                    remaining.remove(fs[1])
        else:                                           # If both of the siblings have been treated before, then merge their clusters.
            part0 = []
            part1 = []
            for part in partition:
                if fs[0] in part:
                    part0 = part
                if fs[1] in part:
                    part1 = part
            if part0 == part1:
                continue
            else:
                part01 = list(set(part0) | set(part1))
                partition.append(part01)
                partition.remove(part0)
                partition.remove(part1)

##### Dealing with parents ##################################
# The goal is to attribute a unique parent to each siblings cluster.
    partition2 = []
    for part in partition:          # For each siblings cluster...
        part = list(set(part))          # Remove repetitions
        ps = [[pe[0] for pe in list_pe if pe[1] == x] for x in part]
        pi = list(set(ps[0]).intersection(*ps[1:]))     # Extract the common parents
        pu = list(set().union(*ps))                     # Extract the union of parents
        partt = []
        if len(pi) == 0 and len(pu) == 0:               # If no parents at all, create a new one.
            partt = ["h"] + part                       
        elif len(pi) == 0:                              # If no common parents, take (arbitrarily) the first one of the union as a common parent.
            if pu[0] in remaining:
                partt = [pu[0]] + part
                remaining.remove(pu[0])
            else:
                for part2 in partition:
                    if pu[0] in part2:
                        partt = part2 + part
                        partition.remove(part2)
        else:                        # If there is something in the intersection of parents, then take the first one as the parent of the group.
            if pi[0] in remaining:
                partt = [pi[0]] + part
                remaining.remove(pi[0])
            else:
                for part2 in partition:
                    if pi[0] in part2:
                        partt = part2 + part
                        partition.remove(part)

        partition2.append(partt)

    for e in list_enf:                                   # Loop over all children that are not siblings.
        p = [pe[0] for pe in list_pe if pe[1] == e]         # Find the list of parents for e.
        if len(p) > 1:
            multi_p = True                                  # If multiple parents detected, then raise a flag.
        if p[0] not in remaining:                           
            for part in partition2:
                if p[0] in part:
                    part.append(e)                          # If the parent of e has already been treated, then add e to its cluster...
                    remaining.remove(e)                     # ... and remove it from remaining.
        else:
            partition2.append([p[0], e])                     # Else create a new partition for e and its parent...
            for item in [p[0], e]:
                remaining.remove(item)                      # ... and remove them from remaining.


    # Dealing with remaining nodes ###############################################################
    for r in remaining:
        partition2.append([r])                           # All the remaining nodes belong to their own cluster.

    if multi_p:
        print("WARNING: Some children had multiple parents!")

    return partition2




def get_pe(list_pe: list[list[str]]) -> tuple[list[str], list[str]]:
    list_e = []
    list_p = []
    for pe in list_pe:
        list_p.append(pe[0])
        list_e.append(pe[1])
    return list_e, list_p

