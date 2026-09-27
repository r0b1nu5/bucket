# Lumo-generated (27.9.26), from Can Chen's code
using Combinatorics  # For combinations()

# =====================================================================
# MAIN FUNCTION: Brute-force optimal driver selection
# =====================================================================
function brute_force_optimal_drivers(E, n::Int)
    """
    Enumerate driver subsets in increasing cardinality.
    Returns first set that satisfies accessibility + no-dilation.
    
    Args:
        E: Edge matrix, each row = [head, tail_1, tail_2, tail_3]
        n: Number of state nodes
    
    Returns:
        Dopt: Optimal driver set (vector)
        optSize: Size of optimal set
        allOptimal: All optimal sets found (matrix)
    """
    Dopt = Int[]
    allOptimal = Matrix{Int}(undef, 0, 0)

    for r in 1:n
        combinations_list = collect(combinations(1:n, r))
        successful = Matrix{Int}(undef, 0, r)

        for D in combinations_list
            # Convert vector to array for functions
            D_arr = collect(D)
            
            # Check BOTH conditions
            accessible = walk_reach(E, D_arr, n)
            dilationFree = is_dilation_free(E, D_arr, n)

            if all(accessible) && dilationFree
                successful = vcat(successful, reshape(D_arr, 1, :))
            end
        end

        if !isempty(successful)
            optSize = r
            allOptimal = successful
            Dopt = successful[1, :]
            return Dopt, optSize, allOptimal
        end
    end

    error("No feasible driver set found.")
end

# =====================================================================
# HYPERGRAPH ACCESSIBILITY CHECK
# =====================================================================
function walk_reach(E, drivers, n::Int)
    """
    Hypergraph accessibility.
    A hyperedge fires when ALL distinct tail nodes are accessible.
    
    Args:
        E: Edge matrix
        drivers: Initial driver nodes (accessible by definition)
        n: Number of state nodes
    
    Returns:
        accessible: Boolean array of which nodes are reachable
    """
    accessible = falses(n)
    accessible[drivers] .= true

    changed = true
    while changed
        changed = false
        
        for e in 1:size(E, 1)
            head = E[e, 1]
            tails = unique(E[e, 2:end])
            
            if all(accessible[tails]) && !accessible[head]
                accessible[head] = true
                changed = true
            end
        end
    end
    
    return accessible
end

# =====================================================================
# DILATION-FREE CHECK
# =====================================================================
function is_dilation_free(E, drivers, n::Int)
    """
    Check matching covers every state node.
    Each hyperedge contributes its head to potential coverage.
    
    Args:
        E: Edge matrix (each row = [head, tails...])
        drivers: Driver node indices
        n: Number of state nodes
    
    Returns:
        dilationFree: Boolean flag
    """
    stateHeads = E[:, 1]
    controlHeads = drivers
    allHeads = vcat(stateHeads, controlHeads)

    matchedState = falses(n)
    
    for v in allHeads
        if !matchedState[v]
            matchedState[v] = true
        end
    end

    return all(matchedState)
end

function verify_structural_controllability(E, drivers, n::Int)
    """
    Full verification of accessibility + dilation-free conditions.
    Useful for testing arbitrary driver sets.
    """
    accessible = walk_reach(E, drivers, n)
    dilationFree = is_dilation_free(E, drivers, n)
    
    println("\n=== CONTROLLABILITY TEST ===")
    println("Accessible nodes: $(findall(accessible))")
    println("All nodes accessible: $(all(accessible))")
    println("Dilation-free: $(dilationFree)")
    
    return all(accessible) && dilationFree
end

function count_all_optimal_sets(E, n::Int, max_size::Int=5)
    """
    Count ALL optimal driver sets up to a maximum size.
    Useful when there are multiple equally good solutions.
    """
    for r in 1:max_size
        Dopt, optSize, allOptimal = try
            brute_force_optimal_drivers(E, n)
        catch e
            continue
        end
        
        if optSize == r
            return optSize, allOptimal
        end
    end
    return 0, Matrix{Int}(undef, 0, 0)
end
