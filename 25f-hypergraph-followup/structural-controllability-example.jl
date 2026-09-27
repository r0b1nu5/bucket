# 1. Define your degree-3 hypergraph
# Each row: [head, tail1, tail2, tail3]
E = [
    6, 1, 4, 5;
    8, 2, 5, 7;
    10, 3, 6, 9;
    # Add more edges as needed
]

n = 10  # Number of state nodes

# 2. Find minimum driver set
Dopt, optSize, allOptimal = brute_force_optimal_drivers(E, n)

println("\n========== RESULTS ==========")
println("Minimum drivers needed: $(optSize)")
println("Driver nodes: $(Dopt)")
println("Total optimal sets found: $(size(allOptimal, 1))")
println("============================\n")
