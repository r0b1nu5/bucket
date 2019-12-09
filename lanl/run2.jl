using Distributed

@everywhere include("final_plot_cluster.jl")

pmap(plts,tups)



