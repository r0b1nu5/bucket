using Distributed

@everywhere include("final_plot_cluster.jl")

pmap(plt,tups)



