using Distributed

addprocs()


@everywhere include("journal_analysis.jl")
@everywhere include("journals.jl")

pmap(journal_analysis_parallel,tups)





