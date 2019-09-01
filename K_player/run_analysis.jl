using Distributed

n_thr = 5

if nworkers() < n_thr
	addprocs(n_thr - nworkers())
end


@everywhere include("journal_analysis.jl")
@everywhere include("journals.jl")

pmap(journal_analysis_parallel,tups)





