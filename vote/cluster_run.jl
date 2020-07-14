using Distributed

n_thr = 9

if nworkers() < n_thr
	addprocs(n_thr - nworkers())
end

@everywhere include("run_country.jl")

t0 = time()

ids = (11,30)
epss = Array(LinRange(.2,.8,20))
Ns = [rand(200:1000,9) for i in 1:3]
Ms = [2*round.(Int64,Ns[i]/100) .+ 1 for i in 1:3]
MUs = [m*ones(9) for m in [.3,.5,.7]]
DEs = [m*ones(9) for m in [0.,.1,.2]]
biass = [m*ones(9) for m in [0.,.05,.1]]
SIs = [m*ones(9) for m in [.1,.2,.3]]

run_country(ids,epss,Ns,Ms,MUs,DEs,biass,SIs)

t1 = time()

@info "Everything done in $((t1-t0)/(60*60)) hours."



