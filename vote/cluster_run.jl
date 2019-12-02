using Distributed, DelimitedFiles, Distributions

include("scripts.jl")
include("big_rand.jl")

@everywhere include("cluster_fcts.jl")

# Generate the tuples in order to compute effort wrt eps (\in [emi,ema]), and run the computation in parallel.

nx = 3 # Number of natural opinion to consider
n0 = 20 # Index of the realization
emi = .05 # Minimal eps value
ema = .6 # Maximal eps value
ne = 30 # Resolution of eps values
epss = Array(LinRange(emi,ema,ne)) # List of eps values
n_run = 20 # Number of runs for the random strategy
d = 1. # Distance between the modes in the distribution of x0
sig = .2 # Standard deviation in the distribution of x0
n1 = 1000 # Number of agents in the left mode
n2 = 1001 # Number of agents in the right mode
n_modes = [2,3,4] # Eigenmodes to consider in the Fiedler strategy

xess = Array{Tuple{Array{Float64,1},Float64,String,Int64},1}()

for i in n0+1:n0+nx
	x0 = zeros(n1+n2)
	dx0 = [1000.,]
	while maximum(dx0) > emi
		x0 = sort([rand(Normal(-d/2,sig),n1);rand(Normal(d/2,sig),n2)])
		dx0 = x0[2:end] - x0[1:end-1]
	end
	writedlm("data/x$i.csv",x0,',')
	for j in 1:ne
		eps = epss[j]
		for st in ["rand","fiedler","mini"]
			if st == "rand"
				for k in 1:n_run
					push!(xess,(x0,eps,"rand$k",i))
				end
			else
				push!(xess,(x0,eps,st,i))
			end
		end
	end
end

   

pmap(eff,xess)
				    







