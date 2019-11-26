using Distributed, DelimitedFiles, Distributions

include("scripts.jl")
include("big_rand.jl")

@everywhere include("cluster_fcts.jl")

nx = 3
n0 = 23
#n0 = 4
emi = .05
ema = .6
#ema = .4
ne = 30
epss = Array(LinRange(emi,ema,ne))
n_run = 20
d = 0.
#d = 0.
sig = .2
n1 = 1000
n2 = 1001
n_modes = [2,3,4]

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
				    







