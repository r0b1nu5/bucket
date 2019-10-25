using Distributed, DelimitedFiles, Distributions

include("scripts.jl")
include("big_rand.jl")

@everywhere include("cluster_fcts.jl")

nx = 2
emi = 0.
ema = 1.
ne = 4
epss = Array(LinRange(emi,ema,ne))
n_run = 2
d = 1.
sig = .2
n1 = 100
n2 = 101

xess = Array{Tuple{Array{Float64,1},Float64,String,Int64},1}()

for i in 1:nx
	x0 = [rand(Normal(-d/2,sig),n1);rand(Normal(d/2,sig),n2)]
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
				    







