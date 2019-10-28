using Distributed, DelimitedFiles, Distributions

include("scripts.jl")
include("big_rand.jl")

@everywhere include("cluster_fcts.jl")

nx = 4
n0 = 0
#n0 = 4
emi = 0.
ema = 1.
#ema = .4
ne = 40
epss = Array(LinRange(emi,ema,ne))
n_run = 100
d = 1.
#d = 0.
sig = .2
n1 = 1000
n2 = 1001

xess = Array{Tuple{Array{Float64,1},Float64,String,Int64},1}()

for i in n0+1:n0+nx
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
				    







