using DelimitedFiles

include("centrality.jl")
include("kf.jl")
include("kuramoto_pert.jl")
include("res_dist.jl")

Asp = readdlm("uk_adj_mat.csv",',') .+ 1.0
n = Int(maximum(Asp))
m = Int(size(Asp)[2]/2)

A = zeros(Int,n,n)
for i in 1:330
	A[Int(Asp[i,1]),Int(Asp[i,2])] = 1.0
end

line_list = Array{Array{Int,1},1}()
for i in 1:m
	push!(line_list,Array{Int,1}(vec(Asp[2*i-1,[1,2]])))
end

L = diagm(0 => vec(sum(A, dims = 2))) - A

C = centralities(L)
Om = res_dist(L)

dKf1 = Array{Float64,1}()

for l in line_list
	b = -L[l[1],l[2]]
	Ome = Om[l[1],l[2]]
	push!(dKf1,b*Ome/(1-b*Ome))
end

ranked_lines = line_list[Array{Int,1}(vec(sortslices([dKf1 1:m],dims=1)[:,2]))]

# TODO Compare performance measures when removing lines...


