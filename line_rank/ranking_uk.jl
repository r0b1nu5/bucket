using DelimitedFiles

include("centrality.jl")
include("kf.jl")
include("kuramoto_pert.jl")
include("res_dist.jl")

Asp = readdlm("uk_adj_mat.csv",',') .+ 1.0
n = Int(maximum(Asp))
A = zeros(Int,n,n)

for i in 1:330
	A[Int(Asp[i,1]),Int(Asp[i,2])] = 1.0
end

L = diagm(0 => vec(sum(A, dims = 2))) - A

C = centralities(L)






