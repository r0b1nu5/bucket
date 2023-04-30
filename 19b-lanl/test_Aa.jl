using PyPlot, DelimitedFiles, LinearAlgebra

include("final.jl")
include("final_Aa.jl")

ex = 1
Xs = readdlm("data_melvyn/ieee57/ieee57_ex$(ex)_Xs.csv",',')

ntw = "ieee57"
ks = 1:160
ls = 1:57
τ = .1
fs = 17
ff = 1.85/2π

nn,N = size(Xs)
n = Int64(nn/2)
T = (N-1)*τ

adj = readdlm("data_melvyn/ieee57/ieee57_adj.csv",',')
A = zeros(n,n)
for k in 1:size(adj)[1]
	i = Int64(adj[k,1])
	j = Int64(adj[k,2])
	A[i,j] = adj[k,3]
end
D = diagm(0 => vec(sum(A,dims=1)))
L = A - D
a = ones(n) 

#L0 = run_l0(Xs,τ,Vector(ls),Vector(ks),false,true)
#L0_Aa = run_l0_Aa(Xs,τ,L,a,Vector(ls),Vector(ks),true)
#L1 = run_l1(Xs,τ,Vector(ks),false,true)
L1_Aa = run_l1_Aa(Xs,τ,L,a,Vector(ks),true)


