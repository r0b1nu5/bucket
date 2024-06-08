using PyPlot

include("lv.jl")
include("tools.jl")

zer0 = 1e-10

S = 157

#N0 = N[:,end]
# #=
N0 = .1*rand(S)
A = randn(S,S)
# =#

#A = readdlm("data/ex-limit-cycle-01-A.csv",',')

κ = 1.
μ = 5.
σ = 2.7

n_iter = 150_000
N = lv_bunin(N0,A,κ,μ/S,σ/sqrt(S),n_iter,n_iter,1e-3,zer0)

survivors = vec(1:S)[N[:,end] .> zer0]
Ss = length(survivors)

figure()
for s in survivors
	PyPlot.plot(N[s,:])
end
title("# survivors: $(Ss)")





