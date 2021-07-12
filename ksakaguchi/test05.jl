using PyPlot

include("iterations.jl")

L = readdlm("ntw_data/ntw9_L.csv",',')
b,w = L2B(L)
B = [b -b]
Bout = B.*(B .> 0)
n,m = size(B)
P = dir_cycle_proj(B,ones(m))

Δ = diagm(0 => rand(m))

M1 = P'*P*Δ*P

r = Array{Float64,1}()

for i in 1:1000000
	x = 2*rand(m) .- 1
	push!(r,(x'*M1*x)[1])
end

PyPlot.hist(r,1000)
title("Min: $(minimum(r))")



