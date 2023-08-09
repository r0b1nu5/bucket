using PyPlot, DelimitedFiles

include("hyper_inf.jl")

t = .01

Xs = Matrix(readdlm("data/coupled_lorenz_solution.txt")')
Ys = (Xs[:,2:end] - Xs[:,1:end-1])./t

n,T = size(Ys)

ooi = [2,3]

xxx = hyper_inf(Xs[:,1:end-1],Ys,ooi,3,-.1)

A2,AA2 = inferred_adj_2nd(xxx[1][2],n)
A3,AA3 = inferred_adj_3rd(xxx[1][3],n)

d = 3




