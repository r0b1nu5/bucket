using DelimitedFiles

include("../system_identification.jl")

L = readdlm("ntw10_lap_mat.csv",',')
n = size(L)[1]

m = ones(n)
d = ones(n)

dt = .1
T1 = 10000
T2 = 30000

Xs = readdlm("ntw10_$(T1)_$(dt).csv",',')

Xf = readdlm("ntw10_forced_0.1_$(T2)_$(dt).csv",',')




