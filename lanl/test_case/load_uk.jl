using DelimitedFiles

include("../system_identification.jl")

L = readdlm("uk_lap_mat.csv",',')
n = size(L)[1]

m = ones(n)
d = ones(n)

dt = .1
T1 = 2000
T2 = 10000

Xs = readdlm("uk_$(T1)_$(dt).csv",',')

Xf = readdlm("ntw10_forced_0.0015_$(T2)_$(dt).csv",',')




