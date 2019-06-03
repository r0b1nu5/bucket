using DelimitedFiles,Statistics,PyPlot,SparseArrays,LinearAlgebra

include("kuramoto.jl")

ntw = "ieee118"
Lsp = readdlm(ntw*"_data/"*ntw*"_lap_mat.csv",',')
L = sparse(Lsp[:,1],Lsp[:,2],Lsp[:,3])

Om = res_dist(L)

n = size(L)[1]
m = round(Int,size(Lsp)[1]/2)

P0 = .001
ratedP = vec(readdlm(ntw*"_data/P_"*ntw*".csv",','))
P = P0*ratedP
P .-= mean(P)

g = .5
# For inertias, see Laurent's Plos One, Appendix 2, Eq. (S2)
H = 2*pi .+ 2*rand(n) .- 1.
omega0 = 50*2*pi
rP = 2*max.(abs.(ratedP),maximum(abs.(ratedP))/100*ones(n))
M = 2*H.*rP./omega0
D = g*M
LD = diagm(0 => D)^(-1/2)*L*diagm(0 => D)^(-1/2)

# Lines to cut: largest bij/mi, median bij/mi, smallest bij/mi,...
# Distribution P: large variance on slow modes, large variance on fast modes,...







