using PyPlot, SparseArrays, Statistics, DelimitedFiles, LinearAlgebra, Distributions, Dates

include("kuramoto.jl")
include("res_dist.jl")
include("rm_line.jl")
include("L2B.jl")
include("isconnected.jl")

thr = 6
h = 1e-4

n_simu = 100

ntw = "ieee118"
Lsp = readdlm(ntw*"_data/"*ntw*"_lap_mat.csv",',')
L = sparse(Lsp[:,1],Lsp[:,2],Lsp[:,3])
Li = pinv(Array(L))

Bsp = readdlm(ntw*"_data/"*ntw*"_inc_mat.csv",',')

n = size(L)[1]
m = round(Int,(size(Lsp)[1]-n)/2)

P0 = .01
ratedP = vec(readdlm(ntw*"_data/P_"*ntw*".csv",','))
P = P0*ratedP
P .-= mean(P)

H = 2*pi .+ 2*rand(n) .- 1.
omega0 = 50*2*pi
rP = 2*max.(abs.(ratedP),maximum(abs.(ratedP))/100*ones(n))
M = 2*H.*rP./omega0

H = 2*pi .+ 2*rand(n) .- 1.
omega0 = 50*2*pi
rP = 2*max.(abs.(ratedP),maximum(abs.(ratedP))/100*ones(n))
D = H.*rP./omega0

cutable_lines = Array{Tuple{Int64,Int64},1}()
for i in 1:n-1
	for j in i+1:n
		if abs(L[i,j]) > 1e-4 && isconnected(rm_line(L,(i,j)))
			if M[i] < M[j]
				push!(cutable_lines,(i,j))
			else
				push!(cutable_lines,(j,i))
			end
		end
	end
end

E = eigen(Array(L))
ls = E.values
us = E.vectors

varP = abs.(P)/3

Ps = zeros(n,n_simu)
for i in 1:n
	Ptemp = rand(Normal(P[i],varP[i]),1,n_simu)
	Ps[i,:] = Ptemp
end
for i in 1:n_simu
	Ps[:,i] .-= mean(Ps[:,i])
end

EP = sum(Ps,dims=2)./n_simu
covP = cov(Ps,dims=2)





