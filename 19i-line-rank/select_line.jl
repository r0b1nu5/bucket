using DelimitedFiles,SparseArrays

include("rm_line.jl")
include("isconnected.jl")

ntw = "ieee118"

Lsp = readdlm(ntw*"_data/"*ntw*"_lap_mat.csv",',')
L = sparse(Lsp[:,1],Lsp[:,2],Lsp[:,3])

n = size(L)[1]
m = round(Int,size(Lsp)[1]/2)

g = .5
H = 2*pi .+ 2*rand(n) .- 1
omega0 = 50*2*pi
ratedP = vec(readdlm(ntw*"_data/P_"*ntw*".csv",','))
rP = 2*max.(abs.(ratedP),maximum(abs.(ratedP))/100*ones(n))
M = 2*H.*rP./omega0

crap = Array{Float64,2}(undef,0,4)

for i in 1:n
	for j in 1:n
		global L,crap,M
		if i !== j && abs(L[i,j]) > 1e-4
			crap = [crap;abs(L[i,j])/M[i] i j isconnected(rm_line(L,(i,j)))]
		end
	end
end





