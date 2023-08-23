using PyPlot, DelimitedFiles

include("hyper_inf.jl")

t = .01
d = 3

Xs = Matrix(readdlm("data/coupled_lorenz_solution.txt")')
Ys = (Xs[:,2:end] - Xs[:,1:end-1])./t

N,T = size(Ys)
n = Int64(N/d)

ooi = [2,3]

xxx = hyper_inf(Xs[:,1:end-1],Ys,ooi,3,d,-.1,)


#A2,AA2 = inferred_adj_2nd(xxx[1][2],n)
#A3,AA3 = inferred_adj_3rd(xxx[1][3],n)

x = readdlm("lorenz-edges.txt",'}')
E2 = Vector{Int64}[]
for s in x[1,:]
	global e = Int64[]
	for c in s
		if !(c in [',','[',']','{',' '])
			push!(e,parse(Int64,c)+1)
		end
	end
	if length(e) > 0
		push!(E2,e)
	end
end
E3 = Vector{Int64}[]
for s in x[2,:]
	global e = Int64[]
	for c in s
		if !(c in [',','[',']','{',' '])
			push!(e,parse(Int64,c)+1)
		end
	end
	if length(e) > 0
		push!(E3,e)
	end
end

A2 = zeros(n,n)
for ij in E2
	i,j = ij
	A2[i,j] = 1.
	A2[j,i] = 1.
end
A3 = zeros(n,n,n)
for ijk in E3
	i,j,k = ijk
	A3[i,j,k] = 1.
	A3[j,k,i] = 1.
	A3[k,i,j] = 1.
	A3[i,k,j] = 1.
	A3[j,i,k] = 1.
	A3[k,j,i] = 1.
end






