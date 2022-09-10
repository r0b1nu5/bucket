using DelimitedFiles, Graphs, LinearAlgebra

include("tools.jl")

L = readdlm("ntw_data/rts96_L.csv",',')
A = diagm(0 => diag(L)) - L

g = SimpleGraph(A)
cb = cycle_basis(g)

str = "Σ = ["
for c in cb
	global str = str*"["
	for i in 1:length(c)
		str = str*"$(c[i]),"
	end
	str = str*"],"
end
str = str*"]"

B,w = L2B(L)
n,m = size(B)

C = zeros(length(cb),m)

for j in 1:length(cb)
	cc = [cb[j];cb[j][1]]
	σ = zeros(m)

	for i in 1:(length(cc)-1)
		x,e = findmin(B[cc[i],:].*B[cc[i+1],:])
		if B[cc[i],e] > 0.
			σ[e] = 1.
		else
			σ[e] = -1.
		end
	end
	C[j,:] = σ
end

a,b = size(C)
str2 = "C = ["
for i in 1:a-1
	global str2
	for j in 1:b
		str2 = str2*"$(C[i,j]) "
	end
	str2 = str2*";"
end
for j in 1:b
	global str2 = str2*"$(C[a,j]) "
end
str2 = str2*"]"


