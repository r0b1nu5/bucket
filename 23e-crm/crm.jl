using DelimitedFiles, PyPlot, LinearAlgebra, Graphs, SimpleWeightedGraphs, Statistics

include("tools.jl")

buses = readdlm("buses.csv",',')
lines = readdlm("lines.csv",',')
gens = readdlm("generators.csv",',')

rm_ids = [243,247,575]

ids = Int64[]
names = String[]
x = Float64[]
y = Float64[]
for i in 2:length(buses[:,1])
	if buses[i,6] == "CH" && !(buses[i,1] in rm_ids)
		push!(ids,Int64(buses[i,1]))
		push!(names,String(buses[i,2]))
		push!(x,Float64(buses[i,8]))
		push!(y,Float64(buses[i,9]))
	end
end

#i2n = Dict{Int64,Int64}(ids[i] => i for i in 1:length(ids))
i2n = Dict{Int64,Int64}()
i2x = Dict{Int64,Float64}()
i2y = Dict{Int64,Float64}()
namr = String[]
X = Float64[]
Y = Float64[]
k = 0
for i in 1:length(ids)
	id = ids[i]
	if !(id in keys(i2n))
		global k += 1
		i2n[id] = k
		i2x[id] = x[i]
		i2y[id] = y[i]
		push!(namr,names[i][1:end-4])
		push!(X,x[i])
		push!(Y,y[i])
		for j in i+1:length(ids)
			if x[i] == x[j] && y[i] == y[j]
				i2n[ids[j]] = k
				i2x[ids[j]] = x[i]
				i2y[ids[j]] = y[i]
			end
		end
	end
end

# Coordinates tweak
X[86] = 8.85

ls = zeros(0,4)
for i in 2:length(lines[:,1])
	if (lines[i,1] in ids) && (lines[i,2] in ids)
		global ls = [ls;lines[[i,],[1,2,3,4]]]
	end
end

n = length(union(collect(values(i2n))))

A0 = zeros(n,n)
A1 = zeros(n,n)
A2 = zeros(n,n)

for k in 1:length(ls[:,1])
	l = ls[k,:]
	i = i2n[l[1]]
	j = i2n[l[2]]
	A0[i,j] = 1.
	A0[j,i] = 1.
	A1[i,j] += norm(l[3:4])
	A1[j,i] += norm(l[3:4])
	A2[i,j] += norm(1 ./l[3:4])
	A2[j,i] += norm(1 ./l[3:4])
end
d0 = vec(sum(A0,dims=1))
D0 = diagm(0 => d0)
L0 = D0 - A0
d1 = vec(sum(A1,dims=1))
D1 = diagm(0 => d1)
L1 = D1 - A1
d2 = vec(sum(A2,dims=1))
D2 = diagm(0 => d2)
L2 = D2 - A2

i2P = Dict{Int64,Float64}(gens[i,1] => Float64(gens[i,3]) for i in 2:size(gens)[1])
P = zeros(length(X))
for id in keys(i2P)
	if id in keys(i2n)
		P[i2n[id]] += min(i2P[id],800.)
	end
end
P .-= mean(P)


 #=
for i in 1:n
	for j in i+1:n
		if A0[i,j] > .1
			PyPlot.plot(X[[i,j]],Y[[i,j]],"-k")
		end
	end
	PyPlot.plot(X[i],Y[i],"ok")
#	PyPlot.text(X[i],Y[i],"$(ids[i])")
	PyPlot.text(X[i],Y[i],namr[i])
end
# =#

g = Graph(A0)

figure("Graph")
plot_ch()
plot_A(A0,X,Y)

 #=
figure("Graph, degree")
v = d0
#vn = normalized(v)
plot_ch()
plot_vscale(A0,X,Y,v,"plasma",cb=true)
str = get_list(v,namr,20)
PyPlot.plot(maximum(X)+.1,minimum(y),"xk")
PyPlot.text(maximum(X)+.1,minimum(Y),str)

figure("Graph, betweenness c.")
v = betweenness_centrality(g)
#vn = normalized(v)
plot_ch()
plot_vscale(A0,X,Y,v,"plasma",cb=true)
str = get_list(v,namr,20)
PyPlot.plot(maximum(X)+.1,minimum(y),"xk")
PyPlot.text(maximum(X)+.1,minimum(Y),str)
# =#

 #=
V = pinv(L0)*P
dV = V*ones(1,n) - ones(n)*V'
F = dV.*A0

figure("Graph, DC flows")
plot_ch()
plot_escale(F,X,Y,"plasma",cb=true,cbl="DC flow")
# =#

# #=
V = pinv(L0)*P
dV = V*ones(1,n) - ones(n)*V'
I = dV.*A0
figure("Graph, DC flows")
plot_ch()
plot_vescale(abs.(I),X,Y,P,"coolwarm","rainbow",cbv=true,cbvl="Power",cbe=true,cbel="DC flow")
# =#




