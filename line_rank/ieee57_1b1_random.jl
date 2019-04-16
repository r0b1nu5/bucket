using DelimitedFiles,Statistics,LinearAlgebra,PyPlot,SparseArrays,Dates

include("rm_line.jl")
include("isconnected.jl")
include("kuramoto.jl")
include("res_dist.jl")

#P0 = 5.0
M0 = .2
D0 = .1

eps = 1e-6
max_iter = 100000
h = .1

Asp = readdlm("ieee57_adj_mat.csv",',')

n = Int(maximum(Asp))
m = Int(size(Asp)[1]/2)

A = sparse(vec(Asp[:,1]),vec(Asp[:,2]),ones(size(Asp)[1]))
L = spdiagm(0 => vec(sum(A,dims=2))) - A

Om = res_dist(L)

ll = Array{Tuple{Int64,Int64},1}()
dist = Array{Float64,1}()
for i in 1:m
	l = (Int64(Asp[2*i-1,1]),Int64(Asp[2*i-1,2]))
	push!(ll,l)
	push!(dist,Om[l[1],l[2]])
end

M = M0*ones(n)
D = D0*ones(n)

P = zeros(n)
P = P0*vec(readdlm("P_57",','))
P .-= mean(P)

x1 = vec(readdlm("sync_states/ieee57_sync_$P0.csv",','))

#=
th0 = zeros(n)
omeg0 = zeros(n)
x0 = [th0;omeg0]

xs1,dxs1 = kuramoto2(L,M,D,P,x0[1:n],x0[(n+1):(2*n)])
x1 = vec(xs1[:,end])
=#

rank3 = Array{Int64,1}(sortslices([dist 1:m],dims=1,rev=true)[:,2])
ranks3 = Array{Array{Int64,1},1}([rank3,])
run3 = true
rmvd3 = Array{Int64,1}()
cut3 = Array{Int64,1}()
for i in 1:m
	Lt = copy(L)
	l = ll[i]
	Lt[[l[1],l[2]],[l[1],l[2]]] -= Lt[l[1],l[2]]*[-1 1;1 -1]
	if !isconnected(Lt)
		push!(cut3,i)
	end
end
cuts3 = Array{Array{Int64,1},1}([cut3,])
L3 = copy(L)

count = 0

if x1[1] == "nope"
	global cuts3,ranks3,rmvd3,run3
	@info "$(now()) -- IEEE57: No sync possible!"
	cuts3 = Array{Array{Int64,1},1}()
	ranks3 = Array{Array{Int64,1},1}()
	run3 = false
end

while run3 && count < m - n + 1
	global count,L3,run3
	count += 1
	@info("$(now()) -- IEEE57: round $count")
	
	cut3 = cuts3[end]
	k = 1
	i = rand(setdiff((1:m),cut3,rmvd3))
	Lt = rm_line(L3,ll[i])
	while !isconnected(Lt)
		push!(cut3,i)
		k += 1
		i = rand(setdiff((1:m),cut3,rmvd3))
		Lt = rm_line(L3,ll[i])
	end
	L3 = copy(Lt)
	push!(rmvd3,i)
	push!(cuts3,cut3)
	l = ll[i]
	Om = res_dist(L3)
	
	dist = Array{Float64,1}()
	for i in setdiff(1:m,rmvd3)
		l = ll[i]
		push!(dist,Om[l[1],l[2]])
	end
	rank = Array{Int64,1}(sortslices([dist setdiff(1:m,rmvd3)],dims=1,rev=true)[:,2])
	push!(ranks3,rank)
	
	xs,dxs,n_iter = kuramoto2(L3,M,D,P,x1[1:n],x1[(n+1):(2*n)])
	if n_iter >= max_iter
		run3 = false
		@info("$(now()) -- IEEE57: No sync anymore.")
	end
	
end

writedlm("data/ieee57_ranks_random_$(P0)_rev.csv",ranks3,',')
writedlm("data/ieee57_rmvd_random_$(P0)_rev.csv",rmvd3,',')
writedlm("data/ieee57_cuts_random_$(P0)_rev.csv",cuts3,',')
 



