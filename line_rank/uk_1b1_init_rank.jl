using DelimitedFiles,Statistics,LinearAlgebra,PyPlot,SparseArrays

include("rm_line.jl")
include("isconnected.jl")
include("kuramoto.jl")
include("uk_gen_idx.jl")
include("res_dist.jl")

#P0 = 5.0
M0 = .2
D0 = .1

eps = 1e-6
max_iter = 100000
h = .1

Asp = readdlm("uk_adj_mat.csv",',') .+ 1

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
P[gen_idx] = .2*ones(length(gen_idx))
P .-= mean(P)

x1 = vec(readdlm("uk_sync_$P0.csv",','))

#=
th0 = zeros(n)
omeg0 = zeros(n)
x0 = [th0;omeg0]

xs1,dxs1 = kuramoto2(L,M,D,P,x0[1:n],x0[(n+1):(2*n)])
x1 = vec(xs1[:,end])
=#

rank1 = Array{Int64,1}(sortslices([dist 1:m],dims=1)[:,2])
ranks1 = Array{Array{Int64,1},1}([rank1,])
run1 = true
rmvd1 = Array{Int64,1}()
cut1 = Array{Int64,1}()
for i in 1:m
	Lt = copy(L)
	l = ll[i]
	Lt[[l[1],l[2]],[l[1],l[2]]] -= Lt[l[1],l[2]]*[-1 1;1 -1]
	if !isconnected(Lt)
		push!(cut1,i)
	end
end
cuts1 = Array{Array{Int64,1},1}([cut1,])
L1 = copy(L)

count = 0

while run1 && count < m - n + 1
	global count,L1
	count += 1
	@info("round $count")
	
	k = 1
	i = setdiff(rank1,cuts1[end],rmvd1)[k]
	Lt = rm_line(L1,ll[i])
	cut1 = cuts1[end]
	while !isconnected(Lt)
		push!(cut1,i)
		k += 1
		i = setdiff(rank1,cuts1[end],rmvd1)[k]
		Lt = rm_line(L1,ll[i])
	end
	L1 = copy(Lt)
	push!(rmvd1,i)
	push!(cuts1,cut1)
	l = ll[i]
	Om = res_dist(L1)
	
	dist = Array{Float64,1}()
	for i in setdiff(1:m,rmvd1)
		l = ll[i]
		push!(dist,Om[l[1],l[2]])
	end
	rank = Array{Int64,1}(sortslices([dist setdiff(1:m,rmvd1)],dims=1)[:,2])
	push!(ranks1,rank)
	
	xs,dxs,n_iter = kuramoto2(L1,M,D,P,x1[1:n],x1[(n+1):(2*n)])
	if n_iter >= max_iter
		run1 = false
	end
	
end

writedlm("data/ranks_init_rank_$P0.csv",ranks1,',')
writedlm("data/rmvd_init_rank_$P0.csv",rmvd1,',')
writedlm("data/cuts_init_rank_$P0.csv",cuts1,',')



