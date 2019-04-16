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
# =#

rank2 = Array{Int64,1}(sortslices([dist 1:m],dims=1)[:,2])
ranks2 = Array{Array{Int64,1},1}([rank2,])
run2 = true
rmvd2 = Array{Int64,1}()
cut2 = Array{Int64,1}()
for i in 1:m
	Lt = copy(L)
	l = ll[i]
	Lt[[l[1],l[2]],[l[1],l[2]]] -= Lt[l[1],l[2]]*[-1 1;1 -1]
	if !isconnected(Lt)
		push!(cut2,i)
	end
end
cuts2 = Array{Array{Int64,1},1}([cut2,])
L2 = copy(L)

count = 0

if x1[1] == "nope"
	global cuts2,ranks2,rmvd2,run2
	@info "$(now()) -- IEEE57: No sync possible!"
	cuts2 = Array{Array{Int64,1},1}()
	ranks2 = Array{Array{Int64,1},1}()
	run2 = false
end

while run2 && count < m - n + 1
	global count,L2,rank2,run2,cuts2,ranks2,rmvd2
	count += 1
	@info("$(now()) -- IEEE57: round $count")
	
	k = 1
	i = setdiff(rank2,cuts2[end],rmvd2)[k]
	Lt = rm_line(L2,ll[i])
	cut2 = cuts2[end]
	while !isconnected(Lt)
		push!(cut2,i)
		k += 1
		i = setdiff(rank2,cuts2[end],rmvd2)[k]
		Lt = rm_line(L2,ll[i])
	end
	L2 = copy(Lt)
	push!(rmvd2,i)
	push!(cuts2,cut2)
	l = ll[i]
	Om = res_dist(L2)
	
	dist = Array{Float64,1}()
	for i in setdiff(1:m,rmvd2)
		l = ll[i]
		push!(dist,Om[l[1],l[2]])
	end
	rank2 = Array{Int64,1}(sortslices([dist setdiff(1:m,rmvd2)],dims=1)[:,2])
	push!(ranks2,rank2)
	
	xs,dxs,n_iter = kuramoto2(L2,M,D,P,x1[1:n],x1[(n+1):(2*n)])
	if n_iter >= max_iter
		run2 = false
		@info("$(now()) -- IEEE57: No sync anymore.")
	end
	
end

writedlm("data/ieee57_ranks_detadpu_rank_$(P0)_rev.csv",ranks2,',')
writedlm("data/ieee57_rmvd_detadpu_rank_$(P0)_rev.csv",rmvd2,',')
writedlm("data/ieee57_cuts_detadpu_rank_$(P0)_rev.csv",cuts2,',')

