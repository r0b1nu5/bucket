using DelimitedFiles,Statistics

include("kuramoto.jl")
include("uk_gen_idx.jl")

P0 = 1.4
@info P0

Asp = readdlm("uk_adj_mat.csv",',') 

n = Int(maximum(Asp))
m = Int(size(Asp)[1]/2)

A = sparse(vec(Asp[:,1]),vec(Asp[:,2]),ones(size(Asp)[1]))
L = spdiagm(0 => vec(sum(A,dims=2))) - A

mm = 0.2*ones(n)
dd = 0.1*ones(n)

P = zeros(n)
P[gen_idx] = P0*ones(length(gen_idx))
P .-= mean(P)

x0 = zeros(2*n)

x1,dx1 = kuramoto2(L,mm,dd,P,x0[1:n],x0[(n+1):(2*n)])

writedlm("uk_sync_$P0.csv",x1,',')



