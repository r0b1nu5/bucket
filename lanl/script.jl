include("generate_time_series.jl")
include("locate_forced.jl")

n = 120
Lsp = Array{Int64,2}(readdlm("data/adj_mat.csv",',').+1)
L = zeros(n,n)
for i in 1:330
	L[Lsp[i,1],Lsp[i,2]] = 1
end
L = diagm(0 => vec(sum(L,dims=1))) - L

m = ones(n)
d = ones(n)

dt = 1e-2

X = generate_time_series(L,m,d,2000,dt)

sol = locate_forced(X[:,101:201],m,d,dt)







