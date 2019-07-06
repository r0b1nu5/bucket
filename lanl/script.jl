include("generate_time_series.jl")
include("locate_forced.jl")

#=
n = 120
Lsp = Array{Int64,2}(readdlm("data/adj_mat.csv",',').+1)
L = zeros(n,n)
for i in 1:330
	L[Lsp[i,1],Lsp[i,2]] = 1
end
L = diagm(0 => vec(sum(L,dims=1))) - L
=# 

n = 5
L = 	Array{Float64,2}([3 -1 -1 0 -1;
	-1 2 -1 0 0;
	-1 -1 3 -1 0;
	0 0 -1 2 -1;
	-1 0 0 -1 2])
m = ones(n)
d = ones(n)

err = Array{Float64,1}()

Ts = [100,200,400,800,1600,3200,6400,12800]
dts = [1e-3,2e-3,4e-3,8e-3,1.6e-2,3.2e-2,6.4e-2,1.28e-1,2.56e-1,5.12e-1]

for dt in dts
	global err,L,m,d
	
#	dt = 1e-1
	T = 1500

@info "Generate time series...δt = $dt "	
	X = generate_time_series(L,m,d,T,dt)
	
	sol = locate_forced(X,m,d,dt)

	push!(err,maximum(abs.(L-sol.L)))
end

PyPlot.plot(dts,err,"-o")






