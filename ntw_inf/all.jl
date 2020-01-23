using PyPlot, DelimitedFiles, LinearAlgebra, SparseArrays

include("kuramoto.jl")
include("ntw_inf.jl")

ntw = "uk_w"
co = "C3"

a0,w0,p0 = readdlm("ntws_data/"*ntw*"_probe.csv",',')
@info "a0 = $a0, w0 = $w0, ψ0 = $p0"

Lsp = readdlm("ntws_data/"*ntw*"_lap_mat_sp.csv",',')
L = sparse(Lsp[:,1],Lsp[:,2],Lsp[:,3])

T = 1000
Ttot = 20000
Ttot = 20100
h = .1
ep = 1e-8

n = size(L)[1]

Ldh = zeros(n,n)

for i in 1:n
	@info "i = $i"

	a = zeros(n)
	a[i] = a0
	w = zeros(n)
	w[i] = w0
	p = zeros(n)
	p[i] = p0

	thij = kuramoto_sine(L,zeros(n),zeros(n),a,w,p,T,Ttot,ep,h)
	for j in i:n
		Ldh[i,j] = ntw_inf_sine(thij[j,:],n,a0,w0,p0,h)
		Ldh[j,i] = Ldh[i,j]
	end
end

Lh = pinv(Ldh)

Id = diagm(0 => ones(n))

figure(33,(5,5))
PyPlot.plot([minimum(L.*(1 .- Id))-.1,maximum(L.*(1 .- Id))+.1],[minimum(L.*(1 .- Id))-.1,maximum(L.*(1 .- Id))+.1],"--k")

for i in 1:n
	for j in i+1:n
		figure(33)
		PyPlot.plot(L[i,j],Lh[i,j],"o",color=co)
	end
end

xlabel("Real value")
ylabel("Inferred value")
title("Network: "*ntw)


