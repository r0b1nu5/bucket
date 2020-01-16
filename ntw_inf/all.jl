using PyPlot, DelimitedFiles, LinearAlgebra

include("kuramoto.jl")
include("ntw_inf.jl")

ntw = "uk"

a0,w0 = readdlm(ntw*"_probe.csv",',')
@info "a0 = $a0, w0 = $w0"

L = readdlm(ntw*"_lap_mat.csv",',')

n = size(L)[1]

Ldh = zeros(n,n)

for i in 1:n-1
	@info "i = $i"

	a = zeros(n)
	a[i] = a0
	w = zeros(n)
	w[i] = w0

	thij = kuramoto_sine(L,zeros(n),zeros(n),a,w,zeros(n),true,10000,1e-8,.1)
	for j in i+1:n
		Ldh[i,j] = ntw_inf_sine(thij[j,:],n,a0,w0,0.,.1)
		Ldh[j,i] = Ldh[i,j]
	end
end

Ldh = Ldh - diagm(0 => vec(sum(Ldh,dims=2)))

Lh = pinv(Ldh)

for i in 1:n
	for j in 1:n
		PyPlot.plot(L[i,j],Lh[i,j],"o",color="C0")
	end
end



