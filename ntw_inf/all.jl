using PyPlot, DelimitedFiles, LinearAlgebra

include("kuramoto.jl")
include("ntw_inf.jl")

ntw = "uk"

a0,w0 = readdlm(ntw*"_probe.csv",',')
@info "a0 = $a0, w0 = $w0"

L = readdlm(ntw*"_lap_mat.csv",',')

n = size(L)[1]

Ldh = zeros(n,n)

for i in 1:n
	@info "i = $i"

	a = zeros(n)
	a[i] = a0
	w = zeros(n)
	w[i] = w0

	thij = kuramoto_sine(L,zeros(n),zeros(n),a,w,zeros(n),true,1000,1e-8,.1)
	for j in i:n
		Ldh[i,j] = ntw_inf_sine(thij[j,:],n,a0,w0,0.,.1)
		Ldh[j,i] = Ldh[i,j]
	end
end

Lh = pinv(Ldh)

for i in 1:n
	for j in 1:n
		PyPlot.plot(max(1e-8,abs(L[i,j])),abs(Lh[i,j]),"o",color="C0")
	end
end



