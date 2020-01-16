using PyPlot, DelimitedFiles

include("kuramoto.jl")
include("ntw_inf.jl")

#run = "noise"
run = "sine"

L = readdlm("uk_lap_mat.csv",',')
n = size(L)[1]
P = zeros(n)

if run == "noise"
	dP0 = 1.
	
	@info "Generate time series..."
	ths = kuramoto_noise(L,P,zeros(n),dP0,true,10000)

	@info "Estimate network"
	Lh = ntw_inf_noise(ths,dP0)

	matshow(L)
	matshow(Lh)
	matshow(L./Lh)

elseif run == "sine"
	a0 = .2
	w0 = .01
	p0 = 0.
	
	ths = Array{Array{Float64,2},1}()

	for i in 1:n
		@info "Probing at $i"
		a = zeros(n)
		a[i] = a0
		w = zeros(n)
		w[i] = w0
		p = zeros(n)
		p[i] = p0

		writedlm("data1/ths_p$i.csv",kuramoto_sine(L,P,zeros(n),a,w,p,true,10000),',')
	end

	Ldh = zeros(n,n)
	for i in 1:n
		ths = readdlm("data1/ths_p$i.csv",',')
		for j in 1:n
			Ldh[i,j] = ntw_inf_sine(ths[j,:],n,a,w,p)
		end
	end

	Lh = pinv(Ldh)

	matshow(L)
	matshow(Lh)
	matshow(L./Lh)
end








# generate time series with noise
#
# analyze time series with noise
#
# generate time series with sine probing
#
# analyze time series with sine probing





