using PyPlot, FFTW

include("tools.jl")

function N_inf(thi, thj, a0, w0, N, plt::Bool=false)
	T = length(thi)

	ti = thi .- thi[1]
	tj = thj .- thj[1]
	tij = ti - tj

	Ts = setdiff((1:T-1).*(tij[1:end-1].*tij[2:end] .< 0.),[0.,])

	dhs = .25*sum([ti[Ts] ti[Ts.+1] tj[Ts] tj[Ts.+1]],dims=2)

	nhs = (2*a0)./(w0*dhs)
	if abs(nhs[1]-N) < nhs[2]-N
		Nhs = nhs[1:2:length(nhs)]
	else
		Nhs = nhs[2:2:length(nhs)]
	end

	Nhi = (2*a0)/(w0*maximum(ti))
	Nhj = (2*a0)/(w0*maximum(tj))

	if plt
		PyPlot.plot(N*ones(length(Nhs)),Nhs,"o",N,Nhi,"o",N,Nhj,"o")
	end

	return Nhs,Nhi,Nhj
end

function compare_fft(thi::Array{Float64,1}, thj::Array{Float64,1}, a0::Float64, w0::Float64, N::Int64=0, l2::Float64=0.)
	Nhs,Nhi,Nhj = N_inf(thi,thj,a0,w0,N)
	
	figure()
	subplot(1,2,2)
	PyPlot.plot(N*ones(length(Nhs)),Nhs,"o",N,Nhi,"o",N,Nhj,"o")

	ti = thi .- thi[1]
	tj = thj .- thj[1]
	T = length(ti)

	subplot(2,2,1)
	PyPlot.plot(.01*(1:T),ti,.01*(1:T),tj)

	fi = norm.(fft(ti))
	fj = norm.(fft(tj))
	
	T2 = round(Int,T/2)
	T2 = 100
	subplot(2,2,3)
	PyPlot.plot(2pi/(.01*T)*(0:T2),fi[1:T2+1],2pi/(.01*T)*(0:T2),fj[1:T2+1])
end















