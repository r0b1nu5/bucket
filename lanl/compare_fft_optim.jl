using FFTW, PyPlot

include("objective.jl")

function compare(Xs::Array{Float64,2}, dt::Float64, Lm::Array{Float64,2}, dm::Array{Float64,1}, a::Array{Float64,1}, f::Array{Float64,1}, phi::Array{Float64,1})
	nn,T = size(Xs)
	n = Int(nn/2)

	f0 = f[1]
	fs = Array(0.: .05*f0 : 5*f0)

	obs = Array{Float64,1}()
	for ff in fs
		push!(obs,objective1(Xs,Lm,dm,a,[ff,0,0,0,0],phi,dt,0.,0.))
	end

	fX = norm.(fft(Xs[6,:]))

	figure(66)
	subplot(121)
	PyPlot.plot((0:T-1)./(dt*T),fX,label="a0 = $(a[1])")
	legend()

	subplot(122)
	PyPlot.plot(fs,obs,label="a0 = $(a[1])")
	legend()
end



