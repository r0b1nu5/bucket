using LinearAlgebra, Statistics

function ntw_inf_noise(ths::Array{Float64,2}, dP0::Float64)
	n,T = size(ths)
	
	mths = vec(sum(ths,dims=2))/T

	dths = ths - repeat(mths,1,T)

	dP2Ld = dths*Array(dths')

	Lh = pinv(dP2Ld/(dP0^2))

	return Lh
end

# th_ij: angle time series at node j when node i is submitted to a sine perturbation.

function ntw_inf_sine(th_ij::Array{Float64,1}, n::Int64, a0::Float64, w0::Float64, p0::Float64, h::Float64=.1)
	Lijds = (th_ij - a0*sin.(w0*h*Array(0:(length(th_ij)-1)) .+ p0)./(n*w0))./(a0*cos(w0*h*Array(0:(length(th_ij)-1)) .+ p0))

	Lijdh = mean(Lijds)

	return Lijdh
end



