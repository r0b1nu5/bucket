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
#	t = round(Int64,(-2*p0)/(2*w0*h))
#	th_ij = th_ij .- mean(th_ij[:,t+1]) .+ a0*cos(w0*h*(t+1) + p0)/(n*w0)
#	th_ij = th_ij .- th_ij[t+1] .- a0/(n*w0)*(1 - cos(w0*h*(t+1) + p0))

#	num = (th_ij + a0*cos.(w0*h*Array(0:(length(th_ij)-1)) .+ p0)./(n*w0))
	num = (th_ij - a0/(n*w0)*(1 .- cos.(w0*h*Array(0:(length(th_ij)-1)) .+ p0)))
	denom = a0*sin.(w0*h*Array(0:(length(th_ij)-1)) .+ p0)
	ids = setdiff((1:length(denom)).*(abs.(denom) .> 1e-4),[0,])

	Lijds = num[ids]./denom[ids]
	
	Lijdh = median(Lijds)

	return Lijdh
end



