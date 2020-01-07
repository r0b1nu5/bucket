using LinearAlgebra

function ntw_inf_noise(ths::Array{Float64,2}, dP0::Float64)
	n,T = size(ths)
	
	mths = vec(sum(ths,dims=2))/T

	dths = ths - repeat(mths,1,T)

	dP2Ld = dths*Array(dths')

	Lh = pinv(dP2Ld/(dP0^2))

	return Lh
end


function ntw_inf_sine(ths::Array{Float64,2}



