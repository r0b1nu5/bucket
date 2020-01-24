using LinearAlgebra,Statistics


function cnoise(dP::Float64,tau0::Float64)
	f = exp(-1/tau0)
	g = randn()
	
	return f*dP + sqrt(1-f^2)*g
end

