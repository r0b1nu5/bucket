

function mle_pl_cutoff(x::Array{Float64,1})
	side = 2.99
	atemp = 3.
	ltemp = 3.
	
	L = zeros(100,100)
	
	while side > 1e-4
		as = LinRange(max(atemp-side,side/100),atemp+side,100)
		ls = LinRange(max(ltemp-side,side/100),ltemp+side,100)
		
		for i in 1:100
			for j in 1:100
				L[i,j] = -length(x)*log(real(polylog(as[i],Complex(exp(-ls[j]))))) - sum(as[i].*log.(x) + ls[j].*x)
			end
		end
		
		atemp = as[findmax(L)[2][1]]
		ltemp = ls[findmax(L)[2][2]]
		side /= 10
	end
	
	return atemp,ltemp
end


	
	

function polylog(s::Float64,z::Complex{Float64})
	Li = gamma(1-s)/(2pi)^(1-s)*((im)^(1-s)*zeta(1-s,.5+log(-z)/(2pi*im)) + (im)^(s-1)*zeta(1-s,.5-log(-z)/(2pi*im)))
	return Li
end


