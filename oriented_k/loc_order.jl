using LinearAlgebra

function loc_r(th::Array{Float64,1}, A::Array{Float64,2})
	n = length(th)

	d = vec(sum(A,dims=2))
	z = exp.(im*A*x)./d
	
	r = norm.(z)
	phi = angle.(z)

	return z,r,phi
end


			


