using LinearAlgebra

function get_jac(θ::Vector{Float64},
		 A::Matrix{Float64})
	n = length(θ)
	m = size(A)[1]

	b = [2 -1 -1;
	     -1 2 -1;
	     -1 -1 2]

	J = zeros(n,n)
	for l in 1:m
		ijk = Int64.(A[l,1:3])
		a = A[l,4]
		
		J[ijk,ijk] -= a*b*diagm(0 => cos.(b*θ[ijk]))
	end

	return J
end



