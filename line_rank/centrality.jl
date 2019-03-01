using LinearAlgebra

function centralities(L::Array{Float64,2})
	ls = eigvals(L)
	us = eigvecs(L)
	
	n = length(ls)
	
	C = zeros(n)
	for i in 2:n
		C += n*vec(us[:,i].^2)./ls[i] .+ 1/ls[i]
	end
	
	return C
end


