using LinearAlgebra, SparseArrays

function cycle(n::Int64, k::Int64=1)
	A = zeros(n,n)
	for i in 1:k
		A += diagm(i => ones(n-i)) + diagm(-i => ones(n-i)) + diagm(n-i => ones(i)) + diagm(i-n => ones(i))
	end

	return diagm(0 => 2*k*ones(n)) - A
end

function sp_cycle(n::Int64)
	A = spzeros(n,n)
	for i in 1:k
		A += spdiagm(i => ones(n-i)) + spdiagm(-i => ones(n-i)) + spdiagm(n-i => ones(i)) + spdiagm(i-n => ones(i))
	end

	return spdiagm(0 => 2*k*ones(n)) - A
end



