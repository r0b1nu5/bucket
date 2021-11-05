using LinearAlgebra, SparseArrays

function cycle(n::Int64)
	return L = diagm(0 => 2*ones(n)) - diagm(1 => ones(n-1)) - diagm(-1 => ones(n-1)) - diagm(n-1 => ones(1)) - diagm(1-n => ones(1))
end

function sp_cycle(n::Int64)
	return L = spdiagm(0 => 2*ones(n)) - spdiagm(1 => ones(n-1)) - spdiagm(-1 => ones(n-1)) - spdiagm(n-1 => ones(1)) - spdiagm(1-n => ones(1))
end



