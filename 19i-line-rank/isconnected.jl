using LinearAlgebra,SparseArrays

function isconnected(L::Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}},tolerance::Float64=1e-8)
	ls = eigvals(Symmetric(Array(L)))
	ncomp = sum(abs.(ls) .< tolerance)
	connected = ncomp < 2
	return connected
end


