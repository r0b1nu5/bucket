using LinearAlgebra,SparseArrays

function isconnected(L::Union{Array{Float64,2},SparseMatrixCSC{Float64,Int64}},tolerance::Float64=1e-8)
	ls = eigvals(Array(L))
	ncomp = sum(abs.(ls) .< tolerance)
	connected = ncomp < 2
	return connected
end


