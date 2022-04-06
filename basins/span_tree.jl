using LinearAlgebra

function span_tree(eWeights::Vector{Float64}, spanTrees::Vector{Vector{Int64}})
	C = 0.

	for el in spanTrees
		C += prod(eWeights[el])
	end

	return C
end

function cycle_span_trees(n::Int64)
	st = Vector{Vector{Int64}}()

	for i in 1:n
		push!(st,Vector([1:i-1;i+1:n]))
	end

	return st
end



