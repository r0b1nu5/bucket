include("L2B.jl")

function potential_k(Bt::SparseMatrixCSC{Float64,Int64}, ω::Vector{Float64}, θ::Vector{Float64})
	return sum(1 .- cos.(Bt*θ)) - dot(ω,θ)
end

	
