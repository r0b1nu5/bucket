using LinearAlgebra

function L2B(L::Matrix{Float64})
	n = size(L)[1]

	Bout = zeros(n,0)
	Bin = zeros(n,0)
	w = Float64[]

	for i in 1:n
		for j in 1:n
			if i !== j && abs(L[i,j]) > 1e-6
				Bout = [Bout I[1:n,i]]
				Bin = [Bin I[1:n,j]]
				push!(w,-L[i,j])
			end
		end
	end

	return Bout,Bin,w
end

function winding(θ::Vector{Float64}, σ::Vector{Int64})
	Δ = θ[[σ[2:end];σ[1]]] - θ[σ]
	return round(Int64,sum(mod.(Δ .+ π,2π) .- π)/(2π))
end

function gen_cycle_undir(n::Int64)
	return diagm(0 => 2*ones(n)) - diagm(1 => ones(n-1)) - diagm(-1 => ones(n-1)) - diagm(n-1 => ones(1)) - diagm(1-n => ones(1))
end

function gen_cycle_dir(n::Int64)
	return diagm(0 => ones(n)) - diagm(1 => ones(n-1)) - diagm(1-n => ones(1))
end



