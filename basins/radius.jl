using LinearAlgebra

include("splay_states.jl")
include("kuramoto.jl")
include("cycle.jl")

# Explores the basin of the q-twisted state with a straight line in the direction of v.
# Return magnitude of the perturbation that leaves the basin. 
# !!! v is not normalized in the script !!!
function radius(q::Int64, n::Int64, v::Array{Float64,1}, δ::Float64=.01, αmax::Float64=2π)
	θq = splay_q(q,n)

	L = cycle(n)
	ω = zeros(n)

	θ = kuramoto(L,ω,θq)
	α = 0.

	if winding(θ,Array(1:n)) != q
		@info "No consistent initial conditions..."
		return α
	end

	dq = 0

	while dq == 0 && α < αmax
		α += δ
		θ = kuramoto(L,ω,θq + α*v)
		dq = q - winding(θ,Array(1:n))
	end

	return α
end



