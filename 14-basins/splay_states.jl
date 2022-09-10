using Statistics 

function splay_states(n::Int64)
	qmax = round(Int64,n/4)

	splays = Dict{Int64,Array{Float64,1}}()

	for q in -qmax:qmax
		θ = Array(0:n-1)*2π*q/n
		splays[q] = θ .- mean(θ)
	end

	return splays
end

function splay_q(q::Int64, n::Int64)
	θ = Array(0:n-1)*2π*q/n
	
	return mod.(θ .- mean(θ) .+ π,2π) .- π
end

