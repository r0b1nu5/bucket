using LinearAlgebra 

include("kuramoto.jl")
include("splay_states.jl")
include("cycle.jl")


# Explores the basin of the q-twisted state with a random walk.
# Returns the sequence of steps of the random walk.
function rwalk(q::Int64, n::Int64, δ::Float64, T::Int64, ρ::Float64=.5)
	θq = splay_q(q,n)
	
	L = cycle(n)
	ω = zeros(n)

	θs = Matrix{Float64}(undef,n,0)
	θqs = Matrix{Float64}(undef,n,0)

	θ = kuramoto(L,ω,θq)
	θ1 = θ
	θ2 = θ

	if winding(θ,Array(1:n)) != q
		@info "No consistent initial conditions..."
		return θs
	else
		θs = [θs θ]
		θqs = [θqs θ]
	end

	dq = 100

	while dq != 0
		v1 = 2*rand(n) .- 1
		v1 ./= norm(v1)
		v3 = v1 .- mean(v1)
		v3 ./= norm(v3)
		θ1 = θ + δ*v3
		θ2 = kuramoto(L,ω,θ1)
		dq = q - winding(θ2,Array(1:n))
	end

	θ = θ1
	θs = [θs θ]
	θqs = [θqs θ2]

	t = 1
	c = 0

	while t < T
		t += 1
		if t%1000 == 0
			c += 1
			@info("iter = $t, d = $(dist(θ,θq))")
			writedlm("data/th_$c.csv",θs[:,1:end-1],',')
			θs = θs[:,end]
			writedlm("data/thq_$c.csv",θqs[:,1:end-1],',')
			θqs = θqs[:,end]
		end

		dq = 100

		while dq != 0
			v1 = 2*rand(n) .- 1
			v1 ./= norm(v1)
			v2 = (θ - θ2)./norm(θ - θ2)
			v3 = (v1 + ρ*v2) .- mean(v1 + ρ*v2)
			v3 ./= norm(v3)

			θ1 = θ + δ*v3
			θ2 = kuramoto(L,ω,θ1)
			dq = q - winding(θ2,Array(1:n))
		end

		θ = θ1
		θs = [θs θ]
		θqs = [θqs θ2]
	end

	Θs = Matrix{Float64}(undef,n,0)
	Θqs = Matrix{Float64}(undef,n,0)
	for i in 1:c
		Θs = [Θs readdlm("data/th_$i.csv",',')]
		rm("data/th_$i.csv")
		Θqs = [Θqs readdlm("data/thq_$i.csv",',')]
		rm("data/thq_$i.csv")
	end
	Θs = [Θs θs]
	Θqs = [Θqs θqs]

	return Θs,Θqs
end




