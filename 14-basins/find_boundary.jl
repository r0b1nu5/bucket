using PyPlot, LinearAlgebra

include("cycle.jl")
include("kuramoto.jl")
include("splay_states.jl")

n = 6
tol = 1e-15

L = cycle(n)

q0 = 1
θi = splay_q(q0,n)

P = 2*rand(n) .- 1

q = q0 

δ = 1.

while δ > tol
	global θi,q0,n,L,δ

	q = q0
	
	while q == q0
		θi += δ*P
		θ,iter = kuramoto(L,zeros(n),θi,false)
		q = winding(θ,Vector(1:n))
	end

	θi -= δ*P
	δ /= 10
end



