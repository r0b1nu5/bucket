using LinearAlgebra

# B: susceptance matrix
# G: conductance matrix
# P: active power vector (for PV and PQ nodes)
# Q: reactive  power verctor (for PQ nodes)
# V: voltage amplitude vector (for PV nodes)
# type: node types (1: slack (node index 1), 2: PQ (node indices 2:n2+1), 3: PV (node indices n2+2:n2+n3+1))

function NR(B::Array{Float64,2}, G::Array{Float64,2}, P::Array{Float64,1}, Q::Array{Float64,1}, V::Array{Float64,1}, type::Array{Int64,1})
	n2 = Int.(sum(type .== 2))
	n3 = Int.(sum(typ2 .== 3))
	n = n2+n3
	N = n+1

	err = 1000.
	t = zeros(n)
	v = ones(n2)

	dt = [0.;t]*ones(1,N) - ones(N)*[0. t']
	V2 = [1.;v;V]*[1. v' V']

	p = sum(V2.*(G.*cos.(dt) + B.*sin.(dt)),dims=2)
	q = sum(V2.*(G.*sin.(dt) - B.*cos.(dt)),dims=2)

	err = max(maximum(abs.(p[2:end] - P)),maximum(abs.(q[2:n2+1] - Q)))
	iter = 0

	while err > eps && iter < max_iter
		iter += 1

		Jpt = #TODO
		Jpv
		Jqt
		Jqv

		J

		dx = inv(J)*([P;Q] - [p;q])

		t += dx[...]
		v += dx[...]



