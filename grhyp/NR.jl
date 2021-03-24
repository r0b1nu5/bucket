using LinearAlgebra

# B: susceptance matrix
# G: conductance matrix
# P: active power vector (for PV and PQ nodes)
# Q: reactive  power verctor (for PQ nodes)
# V: voltage amplitude vector (for PV nodes)
# type: node types (3: slack (node index 1), 0-1: PQ (node indices 2:n2+1), 2: PV (node indices n2+2:n2+n3+1))

function NR(BB::Array{Float64,2}, GG::Array{Float64,2}, P::Array{Float64,1}, Q::Array{Float64,1}, V::Array{Float64,1}, type::Array{Int64,1}, ϵ::Float64=1e-5, max_iter::Int64=15)
	n2 = Int.(sum(type .== 0)) + Int.(sum(type .== 1))
	n3 = Int.(sum(type .== 2))
	n = n2+n3
	N = n+1

	odn = 1 .- diagm(0 => ones(n))
	odN = 1 .- diagm(0 => ones(N))
	odn2 = 1 .- diagm(0 => ones(n2))
	odn3 = 1 .- diagm(0 => ones(n3))

	B = BB #abs.(BB).*(1 .- diagm(0 => ones(N)))
	G = GG #abs.(GG).*(1 .- diagm(0 => ones(N)))

	err = 1000.
	θ = zeros(n)
	v = ones(n2)

	dθ = [0.;θ]*ones(1,N) - ones(N)*[0. θ']
	V2 = [1.;v;V]*[1. v' V']

	S = sin.(dθ)
	C = cos.(dθ)

	GS = G.*S
	GC = G.*C
	BS = B.*S
	BC = B.*C

	p = sum(V2.*(G.*C + B.*S),dims=2)
	q = sum(V2.*(G.*S - B.*C),dims=2)

	err = max(maximum(abs.(p[2:end] - P)),maximum(abs.(q[2:n2+1] - Q)))
	iter = 0

	while err > ϵ && iter < max_iter
		iter += 1
		@info "iter: $iter, error: $err"

		Jpt = odn.*(V2[2:N,2:N].*(GS[2:N,2:N] - BC[2:N,2:N])) - diagm(0 => (q[2:N] + diag(B)[2:N].*diag(V2)[2:N]))
		Jpv = [odn2;ones(n3,n2)].*(repeat([v;V],1,n2).*(GC[2:N,2:n2+1] + BS[2:N,2:n2+1])) + [diagm(0 => (p[2:n2+1]./v + diag(G)[2:n2+1].*v));zeros(n3,n2)]
		Jqt = [odn2 ones(n2,n3)].*(-V2[2:n2+1,2:N].*(GC[2:n2+1,2:N] + BS[2:n2+1,2:N])) + [diagm(0 => (p[2:n2+1] - diag(G)[2:n2+1].*diag(V2)[2:n2+1])) zeros(n2,n3)]
		Jqv = odn2.*(repeat(v,1,n2).*(GS[2:n2+1,2:n2+1] - BC[2:n2+1,2:n2+1])) + diagm(0 => (q[2:n2+1]./v - diag(B)[2:n2+1].*v))

		J = [Jpt Jpv;Jqt Jqv]

		dx = inv(J)*([P;Q] - [p[2:end];q[2:n2+1]])

		θ += dx[1:N-1]
		v += dx[N:N+n2-1]
	
		dθ = [0.;θ]*ones(1,N) - ones(N)*[0. θ']
		V2 = [1.;v;V]*[1. v' V']
		
		S = sin.(dθ)
		C = cos.(dθ)

		GS = G.*S
		GC = G.*C
		BS = B.*S
		BC = B.*C
	
		p = sum(V2.*(GC + BS),dims=2)
		q = sum(V2.*(GS - BC),dims=2)
	
		err = max(maximum(abs.(p[2:end] - P)),maximum(abs.(q[2:n2+1] - Q)))
	end

	return p,q,[0.;θ],[1.;v;V],err,iter
end
		

function GS(B::Array{Float64,2}, G::Array{Float64,2}, P::Array{Float64,1}, Q::Array{Float64,1}, V::Array{Float64,1}, niter::Int64=3)

end


