using PyPlot

include("kuramoto.jl")
include("cycle.jl")
include("splay_states.jl")
include("jacobian.jl")
include("span_tree.jl")

n = 83
q = 15

res = 20

L = cycle(n)
B,w,Bt = L2B(L)
θ0 = splay_q(q,n)

P1 = [ones(41);zeros(42)]
P2 = [zeros(42);ones(41)]

αs = LinRange(-π,π,res)

STs = cycle_span_trees(n)

Qf = zeros(Int64,res,res)
Λs = Matrix{Vector{Float64}}(undef,res,res)
Cst = zeros(res,res)

for i in 1:res
	for j in 1:res
		@info "(i,j) = ($i,$j)"
		
		global Λs,Qf,αs,θ0,P1,P2,n,L,Bt,STs
		
		α1 = αs[i]
		α2 = αs[j]
		θi = θ0 + α1*P1 + α2*P2
		θf,iter = kuramoto(L,zeros(n),θi, false)

		Qf[i,j] = winding(θf,Vector(1:n))
		
		J = jacobian(L,θi)
		λs = eigvals(J)

		m1,i1 = findmin(abs.(λs))
		Λs[i,j] = λs[[1:i1-1;i1+1:n]]

		ew = cos.(Bt*θi)
		Cst[i,j] = span_tree(ew,STs)
	end
end

