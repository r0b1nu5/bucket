using PyPlot

include("kuramoto.jl")
include("cycle.jl")
include("splay_states.jl")
include("jacobian.jl")
include("span_tree.jl")
include("rand_idx.jl")

n = 83
q = 15

id = 1
#id = 2
#id = 3
#id = 4

res1 = 20
res2 = 20

L = cycle(n)
B,w,Bt = L2B(L)
θ0 = splay_q(q,n)

k = floor(Int64,n/2)
P1 = zeros(n)
P1[rand_idx(n,k)] = ones(Int64,k)
P2 = zeros(n)
P2[rand_idx(n,k)] = ones(Int64,k)

#P1 = [ones(41);zeros(42)]
#P2 = [zeros(42);ones(41)]

αs = LinRange(-π,π,res2)

STs = cycle_span_trees(n)

Qf = zeros(Int64,res1,res2)
Λs = Vector{Matrix{Float64}}()
for i in 1:n-1
	push!(Λs,zeros(res1,res2))
end
Cst = zeros(res1,res2)

for i in 1:res1
	for j in 1:res2
		@info "(i,j) = ($i,$j)"
		
		global Λs,Qf,αs,θ0,P1,P2,n,L,Bt,STs
		
		α1 = αs[(id-1)*res1 + i]
		α2 = αs[j]
		θi = θ0 + α1*P1 + α2*P2
		θf,iter = kuramoto(L,zeros(n),θi, false)

		Qf[i,j] = winding(θf,Vector(1:n))
		
		J = jacobian(L,θi)
		λs = eigvals(J)

		m1,i1 = findmin(abs.(λs))
		ids = [1:i1-1;i1+1:n]
		for k in 1:n-1
			idx = ids[k]
			Λs[k][i,j] = λs[idx]
		end

		ew = cos.(Bt*θi)
		Cst[i,j] = span_tree(ew,STs)
	end
end

 #=
writedlm("temp/Qf_$(id).csv",Qf,',')
writedlm("temp/Cs_$(id).csv",Cst,',')
for i in 1:n-1
	writedlm("temp/L$(i)_$(id).csv",Λs[i],',')
end
# =#
