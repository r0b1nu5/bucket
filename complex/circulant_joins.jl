using LinearAlgebra, PyPlot

include("L2B.jl")
include("cycle.jl")

 #=
ns = [6,6,6,6,6,6]
ks = [1,1,1,1,1,1]
# =#
# #=
ns = [6,5]
ks = [1,1]
# =#

n0 = length(ns)
N = sum(ns)

As = [cycle(ns[i],ks[i]) for i in 1:n0]

#A0 = cycle(n0)
A0 = [0 1;1 0.]

#qs = [1,1,1,1,1,1]
qs = [1,-1]

θs = [Vector(qs[i]*2π*(1:ns[i])/ns[i]) for i in 1:n0]
θ0 = Vector{Float64}()
for i in 1:n0
	global θ0 = [θ0;θs[i]]
end

A = zeros(N,N)
for i in 1:n0-1
	i0 = sum(ns[1:i-1])
	for j in i+1:n0
		j0 = sum(ns[1:j-1])
		A[i0+1:i0+ns[i],j0+1:j0+ns[j]] = A0[i,j]*ones(ns[i],ns[j])
		A[j0+1:j0+ns[j],i0+1:i0+ns[i]] = A0[j,i]*ones(ns[j],ns[i])
	end
	A[i0+1:i0+ns[i],i0+1:i0+ns[i]] = As[i]
end
A[end-ns[end]+1:end,end-ns[end]+1:end] = As[end]
D = diagm(0 => vec(sum(A,dims=1)))
L = D - A
B,w = L2B(L)





dot = B*sin.(B'*θ0)
@info "Max. deriv.: $(maximum(abs.(dot)))"


dθ = θ0*ones(1,N) - ones(N)*θ0'
preJ = A.*cos.(dθ)
DJ = diagm(0 => vec(sum(preJ,dims=1)))
J = preJ - DJ
λs = eigvals(J)
@info "Max. eigval: $(maximum(λs))"


