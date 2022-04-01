using LinearAlgebra, PyPlot

include("L2B.jl")

n1 = 12
k1 = 1
n2 = 12
k2 = 1
n = n1 + n2

q1 = 1
q2 = -1

θ1 = Vector(q1*2π*(1:n1)./n1) .+ 2*rand() .-1
θ2 = Vector(q2*2π*(1:n2)./n2) .+ 2*rand() .-1
#θ2 = copy(θ1) .+ π
θ = [θ1;θ2]

A1 = zeros(n1,n1)
for k in 1:k1
	global A1 += diagm(k => ones(n1-k)) + diagm(-k => ones(n1-k)) + diagm(n1-k => ones(k)) + diagm(k-n1 => ones(k))
end
D1 = diagm(0 => vec(sum(A1,dims=1)))
L1 = D1 - A1

A2 = zeros(n2,n2)
for k in 1:k2
	global A2 += diagm(k => ones(n2-k)) + diagm(-k => ones(n2-k)) + diagm(n2-k => ones(k)) + diagm(k-n2 => ones(k))
end
D2 = diagm(0 => vec(sum(A2,dims=1)))
L2 = D2 - A2


A = [A1 ones(n1,n2);ones(n2,n1) A2]
#A = [A1 diagm(0 => ones(n1));diagm(0 => ones(n1)) A2]
#A = [A1 A1;A2 A2]
D = diagm(0 => vec(sum(A,dims=1)))
L = D - A
B,w = L2B(L)

dot = B*sin.(B'*θ)
@info "Max. deriv.: $(maximum(abs.(dot)))"


dθ = θ*ones(1,n) - ones(n)*θ'
preJ = A.*cos.(dθ)
DJ = diagm(0 => vec(sum(preJ,dims=1)))
J = preJ - DJ
λs = eigvals(J)
@info "Max. eigval: $(maximum(λs))"


