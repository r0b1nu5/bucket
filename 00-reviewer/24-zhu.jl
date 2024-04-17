using PyPlot, Statistics, LinearAlgebra

N = 10 # Number of vertices

m = .0001 # Inertia
Ω = rand(N) # Natural frequencies
Ω .-= mean(Ω)
K = 1000000. # Coupling strength
α = .001 # Frustration

# Building the adjacency matrix
A = zeros(N,N)
p = .5 # Probability that an edge exists
a = 1. # Scaling of the edges
as = Float64[]
for i in 1:N-1
	for j in i+1:N
		if rand() > p
			A[i,j] = a*rand() 
			A[j,i] = A[i,j]
			push!(as,A[i,j])
		end
	end
end
D = diagm(0 => vec(sum(A,dims=1)))
L = D - A

# Quantities depending on the system
DΩ = sum((Ω*ones(1,N) - ones(N)*Ω').^2)
au = maximum(as)
al = minimum(as)

E = sum(A .> 1e-8)/2 # Number of edges
Ec = N*(N-1)/2 - E # Number of absent edges
r = 1 # longest shortest path, so 1 is the best case, i.e., complete graph.

# Computing the energy E1
function E1(θ,ω,m,K)
	N = length(θ)

	θθ = θ*ones(1,N) - ones(N)*θ'
	ωω = ω*ones(1,N) - ones(N)*ω'

	return 2*m*sum(θθ.*ωω) + (1-m*sqrt(K))*sum(θθ.^2) + 2*m^2*sum(ωω.^2)
end

# Bounds take (very) small...
D0 = .01 # in (0,π)
Di = .001 # in (0,π/2), needed small otherwise rhs of (2.4)-1 gets too small.

θ = .001*rand(N)
ω = .001*rand(N)

C1 = ((3*(1+r*Ec)*D0)/(K^(3/2)*al*cos(α)*sin(D0)) + 12*m/sqrt(K))*DΩ + 12*N^2*sqrt(K)*au^2*(sin(α))^2*(((1+r*Ec)*D0)/(al*cos(α)*sin(D0)) + 4*m*K)

@info "Checking conditions..."

function ro(x,d=4)
	return round(x,sigdigits=d)
end

x1 = E1(θ,ω,m,K)
x2 = D0^2/8
if x1 > x2
	@info "(2.3)-1 is not satisfied: $(ro(x1)) > $(ro(x2))"
end
if Di > D0
	@info "(2.3)-4 is not satisfied: $(ro(Di)) > $(ro(D0))"
end
x3 = sin(α)
x4 = (al*cos(α)*cos(Di))/(12*au*(1+r*Ec)*sin(Di))
if x3 > x4
	@info "(2.4)-1 is not satisfied: $(ro(x3)) > $(ro(x4))"
end
if m >= 1.
	@info "(2.4)-2 is not statisfied: $m > 1"
end
x5 = sqrt(K) 
x6 = (1+r*Ec)*max(sin(D0)/(36*D0),1/cos(Di))/(al*cos(α))
if x5 < x6
	@info "(2.5) is not satisfied: $(ro(x5)) < $(ro(x6))"
end
if m*sqrt(K) > 1/8
	@info "(2.6)-1 is not satisfied: $(ro(m*sqrt(K))) > $(1/8)"
end
x7 = m*K
x8 = al*min(sin(D0)/(36*D0),cos(Di)/24)/(au^2*(1+r*Ec)*cos(α))
if x7 > x8
	@info "(2.6)-2 is not satisfied: $(ro(x7)) > $(ro(x8))"
end
x9 = K*sin(α)
x10 = 1/(8*au*sin(Di))
if x9 > x10
	@info "(2.7) is not satisfied: $(ro(x9)) > $(ro(x10))"
end
if C1 > Di^2/16
	@info "(2.8) is not satisfied: $(ro(C1)) > $(ro(Di^2/16))"
end








