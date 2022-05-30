using LinearAlgebra

n = 10
P = diagm(0 => ones(n)) - ones(n)*ones(1,n)/n
V = eigvecs(P)[:,2:n]
R = V'

test = true
c = 0

while test && c < 1000000
	global c,n,P,V,R

	c += 1
	if c%100 == 0
		@info "c = $c"
	end

	θ = π/2*rand(n)
	As = sin.(θ*ones(1,n) - ones(n)*θ')
	Ds = diagm(0 => vec(sum(As,dims=2)))
	Ls = Ds - As
	Ms = (P*Ls + Ls'*P)/2
	
	λs = eigvals(Ls)
	μs = eigvals(Ms)

	test = (minimum(Ds) <= minimum(λs))*(maximum(Ds) >= maximum(λs))*(minimum(Ds) <= minimum(μs))*(maximum(Ds) >= maximum(μs))
end




