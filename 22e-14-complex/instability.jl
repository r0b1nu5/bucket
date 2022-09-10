using PyPlot, LinearAlgebra

Ns = Vector(20:100)
p = 2

qs = Vector{Int64}()

for N in Ns
	q = 0 
	Δ = 2π/N
	λ = -1000.
	μs = ones(1)

	while q < N/2 && λ < 1e-6
		q += 1
		λs = [2*sum([cos(k*p*Δ)*(cos(k*j*Δ) - 1) for k in 1:q]) for j in 0:N-1]
		λ = maximum(λs)
	end

	@info "N = $N is done."
	
	push!(qs,q)
end

figure("Fig 1")
PyPlot.plot(Ns,qs,"o",label="p = $p")



