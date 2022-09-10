using PyPlot

N = 100
Δ = 2π/N

for p in 1:3
	Λs = Vector{Float64}()
	qs = Int64.(1:N/2)
	
	for q in qs
		λs = [2*sum([cos(k*p*Δ)*(cos(k*j*Δ) - 1) for k in 1:q]) for j in 0:N-1]
		push!(Λs,maximum(λs[2:end]))
	end
	
	subplot(1,3,p)
	PyPlot.plot(qs./N,Λs./N,"o",label="N = $N, p = $p")
end






