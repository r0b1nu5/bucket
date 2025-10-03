using PyPlot, Statistics

include("gen-hypergraph.jl")
include("hyper-kuramoto.jl")
include("tools.jl")

n = 10
p = .2
n_iter = 3
n_α = 3
α_max = 5

αs = LinRange(0,α_max,n_α)
Λ2 = zeros(n_α,0)

θ0 = zeros(n)

for i in 1:n_iter
	A,B,E = rand_3_graph(n,p)



	α_star = ...
	


	ω = rand(n); ω .-= mean(ω)
	ω = zeros(n)
	λ2 = Float64[]
	for α in αs
		Θs,dΘs = hyper_k(zeros(0,3),A,α*ω,θ0)
		θ0 = Θs[:,end]
		J = get_jac(θ0,A)
		push!(λ2,eigvals(J)[end-1])
	end
	Λ2 = [Λ2 λ2]

	figure("Condition")
	PyPlot.plot(αs,Λ2[:,end],"--o")
end


figure("Conditions")
PyPlot.plot(αs,zeros(n_α),"--k")
PyPlot.plot(
xlabel("α")
ylabel("λ2")


