using PyPlot, Statistics

include("gen-hypergraph.jl")
include("hyper-kuramoto.jl")
include("tools.jl")

n = 10
p = .2
n_iter = 3
n_α = 10
α_max = 20

αs = LinRange(0,α_max,n_α)
Λ1 = zeros(n_α,0)
Λ2 = zeros(n_α,0)

θ0 = zeros(n)
xmin = 0.
xmax = 0.
ymin = 0.
ymax = 0.
for i in 1:n_iter
	A,B,E = rand_3_graph(n,p)
	W = diagm(0 => repeat(A[:,4],inner=3))
	ω = rand(n); ω .-= mean(ω)
	
	α_star = 1/maximum(inv(W)*pinv(B)*ω)
	@info "$(α_star)"

	λ1 = Float64[]
	λ2 = Float64[]
	for α in αs
		Θs,dΘs = hyper_k(zeros(0,3),A,α*ω,θ0,0.,0.,.01,100000,1e-4)
		global θ0 = Θs[:,end]
		J = get_jac(θ0,A)
		λs = eigvals(J)
		push!(λ1,sort(real.(λs))[end])
		push!(λ2,sort(real.(λs))[end-1])
	end
	global Λ1 = [Λ1 λ1]
	global Λ2 = [Λ2 λ2]

	figure("Condition")
	PyPlot.plot(αs./α_star,Λ2[:,end],"--o")
	global xmin = minimum([xmin;αs./α_star])
	global xmax = maximum([xmax;αs./α_star])
	global ymin = minimum([ymin;λ2])
	global ymax = maximum([ymax;λ2])
end


figure("Condition")
PyPlot.plot([xmin,xmax],[0,0],"--k")
PyPlot.plot([1,1],[ymin,ymax],"--k")
xlabel("α/α_star")
ylabel("λ2")


