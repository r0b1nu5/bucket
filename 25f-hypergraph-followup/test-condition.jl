using PyPlot, Statistics

include("gen-hypergraph.jl")
include("hyper-kuramoto.jl")
include("tools.jl")

n = 10
p = .2
n_iter = 1
n_α = 50
α_max = 20

figure("Condition")
subplot(2,1,1)
PyPlot.plot([0,1.05],[0,0],"--k")
xlabel("α/α_star")
ylabel("λ2")
subplot(2,1,2)
PyPlot.plot([0,1.05],[π/2,π/2],"--k")
xlabel("α/α_star")
ylabel("Δ")

αs = LinRange(0,α_max,n_α)
Λ1 = zeros(n_α,0)
Λ2 = zeros(n_α,0)
Δs = zeros(n_α,0)

θ0 = zeros(n)
xmin = 0.
xmax = 0.
ymin = 0.
ymax = 0.

for i in 1:n_iter
	A,B,B2,E = rand_3_graph(n,p)
	W = diagm(0 => repeat(A[:,4],inner=3))
	ω = rand(n); ω .-= mean(ω)

        b = hyper2edge(A)
	
	α_star = 1/maximum(inv(W)*pinv(B)*ω)
        α_bis = 1/maximum(abs.(b'*pinv(b*b')*ω))

	@info "α* = $(α_star)"
	global αs = LinRange(0,α_star*1.01,n_α)

	λ1 = Float64[]
	λ2 = Float64[]
	Δ = Float64[]
	for α in αs
		@info "α/α* = $(α/α_star)"
		Θs,dΘs = hyper_k(zeros(0,3),A,α*ω,θ0,0.,0.,.01,100000,1e-3)
		push!(Δ,maximum(abs.(B2'*Θs[:,end])))

		global θ0 = Θs[:,end]
		J = get_jac(θ0,A)
		λs = eigvals(J)
		push!(λ1,sort(real.(λs))[end])
		push!(λ2,sort(real.(λs))[end-1])

	end
	global Λ1 = [Λ1 λ1]
	global Λ2 = [Λ2 λ2]
	global Δs = [Δs Δ]

	figure("Condition")
	subplot(2,1,1)
#	PyPlot.plot(αs./α_star,Λ1[:,end],":k")
	PyPlot.plot(αs./α_star,Λ2[:,end],"--o",color="C$i")
        PyPlot.plot(α_bis/α_star*[1,1],[-10,10],"--",color="C$i")
	global xmin = minimum([xmin;αs./α_star])
	global xmax = maximum([xmax;αs./α_star])
	global ymin = minimum([ymin;λ2])
	global ymax = maximum([ymax;λ2])
	subplot(2,1,2)
	PyPlot.plot(αs./α_star,Δ,"--o")
end


figure("Condition")
subplot(2,1,1)
PyPlot.plot([1,1],[ymin,ymax],"--k")


