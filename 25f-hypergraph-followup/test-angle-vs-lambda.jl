using PyPlot, Statistics

include("gen-hypergraph.jl")
include("hyper-kuramoto.jl")
include("tools.jl")

n = 10
p = .2
n_iter = 1
n_α = 50
α_max = 20

figure("Δ vs. λ")
PyPlot.plot([0,π],[0,0],"--k")
PyPlot.plot([π/2,π/2],[-10,2],"--k")
xlabel("Δ")
ylabel("λ2")

αs = LinRange(0,α_max,n_α)
Λ2 = Float64[]

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
	@info "α* = $(α_star)"
	global αs = LinRange(0,α_star*1.2,n_α)

        λ2old = -Inf
        λ2 = 0.
	for α in αs
		@info "α/α* = $(α/α_star)"
		Θs,dΘs = hyper_k(zeros(0,3),A,α*ω,θ0,0.,0.,.01,100000,1e-3)
		
                global θ0 = Θs[:,end]
		J = get_jac(θ0,A)
		λs = eigvals(J)
                λ2 = sort(real.(λs))[end-1]

                if λ2[end] > 0 || λ2 < λ2old
                    continue
                end
                λ2old = copy(λ2)
	end
        push!(Λ2,λ2)

	figure("Δ vs. λ")
        PyPlot.plot(maximum(abs.(B2*Θs[:,end])),λ2,"o",color="C3")

end



