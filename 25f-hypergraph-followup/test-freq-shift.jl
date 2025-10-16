using PyPlot, Statistics

include("gen-hypergraph.jl")
include("hyper-kuramoto.jl")
include("tools.jl")

n = 10
p = .2
n_α = 50
α_max = 20

αs = LinRange(0,α_max,n_α)

A,B,B2,E = rand_3_graph(n,p)
W = diagm(0 => repeat(A[:,4],inner=3))
ω = rand(n); ω .-= mean(ω)

α_star = 1/maximum(inv(W)*pinv(B)*ω)
@info "α* = $(α_star)"
global αs = LinRange(0,α_star*1.2,n_α)

λ2old = -Inf
λ2 = -Inf
θ0old = zeros(n)
θ0 = zeros(n)

αr = 0.

for α in αs
    @info "α/α* = $(α/α_star)"
    global Θs,dΘs,iter = hyper_k(zeros(0,3),A,α*ω,θ0,0.,0.,.01,100000,1e-3)

    global λ2old = copy(λ2)
    global θ0old = copy(θ0)
    global θ0 = Θs[:,end]
    J = get_jac(θ0,A)
    λs = eigvals(J)
    global λ2 = sort(real.(λs))[end-1]

    @info "λ2(t-1) = $λ2old, λ2(t) = $λ2"
    if λ2[end] > 0 || λ2 < λ2old || iter == 100000
        global αr = α
        break
    end
end

αs = Float64[]
ρs = LinRange(-10,10,100)
for ρ in ρs
    x = 1/maximum(inv(W)*pinv(B)*(ω .+ ρ))
    push!(αs,x)
end

figure("fig")
PyPlot.plot([-10,10],[αr,αr])
PyPlot.plot(ρs,αs)
xlabel("ρ")
ylabel("α")





