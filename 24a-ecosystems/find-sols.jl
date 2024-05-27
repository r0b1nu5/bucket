using LinearAlgebra, PyPlot

include("tools.jl")
include("lv.jl")

S = 21


κ = 1.
μ = 5.
σ = 2.7
Id = diagm(0 => ones(S))
N0 = -ones(S)
A = zeros(S,S)

while minimum(N0) < -1e-11
	global A = randn(S,S)
	global N0 = inv(Id + μ/S*ones(S,S) + σ/sqrt(S)*A)*κ*ones(S)
end

#f = f_lv_bunin(N0,A,κ*ones(S),μ,σ)
#@info "max|Ṅ| = $(maximum(abs.(f)))"

J,λ,us = analyze_limit_cycle([N0 N0],A)
@info "max(λ) = $(maximum(real(λ)))"
@info "# pos. eig. vals: $(sum(real.(λ) .> 0.))"

N = lv_bunin(N0 + 1e-4*randn(S),A,κ,μ/S,σ/sqrt(S),200_000,100_000,1e-3)


figure()
for s in 1:S
	PyPlot.plot(N[s,:])
end






