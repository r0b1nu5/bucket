using PyPlot, Statistics, LinearAlgebra

ex = "02"

include("lv.jl")
include("tools.jl")

zer0 = 1e-10

N0 = vec(readdlm("data/ex-limit-cycle-"*ex*"-N1.csv",','))
S = length(N0)

A = readdlm("data/ex-limit-cycle-"*ex*"-A.csv",',')

κ = ones(S)

μ = 5.

σ = 2.7

n_iter = 200_000
for i in 1:10
	@info "i=$i"
	global N = lv_bunin(N0,A,κ,μ/S,σ/sqrt(S),n_iter,n_iter,1e-3)
	global N0 = N[:,end]
end
surv = vec(1:S)[N[:,end] .> zer0]

figure()
for s in surv
	PyPlot.plot(N[s,:])
end

Nm = vec(mean(N,dims=2))

J,λ,U = analyze_limit_cycle(N,A)
Jr,λr,Ur = analyze_limit_cycle(N[surv,:],A[surv,surv])






