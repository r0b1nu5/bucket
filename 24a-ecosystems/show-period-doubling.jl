using LinearAlgebra, PyPlot, DelimitedFiles, FFTW

include("tools.jl")
include("lv.jl")

S = 21
Si = 157
A = readdlm("data/ex-limit-cycle-04-A.csv",',')
N0 = vec(readdlm("data/ex-limit-cycle-04-N0.csv",','))

zer0 = 1e-15

κ = 1.
μ = 5.
σ = 2.7
nσ = 20
σs = LinRange(2.7,2.8,nσ)
Id = diagm(0 => ones(S))

#N0 = inv(Id + μ/Si*ones(S,S) + σ/sqrt(Si)*A)*κ*ones(S)
#J,λs,us = analyze_jac(N0,A,κ,μ/Si,σ/sqrt(Si))
#λ = associate_eigvals(λs[:,end],λ)

n_iter = 2_000_000
δt = 1e-3

for i in 1:nσ
	σ = σs[i]
	N = lv_bunin(N0,A,κ,μ/Si,σ/sqrt(Si),n_iter,n_iter,δt)
	F = fft(N[:,100_001:n_iter])
	writedlm("data/ex-limit-cycle-04-s$i-fft.csv",abs.(F[:,1:200]),',')
	global N0 = N[:,end]
end

include("plot-period-doubling.jl")




