using LinearAlgebra, PyPlot, DelimitedFiles, FFTW

include("tools.jl")
include("lv.jl")

 #=
ex = "04"
S = 21
Si = 157
A = readdlm("data/ex-limit-cycle-04-A.csv",',')
N0 = vec(readdlm("data/ex-limit-cycle-04-N0.csv",','))
N0 = vec(readdlm("data/ex-limit-cycle-04-s134-N0.csv",','))

κ = 1.
μ = 5.
nσ = 191
σs = LinRange(2.7,3.7,nσ)
Id = diagm(0 => ones(S))
# =#

# #=
ex = "pj3"
S = 28
Si = 157
A = readdlm("data/ex-limit-cycle-pj3-A.csv",',')
N0 = vec(readdlm("data/ex-limit-cycle-pj3-N0.csv",','))
N0 = vec(readdlm("data/ex-limit-cycle-pj3-s7-N0.csv",','))

κ = 1.
μ = 5.
#nσ = 191
#σs = LinRange(2.7,3.4,nσ)
#nσ = 9
#σs = [2.,2.32,2.35,2.384,2.417,2.45,2.7,2.776,2.857]
nσ = 14
σs = [2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.65,2.7,2.75,2.8,2.85,2.9,3.0]
Id = diagm(0 => ones(S))
# =#

n_iter = 2_000_000
δt = 1e-3
Δt = 1e-1
n_step = round(Int64,n_iter*δt/Δt)

#for i in 1:nσ
for i in 7:nσ
	σ = σs[i]
	@info "σ = $σ/$(maximum(σs))"
	N = lv_bunin(N0,A,κ,μ/Si,σ/sqrt(Si),n_iter,n_iter,δt,Δt,-1.)
	F = fft(N[1,2000:n_step-1])
	writedlm("data/ex-limit-cycle-"*ex*"-s$i-N0.csv",N0,',')
	writedlm("data/ex-limit-cycle-"*ex*"-s$i-N1.csv",N[1,:],',')
	writedlm("data/ex-limit-cycle-"*ex*"-s$i-fft.csv",abs.(F[1:400]),',')
	global N0 = N[:,end]
end

include("plot-period-doubling.jl")




