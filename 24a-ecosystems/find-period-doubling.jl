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

 #=
ex = "pj4"
S = 28
Si = 157
A = readdlm("data/ex-limit-cycle-pj3-A.csv",',')
N0 = vec(readdlm("data/ex-limit-cycle-pj3-N0.csv",','))
N0 = vec(readdlm("data/ex-limit-cycle-pj4-s17-N0.csv",','))
# =#

 #=
ex = "pj5"
S = 42
Si = 157
A = readdlm("data/ex-limit-cycle-"*ex*"-A.csv",',')
N0 = vec(readdlm("data/ex-limit-cycle-"*ex*"-N0.csv",','))
#N0 = vec(readdlm("data/ex-limit-cycle-"*ex*"-s2-N0.csv",','))
nσ = 12 
σs = [2.2,2.3,2.4,2.5,2.55,2.6,2.65,2.7,2.75,2.8,2.9,3.0]
# =#
 #=
ex = "pj6"
S = 42
Si = 157
A = readdlm("data/ex-limit-cycle-"*ex*"-A.csv",',')
N0 = vec(readdlm("data/ex-limit-cycle-"*ex*"-N0.csv",','))
#N0 = vec(readdlm("data/ex-limit-cycle-"*ex*"-s2-N0.csv",','))
nσ = 23 
σs = [2.2,2.3,2.4,2.5,2.55,2.6,2.62,2.64,2.66,2.68,2.7,2.72,2.74,2.76,2.78,2.8,2.82,2.84,2.86,2.88,2.9,2.95,3.0]
# =#
 #=
ex = "pj7"
S = 42
Si = 157
A = readdlm("data/ex-limit-cycle-"*ex*"-A.csv",',')
N0 = vec(readdlm("data/ex-limit-cycle-"*ex*"-N0.csv",','))
#N0 = vec(readdlm("data/ex-limit-cycle-"*ex*"-s8-N0.csv",','))
nσ = 23 
σs = [2.2,2.3,2.4,2.5,2.55,2.6,2.62,2.64,2.66,2.68,2.7,2.72,2.74,2.76,2.78,2.8,2.82,2.84,2.86,2.88,2.9,2.95,3.0]
# =#
 #=
ex = "pj8"
i1 = 3
S = 42
Si = 157
A = readdlm("data/ex-limit-cycle-"*ex*"-A.csv",',')
if i1 == 1
	N0 = vec(readdlm("data/ex-limit-cycle-"*ex*"-N0.csv",','))
else
	N0 = vec(readdlm("data/ex-limit-cycle-"*ex*"-s$i1-N0.csv",','))
end
nσ = 41
σs = LinRange(2.64,2.84,41)
# =#
# #=
ex = "pj9"
i1 = 16
S = 88
Si = 157
A = readdlm("data/ex-limit-cycle-"*ex*"-A.csv",',')
if i1 == 1
	N0 = vec(readdlm("data/ex-limit-cycle-"*ex*"-N0.csv",','))
else
	N0 = vec(readdlm("data/ex-limit-cycle-"*ex*"-s$i1-N0.csv",','))
end
nσ = 26
σs = collect(2.7:.02:3.2)
# =#

κ = 1.
μ = 5.
#nσ = 191
#σs = LinRange(2.7,3.4,nσ)
#nσ = 9
#σs = [2.,2.32,2.35,2.384,2.417,2.45,2.7,2.776,2.857]
#nσ = 20
#σs = [2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.65,2.7,2.75,2.8,2.85,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6]

Id = diagm(0 => ones(S))


n_iter = 2_000_000
δt = 1e-3
Δt = 1e-1
n_step = round(Int64,n_iter*δt/Δt)

N = zeros(S,0)

for i in i1:nσ
	writedlm("data/ex-limit-cycle-"*ex*"-s$i-N0.csv",N0,',')
	σ = σs[i]
	@info "s$i: σ = $σ/$(maximum(σs))"
	global N = lv_bunin(N0,A,κ,μ/Si,σ/sqrt(Si),n_iter,n_iter,δt,Δt,-1.)
	F = fft(N[1,200:n_step-1])
	writedlm("data/ex-limit-cycle-"*ex*"-s$i-N1.csv",N[1,:],',')
	writedlm("data/ex-limit-cycle-"*ex*"-s$i-fft.csv",abs.(F[1:400]),',')
	global N0 = N[:,end]
end

include("plot-period-doubling.jl")




