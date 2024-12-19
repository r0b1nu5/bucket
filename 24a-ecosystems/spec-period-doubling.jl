using LinearAlgebra, PyPlot

include("tools.jl")
include("lv.jl")

 #=
ex = "pj7"
S = 42
Si = 157
A = readdlm("data/ex-limit-cycle-"*ex*"-A.csv",',')
N0 = vec(readdlm("data/ex-limit-cycle-"*ex*"-N0.csv",','))

T = 99999 - 1999
δt = .1
ωs = (0:T-1)*2π/(T*δt)
nσ = 23 
σs = [2.2,2.3,2.4,2.5,2.55,2.6,2.62,2.64,2.66,2.68,2.7,2.72,2.74,2.76,2.78,2.8,2.82,2.84,2.86,2.88,2.9,2.95,3.0]
# =#
# #=
ex = "pj8"
S = 42
Si = 157
A = readdlm("data/ex-limit-cycle-"*ex*"-A.csv",',')
N0 = vec(readdlm("data/ex-limit-cycle-"*ex*"-N0.csv",','))

T = 99999 - 1999
δt = .1
ωs = (0:T-1)*2π/(T*δt)
nσ = 41
σs = LinRange(2.64,2.84,41)
# =#

zer0 = 1e-15

κ = 1.
μ = 5.
σ = 2.7
Id = diagm(0 => ones(S))

λs = zeros(Complex{Float64},S,1)
Λs = Dict{Float64,Vector{Complex{Float64}}}()
λr = zeros(S,0)
λi = zeros(S,0)
for i in 1:nσ-1
	n0 = vec(readdlm("data/ex-limit-cycle-"*ex*"-s$(i+1)-N0.csv",','))
	ids = (1:S)[n0 .> zer0]
@info "$(length(ids))"
	s = length(ids)
	a = A[ids,ids]
	id = diagm(0 => ones(s))
	σ0 = σs[i]
	N0 = inv(id + μ/Si*ones(s,s) + σ0/sqrt(Si)*a)*κ*ones(s)
	co1 = minimum(N0) < -.01 ? "C2" : "C0"
	co2 = minimum(N0) < -.01 ? "C3" : "C1"
	# @info "min(N0) = $(minimum(N0))"
	J,λ,us = analyze_jac(N0,a,κ,μ/Si,σ0/sqrt(Si))
	#λ = associate_eigvals(λs[:,end],λ)
	#global λs = [λs λ]
	#global λr = [λr real(λ)]
	#global λi = [λi imag(λ)]
	#λt = sortslices([imag.(λ) real.(λ)],dims=1)
	#global λr = [λr λt[:,2]]
	#global λi = [λi λt[:,1]]
	Λs[σ0] = λ
	λ1 = λ[real.(λ) .< 0]
	λ2 = λ[real.(λ) .> 0]

	figure("Spectra")
	subplot(ceil(Int64,nσ/3),3,i)
	PyPlot.plot([0,0],[-100,100],"--k")
	PyPlot.plot(real.(λ1),imag.(λ1),".",color=co1,label="$(round(σ0,digits=2))")
	PyPlot.plot(real.(λ2),imag.(λ2),".",color=co2)
	legend()
	xlabel("Re(λ)")
	ylabel("Im(λ)")
	axis([-3,2,-3,3])
end

nλ = maximum([length(Λs[s]) for s in keys(Λs)])
λs = zeros(nλ,nσ-1)
for i in 1:nσ-1
	λ = Λs[σs[i]]
	λs[1:length(λ),i] = imag.(λ)
end

figure()
for i in 1:nλ
	PyPlot.plot(σs[2:end],λs[i,:])
end





#=
figure()
PyPlot.plot([0,0],maximum(imag.(λs))*[-1.2,1.2],"--k")
for i in 1:S
	plot_density(λr[i,:],λi[i,:],"C$(mod(i-1,10))",pα)
#	PyPlot.plot(λr[i,:],λi[i,:],color="C$(mod(i-1,10))")
#	PyPlot.plot(λr[i,1],λi[i,1],">",color="C$(mod(i-1,10))")
	PyPlot.plot(λr[i,end],λi[i,end],"o",color="C$(mod(i-1,10))")
end
=#



