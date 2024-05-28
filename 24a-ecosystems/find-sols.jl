using LinearAlgebra, PyPlot

include("tools.jl")
include("lv.jl")

S = 21
Si = 157

zer0 = 1e-15

κ = 1.
μ = 5.
σ = 2.7
Id = diagm(0 => ones(S))
c = 0

more_ext = Int64[]
stable = Int64[]
oscill = Int64[]

N0 = -ones(S)
A = zeros(S,S)
np = 0
λp = Float64[]

#for k in 1:100
k = 0
while np < 2
	global k += 1
	global N0 = -ones(S)
	global A = zeros(S,S)
	if k%10 == 0
		@info "iter: $k"
	end

	while minimum(N0) < 0.
		A = randn(S,S)
		N0 = inv(Id + μ/Si*ones(S,S) + σ/sqrt(Si)*A)*κ*ones(S)
	end
	
	J,λ,us = analyze_jac(N0,A,κ,μ/Si,σ/sqrt(Si))
	global np = sum(real.(λ) .> -zer0)
	global λp = λ[real.(λ) .> -zer0]

 #=
	if np == 0
		push!(stable,np)
	else
		N = lv_bunin(N0 + 1e-2*randn(S),A,κ,μ/Si,σ/sqrt(Si),150_000,100_000,1e-3)
		surv = vec(1:S)[N[:,end] .> zer0]
		#ampl = [(maximum(N[i,end-10_000:end]) - minimum(N[i,end-10_000:end])) for i in 1:S]

		if length(surv) < S
			push!(more_ext,np)
		else
			push!(oscill,np)
		end
	end
# =#
end

@info "======================="
@info "np = $np"
for λ0 in λp
	@info "$λ0 -- $(imag(λ0)/real(λ0))"
end

N = lv_bunin(N0 + 1e-2*randn(S),A,κ,μ/Si,σ/sqrt(Si),200_000,100_000,1e-3)
surv = vec(1:S)[N[:,end] .> zer0]

figure()
for s in 1:S
	PyPlot.plot(N[s,:])
end
title("np = $np, # survivors = $(length(surv))")





