using LinearAlgebra, PyPlot

include("tools.jl")
include("lv.jl")


# #= 
S = 21
Si = 157

zer0 = 1e-15

κ = 1.
μ = 5.
σ = 2.7
Id = diagm(0 => ones(S))

N0 = -ones(S)
A = zeros(S,S)
np = 0
λp = Float64[]

#for k in 1:100
k = 0
while np < 3
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

end

@info "======================="
@info "np = $np"
for λ0 in λp
	@info "$λ0 -- $(imag(λ0)/real(λ0))"
end

# =#

#σs = LinRange(0,σ,1000)
σs = log.(LinRange(exp(0),exp(σ),1000))

λs = zeros(Complex{Float64},S,1)
λr = zeros(S,0)
λi = zeros(S,0)
for σ0 in σs
	N0 = inv(Id + μ/Si*ones(S,S) + σ0/sqrt(Si)*A)*κ*ones(S)
	J,λ,us = analyze_jac(N0,A,κ,μ/Si,σ0/sqrt(Si))
	λ = associate_eigvals(λs[:,end],λ)
	global λs = [λs λ]
	global λr = [λr real(λ)]
	global λi = [λi imag(λ)]
	#λt = sortslices([imag.(λ) real.(λ)],dims=1)
	#global λr = [λr λt[:,2]]
	#global λi = [λi λt[:,1]]
end
λs = λs[:,2:end]

figure()
PyPlot.plot([0,0],maximum(imag.(λs))*[-1.2,1.2],"--k")
for i in 1:S
	plot_density(λr[i,:],λi[i,:],"C$(mod(i-1,10))")
#	PyPlot.plot(λr[i,:],λi[i,:],color="C$(mod(i-1,10))")
#	PyPlot.plot(λr[i,1],λi[i,1],">",color="C$(mod(i-1,10))")
	PyPlot.plot(λr[i,end],λi[i,end],"o",color="C$(mod(i-1,10))")
end
xlabel("Re(λ)")
ylabel("Im(λ)")
title("σ = $(round(maximum(σs),digits=2))")




