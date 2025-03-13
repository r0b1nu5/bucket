using LinearAlgebra, PyPlot, DelimitedFiles

include("tools.jl")
include("lv.jl")


 #= 
S = 21
Si = 157

zer0 = 1e-15

κ = 1.
μ = 5.
σ = 2.7
#σ = 10.
Id = diagm(0 => ones(S))

N0 = -ones(S)
A = zeros(S,S)
np = 0
λp = Float64[]

k = 0
while np < 2 && k < 1000
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

# #= 
#x = readdlm("data/matr-aij-02.dat")
#S = 42
#Si = 157
#pα = 1

#x = readdlm("data/Matr-240825.dat")
#S = 38
#Si = 157
#pα = 10

#x = readdlm("data/SMatr-240825.dat")
#S = 28
#Si = 157
#pα = 10

 #=
A = zeros(S,S)
for k in 1:size(x)[1]
	i = Int64(x[k,1])
	j = Int64(x[k,2])
	A[i,j] = x[k,3]
end
# =#

 #=
A = readdlm("data/ex-limit-cycle-pj9-A.csv",',')
S = 89
Si = 157
pα = 10
# =#

 #= Fig. 1
zer0 = 1e-14
A = readdlm("data-pj/fig1/A.csv",',')
Si = 157
surv = Int64.(readdlm("data-pj/fig1/survivingspecies_abundancies.dat")[3:end,1])[x[3:end,2] .> zer0]
S = length(surv)
A = A[surv,surv]
pα = 10
# =#

# #=
file2σ = Dict{String,Tuple{Float64,Int64}}("matr2n" => (4.2,3),
					   "matr5n" => (4.02,3),
					   "matr7n" => (3.34,4),
					   "matr8n" => (3.14,3),
					   "matr10n" => (2.49,4),
					   "matr11n" => (3.04,4),
					   "matr14f" => (4.24,5))
file = "matr2n"
zer0 = 1e-12
Si = 157
x = readdlm("data-pj/fig1/"*file*"/matr.dat")
A = zeros(Si,Si)
for l in 1:size(x)[1]
	i = Int64(x[l,1])
	j = Int64(x[l,2])
	A[i,j] = x[l,3]
end
y = readdlm("data-pj/fig1/"*file*"/survivingspecies_abundancies.dat")
k = file2σ[file][2]
surv = Int64.(y[k:end,1])[y[k:end,2] .> zer0]
A = A[surv,surv]
S = length(surv)
pα = 10
# =#

κ = 1.
μ = 5.
σ = 2.7
σ = file2σ[file][1]
Id = diagm(0 => ones(S))


# =#


σs = LinRange(0,σ,1000)
#σs = log.(LinRange(exp(0),exp(1.45*σ),1000))

λs = zeros(Complex{Float64},S,1)
λr = zeros(S,0)
λi = zeros(S,0)
for σ0 in σs
	N0 = inv(Id + μ/Si*ones(S,S) + σ0/sqrt(Si)*A)*κ*ones(S)
	# @info "min(N0) = $(minimum(N0))"
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

figure("fig1-"*file)
PyPlot.plot([0,0],maximum(imag.(λs))*[-1.2,1.2],"--k")
for i in 1:S
	plot_density(λr[i,:],λi[i,:],"C$(mod(i-1,10))",pα)
#	PyPlot.plot(λr[i,:],λi[i,:],color="C$(mod(i-1,10))")
#	PyPlot.plot(λr[i,1],λi[i,1],">",color="C$(mod(i-1,10))")
	PyPlot.plot(λr[i,end],λi[i,end],"o",color="C$(mod(i-1,10))")
end
xlabel("Re(λ)")
ylabel("Im(λ)")
title("σ = $(round(maximum(σs),digits=2))")




