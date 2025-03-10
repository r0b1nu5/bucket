using PyPlot, DelimitedFiles

include("tools.jl")

figure("fig1", figsize=(12,4))

# =====================================================
subplot(1,2,1)

nums = [45,104,109,113,134,173,177,178,187,223,253]
for n in nums
	local x = readdlm("data-pj/fig1/fort.$n")
	PyPlot.plot(x[:,1]*1e-5,x[:,2])
end

axis([0,393,0,6])
xlabel("t [a.u.]")
ylabel("population")


# =====================================================
subplot(1,2,2)
zer0 = 1e-15
A = readdlm("data-pj/fig1/A.csv",',')
Si = 157
#surv = Int64.(readdlm("data-pj/fig1/survivingspecies_abundancies.dat")[3:end,1])[x[3:end,2] .> zer0]
#S = length(surv)
#A = A[surv,surv]
S = size(A)[1]
pα = 10

κ = 1.
μ = 5.
σ = 2.7
Id = diagm(0 => ones(S))

σs = LinRange(0,σ,1000)

λs = zeros(Complex{Float64},S,1)
λr = zeros(S,0)
λi = zeros(S,0)
for σ0 in σs
	local N0 = inv(Id + μ/Si*ones(S,S) + σ0/sqrt(Si)*A)*κ*ones(S)
	J,λ,us = analyze_jac(N0,A,κ,μ/Si,σ0/sqrt(Si))
	λ = associate_eigvals(λs[:,end],λ)
	global λs = [λs λ]
	global λr = [λr real(λ)]
	global λi = [λi imag(λ)]
end
λs = λs[:,2:end]

PyPlot.plot([0,0],maximum(imag.(λs))*[-1.2,1.2],"--k")
for i in 1:S
	plot_density(λr[i,:],λi[i,:],"C$(mod(i-1,10))",pα)
	PyPlot.plot(λr[i,end],λi[i,end],"o",color="C$(mod(i-1,10))")
end
xlabel("Re(λ)")
ylabel("Im(λ)")

