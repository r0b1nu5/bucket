using PyPlot, DelimitedFiles

include("tools.jl")

figure("fig1", figsize=(12,4))

# =====================================================
subplot(1,2,1)

nums = [107,110,113,117,122,133,141,142,146,153,154,157,162,166,174,176,190,198,206,213,220,225,227,235,239,243,250,254]
for n in nums
	local x = readdlm("data-pj/fort.$n")
	PyPlot.plot(x[:,1],x[:,2])
end

axis([0,100,0,16])
xlabel("t [a.u.]")
ylabel("population")


# =====================================================
subplot(1,2,2)
x = readdlm("data/SMatr-240825.dat")
S = 28
Si = 157
pα = 10

A = zeros(S,S)
for k in 1:size(x)[1]
	i = Int64(x[k,1])
	j = Int64(x[k,2])
	A[i,j] = x[k,3]
end

zer0 = 1e-15

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

