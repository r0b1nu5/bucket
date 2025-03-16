using PyPlot, DelimitedFiles

include("tools.jl")

file2σ = Dict{String,Tuple{Float64,Int64}}("matr2n/" => (4.2,3),
					   "matr5n/" => (4.02,3),
					   "matr7n/" => (3.34,4),
					   "matr8n/" => (3.14,3),
					   "matr10n/" => (2.49,4),
					   "matr11n/" => (3.04,4),
					   "matr14f/" => (4.24,5))

 #=
file = ""
nums = [45,104,109,113,134,173,177,178,187,223,253]
# =#
 #=
file = "matr2n/"
nums = [102,103,110,118,120,132,135,168]
σ = 4.2
k = 3
ax = [150,300,-.2,3.5]
# =#
# #=
file = "matr5n/"
nums = [102,106,117,125,127,130,181,248]
σ = 4.02
k = 3
ax = [0,100,-.2,3]
# =#

figure("fig1", figsize=(4,6))

# =====================================================
subplot(2,1,1)

for n in nums
	local x = readdlm("data-pj/fig1/"*file*"fig/fort.$n")
	PyPlot.plot(x[:,1],x[:,2])
end

axis(ax)
xticks(fontname="serif")
yticks(fontname="serif")
xlabel("t [a.u.]",fontname="serif")
ylabel("N(t)",fontname="serif")


# =====================================================
subplot(2,1,2)
zer0 = 1e-12
Si = 157
x = readdlm("data-pj/fig1/"*file*"matr.dat")
A = zeros(Si,Si)
for l in 1:size(x)[1]
	i = Int64(x[l,1])
	j = Int64(x[l,2])
	A[i,j] = x[l,3]
end
y = readdlm("data-pj/fig1/"*file*"survivingspecies_abundancies.dat")
surv = Int64.(y[k:end,1])[y[k:end,2] .> zer0]
A = A[surv,surv]
S = length(surv)
pα = 10

κ = 1.
μ = 5.
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
xticks(fontname="serif")
yticks(fontname="serif")
xlabel("Re(ϵ)",fontname="serif")
ylabel("Im(ϵ)",fontname="serif")

