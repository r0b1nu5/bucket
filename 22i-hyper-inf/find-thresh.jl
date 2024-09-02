using Random, Dates, DelimitedFiles

include("hyper_inf.jl")
include("hyper_kuramoto.jl")
include("gen_rand_hyperg.jl")
include("tools_hyper.jl")

@info "############# START: $(now())"

n = 30; iter = 1000

 #=
ntw = "Simplicial-ER-py"
run = "001"
p1 = .001
p2 = .01
# =#
# #=
ntw = "Hyper-ER-py"
runs = ["003","004","005","006","007"]
p1 = .005
p2 = .05
# =#

cmapme = get_cmap("plasma")

for run in runs
A2 = zeros(n,n)
A3 = zeros(n,n,n)
A2l = zeros(0,3)
A3l = zeros(0,4)

if ntw in ["Simplicial-ER-py","Hyper-ER-py"]
	el = readdlm("data/edgelist-n$n-"*run*".csv",',')
	for l in 1:size(el)[1]
#		i,j,k,A2,A2l,A3,A3l
		i,j,k = el[l,:]
		if k == ""
			i,j = sort([i,j] .+ 1)
			A2[i,j] = A2[j,i] = 1.
			A2l = vcat(A2l,[i j 1.;j i 1.])
		else
			i,j,k = sort([i,j,k] .+ 1)
			A3[i,j,k] = A3[i,k,j] = A3[j,i,k] = A3[j,k,i] = A3[k,i,j] = A3[k,j,i] = 1.
			A3l = vcat(A3l,[i j k 1.;j i k 1.;k i j 1.])
		end
	end
end

adj = cat_As(A2,A3)
# ========================================================================

# #= ########## PERFECT MEASUREMENTS ###############
amplitude = 1.
ξ0 = 0.0005
X = amplitude*(rand(n,iter) .- .5)
Y = f_kuramoto_3rd(X,A2,A3,zeros(n),π/4,π/4) + ξ0*randn(size(X))
# =#



# ========================================================================

ooi = [2,3]
dmax = 2
λ = .01

xxx = hyper_inf(X,Y,ooi,dmax,λ)
A2us = xxx[1][2]
A3us = xxx[1][3]
coeff = xxx[2]
d = get_d(n,dmax)

@info "============= WE ARE DONE: $(now()) ================"

#th = LinRange(minimum([abs.(A2us[:,3]);abs.(A3us[:,4])]) - λ/10,maximum([abs.(A2us[:,3]);abs.(A3us[:,4])]) + λ/10,10)
th = LinRange(0,maximum(abs.(coeff)) + λ/10,100)

θ = get_θ(X,dmax)

auc = Float64[]
relerr = Float64[]

for t in th
	a2 = A2us[(abs.(A2us[:,3]) .> t),:]
	a3 = A3us[(abs.(A3us[:,4]) .> t),:]

	tpr,fpr = my_ROC(abs.(a2),A2l,abs.(a3),A3l,n)
	push!(auc,get_auc(tpr,fpr))

	push!(relerr,get_re(θ,Y,coeff.*(abs.(coeff) .> t)))
end

figure("Perf. vs thr.",(10,5))

subplot(2,1,1)
PyPlot.plot(th,auc)
ylabel("AUC")

subplot(2,1,2)
PyPlot.plot(th,relerr)
xlabel("Threshold")
ylabel("Relative error")
end



