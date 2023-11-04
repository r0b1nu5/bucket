using Random

include("hyper_inf.jl")
include("hyper_kuramoto.jl")
include("hyper_ktanh.jl")
include("gen_rand_hyperg.jl")
include("tools_hyper.jl")

# Generating the hypergraph.
n = 7
 #=
ntw = "Hyper-wheel"
p1 = .3
p2 = .3 
p3 = .3
# =#
 #=
ntw = "Wheel"
p1 = 0.
p2 = 0.05
p3 = .3
# =# 
 #=
ntw = "ER"
p1 = 0.
p2 = .99
# =#
# #=
ntw = "Hyper-ER"
p1 = .05
p2 = .4
# =#
 #= 
ntw = "Hyper-ER"
p1 = .5
p2 = .8
# =#

if ntw in ["Wheel", "Hyper-wheel"]
	A2,A3 = gen_rand_hyperwheel(n,p1,p2,p3,true)
	A4 = zeros(n,n,n,n)
elseif ntw in ["ER", "Hyper-ER"]
	A2,A3 = gen_hyper_er(n,p1,p2,true)
	A4 = zeros(n,n,n,n)
end

cmapme = get_cmap("RdPu")
cmaparni = get_cmap("GnBu")
cmapdiff = get_cmap("Greys")
cmapme = get_cmap("magma")

c0 = 0.
c1 = 1.

adj = cat_As(A2,A3)
# ========================================================================

T = 200
ΔT = 5
NT = 80
amplitude = 2.
ξ0 = 0.0005

λs = [1e-3,2e-3,5e-3,1e-2,2e-2,5e-2,1e-1]

ooi = [2,3]
c = 0
for λ in λs
	global c += 1
	@info "Run $c/$(length(λs))"

	X = amplitude*(rand(n,T+1) .- .5)
	Y = f_kuramoto_3rd(X,A2,A3,zeros(n),π/4,π/4) + ξ0*randn(size(X))

	xxx = hyper_inf_sparsity(X,Y,ooi,4,-1e-4,λ)
	A2us = inferred_adj_2nd(xxx[1][2],n)[2]
	A3us = inferred_adj_3rd(xxx[1][3],n)[2]
	adjus = get_adj_3rd(A2us,A3us)[1]
	adju = cat_As(A2us,A3us)
@info "============= WE ARE DONE ================"


	ξ = 1e-10
		
	rocadjus = roc(abs.(adju) + ξ*rand(Float64,size(adju)),adj)
	rocA2us = roc(abs.(A2us) + ξ*rand(n,n),A2)
	rocA3us = roc(abs.(A3us) + ξ*rand(n,n,n),A3)

	figure("ROCs-"*ntw*"-$n",(15,7))
	subplot(1,3,1)
	PyPlot.plot(rocadjus.FPR,rocadjus.TPR,color=cmapme(c/length(λs)))
	subplot(1,3,2)
	PyPlot.plot(rocA2us.FPR,rocA2us.TPR,color=cmapme(c/length(λs)))
	subplot(1,3,3)
	PyPlot.plot(rocA3us.FPR,rocA3us.TPR,color=cmapme(c/length(λs)))
end

figure("ROCs-"*ntw*"-$n",(15,7))
subplot(1,3,1)
ylabel("TPR")
xlabel("FPR")
title("ROC adj, us")
subplot(1,3,2)
xlabel("FPR")
title("ROC A2, us")
subplot(1,3,3)
xlabel("FPR")
title("ROC A3, us")
