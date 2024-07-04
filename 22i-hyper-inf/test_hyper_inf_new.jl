using Random

include("hyper_inf.jl")
include("hyper_kuramoto.jl")
include("hyper_ktanh.jl")
include("gen_rand_hyperg.jl")
include("tools_hyper.jl")

# Generating the hypergraph.
n = 10
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
p1 = .01
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
	A2,A3,A2l,A3l = gen_hyper_er(n,p1,p2,true)
	A4 = zeros(n,n,n,n)
end

cmapme = get_cmap("RdPu")
cmaparni = get_cmap("GnBu")
cmapdiff = get_cmap("Greys")
cmapme = get_cmap("plasma")
cmaparni = get_cmap("plasma")
#cmapdiff = get_cmap("plasma")
c0 = 0.
c1 = 1.

#adj = get_adj_3rd(A2,A3)[1]
#adj = cat_As(A2,A3)
# ========================================================================

# #=
# Testing the efficiency of the inference using the result of the vector field directly.

# Generate the data
 #=
amplitude = .1
X = amplitude*(rand(n,400) .- .5)
Y = f_kuramoto_3rd(X,A2,A3,zeros(n),π/4,π/4) + .01*randn(size(X))
# =#
# #= ########## PERFECT MEASUREMENTS ###############
amplitude = 2.
ξ0 = 0.0005
X = amplitude*(rand(n,400) .- .5)
Y = f_kuramoto_3rd(X,A2,A3,zeros(n),π/4,π/4) + ξ0*randn(size(X))
# =#
 #= ########### TRUNCATE TIME SERIES IN THE δ-BOX #############
T = 10000
δ0 = .1
δ1 = 2.
h = .001
amplitude = 2π
X = zeros(n,0)
Yh = zeros(n,0)
Y = zeros(n,0)
count = 0
while size(X)[2] < T
	global count += 1
	xxx = hyper_k(2*A2,A3,zeros(n),amplitude*(rand(n) .- .5),0.,0.,π/4,π/4,h,Int64(50/h))
	xx = mod.(xxx[1] .+ π,2π) .- π
	xx = xx - repeat(xx[:,end],1,size(xx)[2])
	t = [(δ0 < norm(xx[:,i]) < δ1) for i in 1:(size(xx)[2]-1)]
	idx = vec(1:size(xx)[2]-1)[t]
	yyy = (xxx[1][:,2:end] - xxx[1][:,1:end-1])./h
	global X = [X xx[:,idx]]
	global Yh = [Yh yyy[:,idx]]
	global Y = [Y xxx[2][:,idx]]
end
@info "Number of time series used: $count"
# =#

# =#
 #=
ω = 2*rand(n)
ω .-= mean(ω)
xxx = hyper_k(A2,A3,ω,zeros(n),1.)
X = .2*(rand(n,400) .- .5) + repeat(xxx[1][:,end],1,400)
Y = f_kuramoto_3rd(X,A2,A3,ω)
# =#
 #=
X = amplitude*(rand(n,400) .- .5)
Y = f_ktanh_3rd(X,A2,A3,zeros(n))
# =#
 #=
ω = 2*rand(n)
ω .-= mean(ω)
xxx = hyper_ktanh(A2,A3,ω,zeros(n),1.)
X = .2*(rand(n,400) .- .5) + repeat(xxx[1][:,end],1,400)
Y = f_ktanh_3rd(X,A2,A3,ω)
# =#


# Compute the sensitivity and specificity of the inference for various lengths of time series.
#iters = 10:5:200
#iters = 10:15:150
#iters = 5000:1000:T
#iters = 10:10:200
iters = 200:5:200
ooi = [2,3]
dmax = 2
c = 0
for iter in iters
	global c += 1
	@info "Run $c/$(length(iters))"

	xxx = hyper_inf(X[:,1:iter],Y[:,1:iter],ooi,dmax,1e-1)
#	A2us = inferred_adj_2nd(xxx[1][2],n)[1]
#	A2us = inferred_adj_2nd(xxx[1][2],n)[2]
	A2us = xxx[1][2]
#	A3us = inferred_adj_3rd(xxx[1][3],n)[1]
#	A3us = inferred_adj_3rd(xxx[1][3],n)[2]
	A3us = xxx[1][3]
#	adjus = get_adj_3rd(A2us,A3us)[1]
#	adju = cat_As(A2us,A3us)
@info "============= WE ARE DONE ================"

	ξ = 1e-10
		
	tpr2,fpr2 = my_ROC(A2us,A2l,n)
	tpr3,fpr3 = my_ROC(A3us,A3l,n)

	figure("ROCs-"*ntw*"-$n",(15,4))
	subplot(1,3,1)
#	PyPlot.plot(rocadjus.FPR,rocadjus.TPR,color=cmapme(c0+c1*iter/maximum(iters)))
	subplot(1,3,2)
	PyPlot.plot(fpr2,tpr2,color=cmapme(c0+c1*iter/maximum(iters)))
	subplot(1,3,3)
	PyPlot.plot(fpr3,tpr3,color=cmapme(c0+c1*iter/maximum(iters)))

end


figure("ROCs-"*ntw*"-$n",(15,10))
subplot(1,3,1)
xlabel("FPR")
ylabel("TPR")
title("ROC adj, THIS")
subplot(1,3,2)
xlabel("FPR")
title("ROC A2, THIS")
subplot(1,3,3)
title("ROC A3, THIS")
xlabel("FPR")
ylabel("TPR")


