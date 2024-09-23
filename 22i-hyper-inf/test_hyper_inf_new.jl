using Random, Dates, DelimitedFiles

include("hyper_inf.jl")
include("hyper_kuramoto.jl")
include("hyper_ktanh.jl")
include("gen_rand_hyperg.jl")
include("tools_hyper.jl")

include("arni-reconstruct.jl")
include("arni-reconstruct-3rd.jl")

@info "############# START: $(now())"

# Generating the hypergraph.
#n = 7; T = 150; iters = 10:10:150 		# Takes < 1sec
#n = 30; T = 2500; iters = 500:500:2500 	# Takes ~ 10sec
#n = 60; T = 2500; iters = 1000:500:2500	# Takes ~ 10min
n = 100; T = 2000; iters = 500:500:2000		

save = true
A2 = zeros(n,n)
A3 = zeros(n,n,n)
A2l = zeros(0,3)
A3l = zeros(0,4)


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
 #=
ntw = "Hyper-ER"
p1 = n > 30 ? .01 : .05
p2 = .4
# =#
 #= 
ntw = "Hyper-ER"
p1 = .5
p2 = .8
# =#
# #=
ntw = "Simplicial-ER-py"
run = "001"
p1 = .001
p2 = .01
# =#
 #=
ntw = "Hyper-ER-py"
run = "002"
p1 = .001
p2 = .01
# =#
 #=
ntw = "Simplicial-ER-corr"
p1 = n > 30 ? .01 : .05
p2 = .4
# =#

if ntw in ["Wheel", "Hyper-wheel"]
	A2,A3 = gen_rand_hyperwheel(n,p1,p2,p3)
	A4 = zeros(n,n,n,n)
elseif ntw in ["ER", "Hyper-ER"]
	A2,A3,A2l,A3l = gen_hyper_er(n,p1,p2)
	A4 = zeros(n,n,n,n)
elseif ntw in ["Simplicial-ER","Simplicial-ER-2","Simplicial-ER-corr"]
	A2,A3,A2l,A3l = gen_simplicial_er(n,p1,p2)
	A4 = zeros(n,n,n,n)
elseif ntw in ["Simplicial-ER-py","Hyper-ER-py"]
	el = readdlm("data/edgelist-n$n-"*run*".csv",',')
	for l in 1:size(el)[1]
		global i,j,k,A2,A2l,A3,A3l
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


cmapme = get_cmap("RdPu")
cmaparni = get_cmap("GnBu")
cmapdiff = get_cmap("Greys")
cmapme = get_cmap("plasma")
cmaparni = get_cmap("plasma")
#cmapdiff = get_cmap("plasma")
c0 = 0.
c1 = 1.

#adj = get_adj_3rd(A2,A3)[1]
adj = cat_As(A2,A3)
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
amplitude = 1.
ξ0 = 0.0005
X = amplitude*(rand(n,T) .- .5)
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
	@info "Size of X is $(size(X)[2])/$T."
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
#iters = 300:10:400
#iters = 5000:1000:T
#iters = 10:10:200
#iters = 10:10:150  # Good for n = 7
#iters = 200:5:200
#iters = 1000:1000:5000
#iters = 500:500:2500
#iters = [T,]
ooi = [2,3]
dmax = 2
c = 0
for iter in iters
	global c += 1
	@info "Run $c/$(length(iters))"

	xxx = hyper_inf(X[:,1:iter],Y[:,1:iter],ooi,dmax,1e-1)
#	A2us = inferred_adj_2nd(xxx[1][2],n)[1]
	A2us = xxx[1][2]
#	A3us = inferred_adj_3rd(xxx[1][3],n)[1]
	A3us = xxx[1][3]

	if save
		writedlm("data/kuramoto-"*ntw*"-n$n-iter$iter-A2.csv",A2,',')
		writedlm("data/kuramoto-"*ntw*"-n$n-iter$iter-A3.csv",A3,',')
		writedlm("data/kuramoto-"*ntw*"-n$n-iter$iter-A2this.csv",A2us,',')
		writedlm("data/kuramoto-"*ntw*"-n$n-iter$iter-A3this.csv",A3us,',')
	end
 #=
	A2us_ = zeros(n,n)
	for l in 1:size(A2us)[1]
		i = Int64(A2us[l,1])
		j = Int64(A2us[l,2])
		A2us_[i,j] = A2us[l,3]
	end
# =#
 #=
	A3us_ = zeros(n,n,n)
	for l in 1:size(A3us)[1]
		i = Int64(A3us[l,1])
		j = Int64(A3us[l,2])
		k = Int64(A3us[l,3])
		A3us_[i,j,k] = A3us[l,4]
		A3us_[i,k,j] = A3us[l,4]
	end
	adjus_ = get_adj_3rd(A2us_,A3us_)[1]
	adju_ = cat_As(A2us_,A3us_)
# =#
 #=
	A2us_ = inferred_adj_2nd(xxx[2][2],n)[2]
	A3us_ = inferred_adj_3rd(xxx[2][3],n)[2]
	adjus_ = get_adj_3rd(A2us_,A3us_)[1]
	adju_ = cat_As(A2us_,A3us_)
# =#
@info "============= WE ARE DONE: $(now()) ================"

	ξ = 1e-10

	tpr,fpr = my_ROC(abs.(A2us),A2l,abs.(A3us),A3l,n)
	tpr2,fpr2 = my_ROC(abs.(A2us),A2l,n)
	tpr3,fpr3 = my_ROC(abs.(A3us),A3l,n)

	figure("ROCs-"*ntw*"-$n",(13.5,4))
	subplot(1,3,1)
	PyPlot.plot(fpr,tpr,color=cmapme((iter-minimum(iters))/max(1,(maximum(iters)-minimum(iters)))))
	subplot(1,3,2)
	PyPlot.plot(fpr2,tpr2,color=cmapme((iter-minimum(iters))/max(1,(maximum(iters)-minimum(iters)))))
	subplot(1,3,3)
	PyPlot.plot(fpr3,tpr3,color=cmapme((iter-minimum(iters))/max(1,(maximum(iters)-minimum(iters)))))

 #=
	rocadjus = roc(abs.(adju_) + ξ*rand(Float64,size(adju_)),adj)
	subplot(1,3,1)
	PyPlot.plot(rocadjus.FPR,rocadjus.TPR,"k")
# =# 
 #=
	rocA2us = roc(abs.(A2us_) + ξ*rand(n,n),A2)
	subplot(1,3,2)
	PyPlot.plot(rocA2us.FPR,rocA2us.TPR,":k",alpha=1-.5*((iter-minimum(iters))/max(1,maximum(iters)-minimum(iters))))
# =#
#=
	rocA3us = roc(abs.(A3us_) + ξ*rand(n,n,n),A3)
	subplot(1,3,3)
	PyPlot.plot(rocA3us.FPR,rocA3us.TPR,":k",alpha=1-.5*((iter-minimum(iters))/max(1,maximum(iters)-minimum(iters))))
# =#
end



figure("ROCs-"*ntw*"-$n",)
subplot(1,3,1)
xlabel("FPR")
ylabel("TPR")
title("THIS")
subplot(1,3,2)
xlabel("FPR")
title("THIS, pairwise")
subplot(1,3,3)
title("THIS, triadic")
xlabel("FPR")

@info "############# FINISHED: $(now())"

