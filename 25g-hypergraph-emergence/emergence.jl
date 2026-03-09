using Random, Dates, DelimitedFiles

include("kuramoto.jl")
include("this.jl")
include("graph-tools.jl")
include("coarse-grain.jl")

@info "############# START: $(now())"

# Generating the hypergraph.
n = 8; T = 150; iters = 10:10:150 		# Takes < 1sec
#n = 30; T = 2500; iters = 500:500:2500 	# Takes ~ 10sec
#n = 60; T = 2500; iters = 1000:500:2500	# Takes ~ 10min
#n = 100; T = 2000; iters = 500:500:2000		

save = true # Saving the inference and the ground truth?

p = .4
A,B = gen_rand_graph(n,p)

plot_graph(A)
# ========================================================================

amplitude = 1.
ξ0 = 0.0005
X = amplitude*(rand(n,T) .- .5)
Y = f_kuramoto(X,zeros(n),B,1.,π/4)


#################### INFERENCE ###################################

ooi = [2,3]
dmax = 2

Ainf,coeff,relerr = this(X,Y,ooi,dmax)
figure()
plot_hypergraph(A,Ainf)

# Coarse graining
k = 4
g,X2,Y2 = coarse_grain(A,k,X,Y)
Ainf2,coeff2,relerr2 = this(X2,Y2,ooi,dmax)
figure()
plot_hypergraph(A,Ainf2)



 #=
c = 0
for iter in iters
	global c += 1
	@info "Run $c/$(length(iters))"

	# Run THIS
	xxx = hyper_inf(X[:,1:iter],Y[:,1:iter],ooi,dmax,1e-1)
	A2us = xxx[1][2] # Pairwise interactions
	A3us = xxx[1][3] # Triadic interactions

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

	cmapme = get_cmap("RdPu")
	
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
# =#
