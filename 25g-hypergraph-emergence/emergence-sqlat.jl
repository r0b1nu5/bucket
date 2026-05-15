using Random, Dates, DelimitedFiles

include("kuramoto.jl")
include("this.jl")
include("graph-tools.jl")
include("coarse-grain.jl")
include("gen-lattice.jl")

@info "############# START: $(now())"

n = 16
ks = [1,2]

save = true # Saving the inference and the ground truth?

A,B = gen_square_lattice(n)
m = Int64(nnz(A)/2)

# ========================================================================

amplitude = 1.
ξ0 = 0.0005
X = amplitude*(rand(n^2,T) .- .5)
Y = f_kuramoto(X,zeros(n^2),B,1.,π/4)

#################### INFERENCE ###################################

ooi = [2,3]
dmax = 2
zer0 = 1e-4
λ = 0.1

Ainf,coeff,relerr = this(X,Y,ooi,dmax,λ)

m2max = 1.
m3max = 1.
#m2max = (n*(n-1))
#m3max = (n*(n-1)*(n-2))
m2 = [sum(abs.(Ainf[2][:,3]) .> zer0)/m2max,]
m3 = [sum(abs.(Ainf[3][:,4]) .> zer0)/m3max,]
 
# Contribution of 2-edges to the dynamics
global mag2 = zeros(size(X)[2])
for i in 1:size(Ainf[2])[1]
	a = Int64(Ainf[2][i,2])
	global mag2 += abs.(Ainf[2][i,3]*X[a,:])
end
contribution2 = [median(mag2),]
# Contribution of 3-edges to the dynamics
global mag3 = zeros(size(X)[2])
for i in 1:size(Ainf[3])[1]
	a = Int64(Ainf[3][i,2])
	b = Int64(Ainf[3][i,3])
	global mag3 += abs.(Ainf[3][i,4]*(X[a,:].*X[b,:]))
end
contribution3 = [median(mag3),]

distances = Vector{Int64}[]

for k in ks
	@info "k/kmax = $k/$(maximum(ks))"
	
	A2,B2,X2 = coarse_grain_sqlat_4(n,X)
	A2,B2,Y2 = coarse_grain_sqlat_4(n,Y)

	Ainf2,coeff2,relerr2 = this(X2,Y2,ooi,dmax,λ)

	push!(m2,sum(abs.(Ainf2[2][:,3]) .> zer0))
	push!(m3,sum(abs.(Ainf2[3][:,4]) .> zer0))
	#push!(m2,sum(abs.(Ainf2[2][:,3]) .> zer0)/(nc*(nc-1)))
	#push!(m3,sum(abs.(Ainf2[3][:,4]) .> zer0)/(nc*(nc-1)*(nc-2)))

	# Contribution of 2-edges to the dynamics
	mag2 = zeros(size(X2)[2])
	for i in 1:size(Ainf2[2])[1]
		a = Int64(Ainf2[2][i,2])
		mag2 += abs.(Ainf2[2][i,3]*X2[a,:])
	end
	push!(contribution2,median(mag2))
	# Contribution of 3-edges to the dynamics
	mag3 = zeros(size(X2)[2])
	for i in 1:size(Ainf2[3])[1]
		a = Int64(Ainf2[3][i,2])
		b = Int64(Ainf2[3][i,3])
		mag3 += abs.(Ainf2[3][i,4]*(X2[a,:].*X2[b,:]))
	end
	push!(contribution3,median(mag3))
	
end


# Plot number of edges (normalized)
fig1, (ax11,ax21) = subplots(1,2,figsize=(15,5))

ax11.plot([0;ks], m2, color="C0")
ax11.set_xlabel("k")
ax11.set_ylabel("#2-edges", color="C0")
ax11.set_ylim(-maximum(m2)*0.05,maximum(m2)*1.05)
ax11.tick_params(axis="y", labelcolor="C0")

ax12 = ax11.twinx()   # Share the same X-axis
ax12.plot([0;ks], m3, color="C1")
ax12.set_ylabel("#3-edges", color="C1")
ax12.set_ylim(-maximum(m3)*0.05,maximum(m3)*1.05)
ax12.tick_params(axis="y", labelcolor="C1")

title("λ = $λ")

tight_layout()
show()

# Plot contribution of edges
#fig2, ax21 = subplots(1,2,2)

ax21.plot([0;ks], contribution2, color="C0")
ax21.set_xlabel("k")
ax21.set_ylabel("2-edges contribution", color="C0")
ax21.set_ylim(-maximum(contribution2)*0.05,maximum(contribution2)*1.05)
ax21.tick_params(axis="y", labelcolor="C0")

ax22 = ax21.twinx()
ax22.plot([0;ks], contribution3, color="C1")
ax22.set_ylabel("3-edges contribution", color="C1")
ax22.set_ylim(-maximum(contribution3)*0.05,maximum(contribution3)*1.05)
ax22.tick_params(axis="y", labelcolor="C1")

tight_layout()
show()


# =#


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
#end



