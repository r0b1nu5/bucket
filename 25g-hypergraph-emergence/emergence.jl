using Random, Dates, DelimitedFiles

include("kuramoto.jl")
include("this.jl")
include("graph-tools.jl")
include("coarse-grain.jl")

@info "############# START: $(now())"

# Generating the hypergraph.
n = 30; T = 2500; iters = 10:10:150 		# Takes < 1sec
#n = 30; T = 2500; iters = 500:500:2500 	# Takes ~ 10sec
#n = 60; T = 2500; iters = 1000:500:2500	# Takes ~ 10min
#n = 100; T = 2000; iters = 500:500:2000		

save = true # Saving the inference and the ground truth?

p = .2
A,B = gen_rand_graph(n,p)
m = Int64(nnz(A)/2)

plot_graph(A)
# ========================================================================

amplitude = 1.
ξ0 = 0.0005
X = amplitude*(rand(n,T) .- .5)
Y = f_kuramoto(X,zeros(n),B,1.,π/4)


#################### INFERENCE ###################################

ooi = [2,3]
dmax = 2
zer0 = 1e-4

Ainf,coeff,relerr = this(X,Y,ooi,dmax)
figure()
plot_hypergraph(A,Ainf)
title("k = 0")

 #=
# Coarse graining, one shot
k = 3
g,X2,Y2,A2 = coarse_grain(A,k,X,Y)
Ainf2,coeff2,relerr2 = this(X2,Y2,ooi,dmax)
figure()
plot_hypergraph(A,Ainf2)
# =#

# #=
# Coarse graining, step by step
n,m = size(B)
ids = randperm(m)
ks = 1:1:15

m2max = 1.
m3max = 1.
#m2max = (n*(n-1))
#m3max = (n*(n-1)*(n-2))
m2 = [sum(abs.(Ainf[2][:,3]) .> zer0)/m2max,]
m3 = [sum(abs.(Ainf[3][:,4]) .> zer0)/m3max,]
 
# Contribution of 2-edges to the dynamics
mag2 = zeros(size(X)[2])
for i in 1:size(Ainf[2])[1]
	a = Int64(Ainf[2][i,2])
	global mag2 += abs.(Ainf[2][i,3]*X[a,:])
end
contribution2 = [median(mag2),]
# Contribution of 3-edges to the dynamics
mag3 = zeros(size(X)[2])
for i in 1:size(Ainf[3])[1]
	a = Int64(Ainf[3][i,2])
	b = Int64(Ainf[3][i,3])
	global mag3 += abs.(Ainf[3][i,4]*(X[a,:].*X[b,:]))
end
contribution3 = [median(mag3),]

distances = Vector{Int64}[]

for k in ks
	@info "k/kmax = $k/$(maximum(ks))"

	clusters = [[i,] for i in 1:n]

	id = sort(ids[1:k])
	for i in 1:n
		for e in id
			if abs(B[i,e]) > .1
				j = findmax(abs.([B[1:i-1,e];0;B[i+1:n,e]]))[2]
				push!(clusters[i],j)
			end
		end
	end

	final_clusters = Vector{Int64}[]
	for i in n:-1:1
		x = clusters[i]
		j = i
		test = true
		while test && j > 1
			j -= 1
			if length(intersect(x,clusters[j])) > 0
				clusters[j] = union(x,clusters[j])
				test = false
			end
		end
		if test
			push!(final_clusters,x)
		end
	end

	X2 = zeros(0,T)
	Y2 = zeros(0,T)
	nc = length(final_clusters)

	for c in final_clusters
		X2 = [X2;sum(X[c,:],dims=1)./length(c)]
		Y2 = [Y2;sum(Y[c,:],dims=1)./length(c)]
	end

	Ainf2,coeff2,relerr2 = this(X2,Y2,ooi,dmax)

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
	
	# Tracking where 3-edges come from
	D = fill(Inf,nc,nc)
	D[nc,nc] = 0.
	for i in 1:nc-1
		D[i,i] = 0.
		for j in i+1:nc
			if maximum(A[final_clusters[i],final_clusters[j]]) > zer0
				D[i,j] = 1.
				D[j,i] = 1.
			end
		end
	end

	# Floyd-Warshall algorithm (gets pairwise distances)
	for k in 1:nc
		for i in 1:nc
			for j in 1:nc
				if D[i,k] + D[k,j] < D[i,j]
					D[i,j] = D[i,k] + D[k,j]
				end
			end
		end
	end

	cluster_id = (1:nc)[length.(final_clusters) .> 1]
	dist = Int64[]
	for l in 1:size(Ainf2[3])[1]
		ids = Int64.(Ainf2[3][l,1:3])
		push!(dist,sum(minimum(D[ids,cluster_id],dims=2)))
	end
	push!(distances,dist)

 #=
	figure()
	plot_hypergraph(A,Ainf2)
	title("k = $k")
# =#
end


# Plot number of edges (normalized)
fig1, ax11 = subplots()

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

tight_layout()
show()

# Plot contribution of edges
fig2, ax21 = subplots()

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

# Plot distance from inferred 3-edges to clusters of nodes
figure()
PyPlot.bar(ks,ones(length(ks)),color="C4")
PyPlot.bar(ks,[sum(distances[i] .<= 3) for i in 1:length(distances)]./(m3[2:end]),color="C3")
PyPlot.bar(ks,[sum(distances[i] .<= 2) for i in 1:length(distances)]./(m3[2:end]),color="C2")
PyPlot.bar(ks,[sum(distances[i] .<= 1) for i in 1:length(distances)]./(m3[2:end]),color="C1")
PyPlot.bar(ks,[sum(distances[i] .<= 0) for i in 1:length(distances)]./(m3[2:end]),color="C0")
xlabel("k")

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
