using PyPlot, DelimitedFiles

 #=
cmap = get_cmap("plasma")
colshift1 = .5
colshift2 = .5
cols = [cmap(1-(i+colshift1)/(2+colshift1+colshift2)) for i in 0:2]
# =#
# #=
cols = [(243,111,33)./255,(0,145,194)./255,(161,0,202)./255]
# =#
gra = (.8,.8,.8,1.)

ntw = "ntw20"
ex = 1
Xs = readdlm("data_melvyn/"*ntw*"/"*ntw*"_ex$(ex)_Xs.csv",',')
nn,N = size(Xs)
n = Int64(nn/2)
ls = 1:n
ks = 1:40
δk = 1.
τ = .1
T = (N-1)*τ

ff = .08
fs = 20

L0 = zeros(n,length(ks))
for i in ls
	for j in 1:length(ks)
		L0[i,j] = readdlm("data_melvyn/"*ntw*"/"*ntw*"_ex$(ex)_l0_$(i).$(ks[j])_obj.csv",',')[1]
	end
end
nL0 = (L0[ls,:] .- maximum(L0[ls,:]))./(maximum(L0[ls,:]) - minimum(L0[ls,:]))

xxx = findmin(nL0)
l0 = xxx[2][1]
L0red = nL0[[1:l0-1;l0+1:size(nL0)[1]],:]
L0max = [maximum(L0red[:,i]) for i in 1:size(L0red)[2]]
L0min = [minimum(L0red[:,i]) for i in 1:size(L0red)[2]]

L1 = Vector{Float64}()
γ1 = Vector{Vector{Float64}}()
for j in 1:length(ks)
	push!(L1,readdlm("data_melvyn/"*ntw*"/"*ntw*"_ex$(ex)_l1_$(ks[j])_obj.csv",',')[1])
	push!(γ1,vec(readdlm("data_melvyn/"*ntw*"/"*ntw*"_ex$(ex)_l1_$(ks[j])_gamma.csv",',')))
end
nL1 = (L1 .- maximum(L1))./(maximum(L1) - minimum(L1))

xxx = findmin(nL1)
k1 = xxx[2]
xxx = findmax(γ1[k1])
l1 = xxx[2]
nγ1 = γ1[k1] ./ maximum(γ1[k1])

xy = readdlm("data_melvyn/"*ntw*"/"*ntw*"_xy.csv",',')
A = readdlm("data_melvyn/ntw20/ntw20_A.csv",',')


figure("ntw20",(14,4))

subplot2grid((20,3),(0,0),rowspan=20,colspan=1)

PyPlot.plot([ff,ff],[-.1,1.1],"--",color="C7")
PyPlot.fill([ks[1]-δk;ks;ks[end]+δk;ks[end]+δk;ks[end:-1:1];ks[1]-δk]./T,-[L0max[1];L0max;L0max[end];L0min[end];L0min[end:-1:1];L0min[1]],color=gra)
#for i in 1:length(ks)
#	PyPlot.plot([1,1].*ks[i]/T,[0.,-nL0[l0,i]],color=cols[1])
#end
PyPlot.plot(ks/T,-nL0[l0,:],"-",color=cols[1])
axis([(ks[1]-δk/2)/T,(ks[end]+δk/2)/T,-.1,1.1])
xlabel("freq")
ylabel("rescaled log-likelihood")

subplot2grid((20,3),(0,1),rowspan=10,colspan=1)

PyPlot.plot([ff,ff],[-.1,1.1],"--",color="C7")
#for j in 1:length(ks)
#	PyPlot.plot([ks[j],ks[j]]./T,[0.,-nL1[j]],color=cols[1])
#end
PyPlot.plot(ks./T,-nL1,"-",color=cols[1])
axis([(ks[1]-δk/2)/T,(ks[end]+δk/2)/T,-.1,1.1])
xlabel("freq")
ylabel("rescaled \n log-likelihood")

subplot2grid((20,3),(13,1),rowspan=7,colspan=1)

PyPlot.plot([fs,fs],[-.1,1.1],"--",color="C7")
PyPlot.plot(ls,nγ1,"ok",markersize=4.)
PyPlot.plot(l1,nγ1[l1],"o",color=cols[1],markersize=6.)
axis([0,n+1,-.1,1.1])
xlabel("node index")
ylabel("amplitude")

subplot2grid((20,3),(5,2),colspan=1,rowspan=15)

for i in 1:n-1
	for j in i+1:n
		if A[i,j] > 1e-2
			PyPlot.plot(xy[[i,j],1],-xy[[i,j],2],"k")
		end
	end
end
PyPlot.plot(xy[:,1],-xy[:,2],"ok",markersize=6.)
PyPlot.plot(xy[l0,1],-xy[l0,2],"o",color=cols[1],markersize=8.)
xticks([])
yticks([])




