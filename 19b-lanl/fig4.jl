using PyPlot, DelimitedFiles, FFTW, LinearAlgebra

τ = .1
N = 3000
T = N*τ
n1 = 19
n2 = 20

ls1 = 1:n1
ls2 = 1:n2
ks1 = 5:205
ks2 = 5:205

ff1 = 2.4
ff2 = ff1*.4

ntw = "ntw20"
show_graph = true

L1 = zeros(n1,length(ks1))
for i in ls1
	for j in 1:length(ks1)
		L1[i,j] = readdlm("data_melvyn/"*ntw*"/l0_missing_hidden_20_test__l0_$(i).$(ks1[j])_obj.csv",',')[1]
	end
end
nL1 = (L1[ls1,:] .- maximum(L1[ls1,:]))./(maximum(L1[ls1,:]) - minimum(L1[ls1,:]))

L2 = zeros(n2,length(ks2))
for i in ls2
	for j in 1:length(ks2)
		L2[i,j] = readdlm("data_melvyn/"*ntw*"/l0_2_hidden_20_test__l0_$(i).$(ks2[j])_obj.csv",',')[1]
	end
end
nL2 = (L2[ls2,:] .- maximum(L2[ls2,:]))./(maximum(L2[ls2,:]) - minimum(L2[ls2,:]))

cols = [(243,111,33)./255,(0,145,194)./255,(161,0,202)./255]
gra = (.8,.8,.8,1.)

#=
x1 = findmin(nL1)
sol1 = x1[2][1]
nL1red = nL1[[1:sol1-1;sol1+1:size(nL1)[1]],:]
nL1max = [maximum(nL1red[:,i]) for i in 1:size(nL1red)[2]]
nL1min = [minimum(nL1red[:,i]) for i in 1:size(nL1red)[2]]
=#
sol11 = 13
sol12 = 14
sol13 = 15
nL1red = nL1[[1:12;16:19],:]
nL1max = [maximum(nL1red[:,i]) for i in 1:size(nL1red)[2]]
nL1min = [minimum(nL1red[:,i]) for i in 1:size(nL1red)[2]]

#=
x2 = findmin(nL2)
sol2 = x2[2][1]
nL2red = nL2[[1:sol2-1;sol2+1:size(nL2)[1]],:]
nL2max = [maximum(nL2red[:,i]) for i in 1:size(nL2red)[2]]
nL2min = [minimum(nL2red[:,i]) for i in 1:size(nL2red)[2]]
=#
sol21 = 14
sol22 = 18
nL2red = nL2[[1:13;15:17;19:20],:]
nL2max = [maximum(nL2red[:,i]) for i in 1:size(nL2red)[2]]
nL2min = [minimum(nL2red[:,i]) for i in 1:size(nL2red)[2]]



figure("ntw20",(14,10))


subplot2grid((4,22),(0,0),rowspan=2,colspan=11)
PyPlot.plot([ff1,ff1],[-.1,1.1],"--",color="C7")
PyPlot.fill([1;ks1;ks1[end:-1:1];1]*2π/(N*τ),-[nL1max[1];nL1max;nL1min[end:-1:1];nL1min[1]],color=gra)
PyPlot.plot([1;ks1]*2π/(N*τ),-[nL1[sol11,1];nL1[sol11,:]],color=cols[1])
PyPlot.plot([1;ks1]*2π/(N*τ),-[nL1[sol12,1];nL1[sol12,:]],color=cols[2])
PyPlot.plot([1;ks1]*2π/(N*τ),-[nL1[sol13,1];nL1[sol13,:]],color=cols[3])
axis([ks1[1]*2π/(N*τ),ks1[end]*2π/(N*τ),-.1,1.1])
xlabel("freq")
ylabel("normalized \n log-likelihood")

subplot2grid((4,22),(0,11),rowspan=2,colspan=11)
PyPlot.plot([ff1,ff1],[-.1,1.1],"--",color="C7")
PyPlot.plot([ff2,ff2],[-.1,1.1],"--",color="C7")
PyPlot.fill([1;ks2;ks2[end:-1:1];1]*2π/(N*τ),-[nL2max[1];nL2max;nL2min[end:-1:1];nL2min[1]],color=gra)
PyPlot.plot([1;ks2]*2π/(N*τ),-[nL2[sol21,1];nL2[sol21,:]],color=cols[1])
PyPlot.plot([1;ks2]*2π/(N*τ),-[nL2[sol22,1];nL2[sol22,:]],color=cols[2])
axis([ks2[1]*2π/(N*τ),ks2[end]*2π/(N*τ),-.1,1.1])
xlabel("freq")
ylabel("normalized \n log-likelihood")







f = .3
t = LinRange(-1,11,200)

# Multi-sine

subplot2grid((4,22),(3,0),colspan=6,rowspan=1)
PyPlot.plot(t,sin.(2π*f*t) + sin.(4π*f*t) + sin.(6π*f*t),"k",linewidth=3.)
axis([0.,10.,-2.8,2.8])
xticks([])
yticks([])
#ylabel("forcing")
#xlabel("t")

# l0
ls = Array(1:20)
ks = Array(10:10:200)
δk = 10
L = zeros(20,20)
for i in 1:20
	for j in 1:20
		L[i,j] = readdlm("data/ntw20/ntw20_multisine_l0_$(ls[i]).$(ks[j])_obj.csv",',')[1]
	end
end
nL = (L .- maximum(L))./(maximum(L) - minimum(L))
nLmi,ii = findmin(nL)
iii = ii[1]
nLmax = [maximum(nL[[1:iii-1;iii+1:size(nL)[1]],i]) for i in 1:size(nL)[2]]
nLmin = [minimum(nL[[1:iii-1;iii+1:size(nL)[1]],i]) for i in 1:size(nL)[2]]
freq = round(ks[ii[2]]/100.,digits=3)

subplot2grid((4,22),(2,0),colspan=6,rowspan=1)
PyPlot.plot([f,f],[-.1,1.1],"--",color="C7")
PyPlot.fill([ks[1]-δk;ks;ks[end]+δk;ks[end]+δk;ks[end:-1:1];ks[1]-δk]/100.,-[nLmax[1];nLmax;nLmax[end];nLmin[end];nLmin[end:-1:1];nLmin[1]],color=gra)
#=
for i in 1:length(ks)
	k = ks[i]
	PyPlot.plot([1,1]*k/100.,[0.,-nL[iii,i]],color=cols[1])
end
PyPlot.plot(ks/100.,-nL[iii,:],"o",color=cols[1])
=#
PyPlot.plot(ks/100.,-nL[iii,:],color=cols[1])
xlabel("freq")
ylabel("normalized log-likelihood")
axis([(ks[1]-δk/2)/100.,(ks[end]+δk/2)/100.,-.1,1.1])


# Saw
subplot2grid((4,22),(3,6),colspan=6,rowspan=1)
PyPlot.plot(t,2*mod.(f*t,1.) .- 1,"k",linewidth=3.)
axis([0.,10.,-1.1,1.1])
xticks([])
yticks([])
#ylabel("forcing")
#xlabel("t")


# l0
ls = Array(1:20)
ks = Array(10:10:200)
L = zeros(20,20)
for i in 1:20
	for j in 1:20
		L[i,j] = readdlm("data/ntw20/ntw20_saw_l0_$(ls[i]).$(ks[j])_obj.csv",',')[1]
	end
end
nL = (L .- maximum(L))./(maximum(L) - minimum(L))
nLmi,ii = findmin(nL)
iii = ii[1]
nLmax = [maximum(nL[[1:iii-1;iii+1:size(nL)[1]],i]) for i in 1:size(nL)[2]]
nLmin = [minimum(nL[[1:iii-1;iii+1:size(nL)[1]],i]) for i in 1:size(nL)[2]]
freq = round(ks[ii[2]]/100.,digits=3)

subplot2grid((4,22),(2,6),colspan=6,rowspan=1)
PyPlot.plot([f,f],[-.1,1.1],"--",color="C7")
PyPlot.fill([ks[1]-δk;ks;ks[end]+δk;ks[end]+δk;ks[end:-1:1];ks[1]-δk]/100.,-[nLmax[1];nLmax;nLmax[end];nLmin[end];nLmin[end:-1:1];nLmin[1]],color=gra)
#=
for i in 1:length(ks)
	k = ks[i]
	PyPlot.plot([1,1]*k/100.,[0.,-nL[iii,i]],color=cols[1])
end
PyPlot.plot(ks/100.,-nL[iii,:],"o",color=cols[1])
=#
PyPlot.plot(ks/100.,-nL[iii,:],color=cols[1])
xlabel("freq")
ylabel("normalized log-likelihood")
axis([(ks[1]-δk/2)/100.,(ks[end]+δk/2)/100.,-.1,1.1])


# Step
subplot2grid((4,22),(3,12),colspan=6,rowspan=1)
PyPlot.plot(t,2*mod.(floor.(f*t*2),2).-1,"k",linewidth=3.)
axis([0.,10.,-1.1,1.1])
xticks([])
yticks([])
#ylabel("forcing")
#xlabel("t")

# l0
ls = Array(1:20)
ks = Array(10:10:200)
L = zeros(20,20)
for i in 1:20
	for j in 1:20
		L[i,j] = readdlm("data/ntw20/ntw20_step_l0_$(ls[i]).$(ks[j])_obj.csv",',')[1]
	end
end
nL = (L .- maximum(L))./(maximum(L) - minimum(L))
nLmi,ii = findmin(nL)
iii = ii[1]
nLmax = [maximum(nL[[1:iii-1;iii+1:size(nL)[1]],i]) for i in 1:size(nL)[2]]
nLmin = [minimum(nL[[1:iii-1;iii+1:size(nL)[1]],i]) for i in 1:size(nL)[2]]
freq = round(ks[ii[2]]/100.,digits=3)


subplot2grid((4,22),(2,12),colspan=6,rowspan=1)
PyPlot.plot([f,f],[-.1,1.1],"--",color="C7")
PyPlot.fill([ks[1]-δk;ks;ks[end]+δk;ks[end]+δk;ks[end:-1:1];ks[1]-δk]/100.,-[nLmax[1];nLmax;nLmax[end];nLmin[end];nLmin[end:-1:1];nLmin[1]],color=gra)
#=
for i in 1:length(ks)
	k = ks[i]
	PyPlot.plot([1,1]*k/100.,[0.,-nL[iii,i]],color=cols[1])
end
PyPlot.plot(ks/100.,-nL[iii,:],"o",color=cols[1])
=#
PyPlot.plot(ks/100.,-nL[iii,:],color=cols[1])
xlabel("freq")
ylabel("normalized log-likelihood")
axis([(ks[1]-δk/2)/100.,(ks[end]+δk/2)/100.,-.1,1.1])




if show_graph
	figure(ntw*": ntw")

	xy = readdlm("data_melvyn/"*ntw*"/"*ntw*"_xy.csv",',')
	adj = Int64.(readdlm("data_melvyn/"*ntw*"/"*ntw*"_adj.csv",','))
	
	for i in 1:2:size(adj)[1]
		PyPlot.plot(xy[adj[i,1:2],1],xy[adj[i,1:2],2],"k",linewidth=1.)
	end
	PyPlot.plot(xy[:,1],xy[:,2],"ok",markersize=5.)
	PyPlot.plot(xy[sol1,1],xy[sol1,2],"o",color=cols[1],markersize=8.)
	xticks([])
	yticks([])
end
