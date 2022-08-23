using PyPlot, DelimitedFiles, FFTW, LinearAlgebra

τ = .1
N = 3000
T = N*τ
n1 = 19
n2 = 20

ls1 = 1:n1
ls2 = 1:n2
ks1 = [5:5:100;110:120;125:5:205]
ks2 = [5:5:35;40:60;110:120;125:5:205]

ff1 = 2.4
ff2 = ff1*.4

ntw = "ntw20"
show_graph = true

L1 = zeros(n1,length(ks1))
for i in ls1
	for j in 1:length(ks1)
		L1[i,j] = readdlm("data_melvyn/"*ntw*"/missing_hidden_20_test__l0_$(i).$(ks1[j])_obj.csv",',')[1]
	end
end
nL1 = (L1[ls1,:] .- maximum(L1[ls1,:]))./(maximum(L1[ls1,:]) - minimum(L1[ls1,:]))

L2 = zeros(n2,length(ks2))
for i in ls2
	for j in 1:length(ks2)
		L2[i,j] = readdlm("data_melvyn/"*ntw*"/2_hidden_20_test__l0_$(i).$(ks2[j])_obj.csv",',')[1]
	end
end
nL2 = (L2[ls2,:] .- maximum(L2[ls2,:]))./(maximum(L2[ls2,:]) - minimum(L2[ls2,:]))

cols = [(243,111,33)./255,(0,145,194)./255,(161,0,202)./255]
gra = (.8,.8,.8,1.)

x1 = findmin(nL1)
sol1 = x1[2][1]
nL1red = nL1[[1:sol1-1;sol1+1:size(nL1)[1]],:]
nL1max = [maximum(nL1red[:,i]) for i in 1:size(nL1red)[2]]
nL1min = [minimum(nL1red[:,i]) for i in 1:size(nL1red)[2]]

x2 = findmin(nL2)
sol2 = x2[2][1]
nL2red = nL2[[1:sol2-1;sol2+1:size(nL2)[1]],:]
nL2max = [maximum(nL2red[:,i]) for i in 1:size(nL2red)[2]]
nL2min = [minimum(nL2red[:,i]) for i in 1:size(nL2red)[2]]





figure("ntw20",(10,5))


subplot(2,1,1)
PyPlot.plot([ff1,ff1],[-.1,1.1],"--",color="C7")
PyPlot.fill([1;ks1;ks1[end:-1:1];1]/(N*τ),-[nL1max[1];nL1max;nL1min[end:-1:1];nL1min[1]],color=gra)
PyPlot.plot([1;ks1]/(N*τ),-[nL1[sol1,1];nL1[sol1,:]],color=cols[1])
axis([ks1[1]/(N*τ),ks1[end]/(N*τ),-.1,1.1])
xlabel("freq")
ylabel("normalized \n log-likelihood")

subplot(2,1,2)
PyPlot.plot([ff1,ff1],[-.1,1.1],"--",color="C7")
PyPlot.fill([1;ks2;ks2[end:-1:1];1]/(N*τ),-[nL2max[1];nL2max;nL2min[end:-1:1];nL2min[1]],color=gra)
PyPlot.plot([1;ks2]/(N*τ),-[nL2[sol2,1];nL2[sol2,:]],color=cols[1])
axis([ks2[1]/(N*τ),ks2[end]/(N*τ),-.1,1.1])
xlabel("freq")
ylabel("normalized \n log-likelihood")



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







