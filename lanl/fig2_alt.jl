using PyPlot, DelimitedFiles

cmap = get_cmap("plasma")
colshift = .5
cols = [cmap(1-(i+colshift)/(2+colshift)) for i in 0:2]

# #=
ntw = "ieee57"
n = 57
ls = 1:n
ks0 = 1:50
ks1 = 1:50
T = 200000*2e-3
file = "mysterious_forcing_57"
fs = 2
ff = .01
# =#

L0 = zeros(length(ls),length(ks0))
for i in 1:length(ls)
	for j in 1:length(ks0)
		L0[i,j] = readdlm("data/"*file*"_l0_$(ls[i]).$(ks0[j])_obj.csv",',')[1]
	end
end
nL0 = (L0 .- maximum(L0))./(maximum(L0) - minimum(L0))

 #=
L1 = Array{Float64,1}()
γ1 = Array{Array{Float64,1},1}()

for j in 1:length(ks1)
	push!(L1,readdlm("data/"*file*"_l1_$(ks1[j])_obj.csv",',')[1])
	push!(γ1,vec(readdlm("data/"*file*"_l1_$(ks1[j])_g.csv",',')))
end
L1mi,iii = findmin(L1)
# =#

figure("fig 2 (new)",(19,4.5))

subplot(1,4,1)
xy = readdlm("data_melvyn/ieee57_xy.csv",',')
adj = Int.(readdlm("data_melvyn/ieee57_adj.csv",','))
for i in 1:2:size(adj)[1]
	PyPlot.plot(xy[adj[i,1:2],1],xy[adj[i,1:2],2],"k",linewidth=1.)
end
PyPlot.plot(xy[:,1],xy[:,2],"ok",markersize=5.)
PyPlot.plot(xy[fs,1],xy[fs,2],"o",markersize=10.,color=cols[2])
axis([-.1,1.1,-.1,1.1])
axis("off")

subplot(1,4,2)
L0red = [nL0[1:fs-1,:];nL0[fs+1:end,:]]
L0max = [maximum(L0red[:,i]) for i in 1:size(L0red)[2]]
L0min = [minimum(L0red[:,i]) for i in 1:size(L0red)[2]]
PyPlot.fill([ks0;ks0[end:-1:1]]/T,[L0max;L0min[end:-1:1]],color="C7",alpha=.7)
PyPlot.plot([ff,ff],[.1,-1.1],"--",color="C7")
PyPlot.plot(ks0/T,nL0[fs,:],"-",color=cols[2])
axis([0.,50/T,-1.1,.1])
xlabel("freq")
ylabel("normalized inverse log-likelihodd")


# #=
ntw = "uk"
n = 120
ls = 1:n
ks0 = 1:50
ks1 = 1:50
T = 50000*.01
file = "mysterious_forcing_UK"
fs = 2
ff = .01
# =#

L0 = zeros(length(ls),length(ks0))
for i in 1:length(ls)
	for j in 1:length(ks0)
		L0[i,j] = readdlm("data/"*file*"_l0_$(ls[i]).$(ks0[j])_obj.csv",',')[1]
	end
end
nL0 = (L0 .- maximum(L0))./(maximum(L0) - minimum(L0))

subplot(1,4,3)
xy = readdlm("data_melvyn/uk_xy.csv",',')
xyb = readdlm("data_melvyn/uk_bord.csv",',')
adj = Int.(readdlm("data_melvyn/uk_adj.csv",','))
xshift = minimum(xy[:,1]) - .2
PyPlot.plot(xyb[:,1].-xshift,xyb[:,2],"k",linewidth=.5)
for i in 1:2:size(adj)[1]
	PyPlot.plot(xy[adj[i,1:2],1].-xshift,xy[adj[i,1:2],2],"k",linewidth=1.)
end
PyPlot.plot(xy[:,1].-xshift,xy[:,2],"ok",markersize=5.)
PyPlot.plot(xy[fs,1].-xshift,xy[fs,2],"o",markersize=8.,color=cols[1])
axis([0.,1.,0.,1.])
axis("off")

subplot(1,4,4)
L0red = [nL0[1:fs-1,:];nL0[fs+1:end,:]]
L0max = [maximum(L0red[:,i]) for i in 1:size(L0red)[2]]
L0min = [minimum(L0red[:,i]) for i in 1:size(L0red)[2]]
PyPlot.fill([ks0;ks0[end:-1:1]]/T,[L0max;L0min[end:-1:1]],color="C7",alpha=.7)
PyPlot.plot([ff,ff],[.1,-1.1],"--",color="C7")
PyPlot.plot(ks0/T,nL0[fs,:],"-",color=cols[1])
axis([0.,50/T,-1.1,.1])
xlabel("freq")
ylabel("normalized inverse log-likelihodd")


#=
L1 = Array{Float64,1}()
γ1 = Array{Array{Float64,1},1}()

for j in 1:length(ks1)
	push!(L1,readdlm("data/"*file*"_l1_$(ks1[j])_obj.csv",',')[1])
	push!(γ1,vec(readdlm("data/"*file*"_l1_$(ks1[j])_g.csv",',')))
end
L1mi,iii = findmin(L1)
=#



