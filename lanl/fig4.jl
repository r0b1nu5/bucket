using PyPlot, DelimitedFiles, RecipesBase, Shapefile, Colors

include("processing_tools.jl")

cmap = get_cmap("plasma")
colshift1 = .5
colshift2 = .5
cols = [cmap(1-(i+colshift1)/(2+colshift1+colshift2)) for i in 0:2]

ntw = "utk"
n = 99
ls = 1:n
ks = 1:1500
T = 300.
file = "mysterious_forcing"
fs1 = 85

L0 = zeros(length(ls),length(ks))
for i in 1:length(ls)
	@info "i = $i"
	for j in 1:length(ks)
		L0[i,j] = readdlm("data/utk/"*file*"_l0_$(ls[i]).$(ks[j])_obj.csv",',')[1]
	end
end
nL0 = (L0 .- maximum(L0))./(maximum(L0) - minimum(L0))

@info "L0 is loaded."

##############################################################

figure("UTK",(14,10))

subplot2grid((4,6),(0,0),colspan=3,rowspan=2)

# #=
s_xy, red_i, name, l_xy, lids, lakes = utk_preprocess_boundaries()

for i in 1:length(s_xy)
	s = s_xy[i]
	for j in 1:length(red_i[i])
		PyPlot.fill(s[1][red_i[i][j]],s[2][red_i[i][j]],color=(.9,.9,.9,1.))
		PyPlot.plot(s[1][red_i[i][j]],s[2][red_i[i][j]],"-k",linewidth=.5)
	end
end

for i in 1:length(l_xy)
	l = l_xy[i]
	ids = lids[i]
	PyPlot.fill(l[1][ids],l[2][ids],color=(1.,1.,1.,1.))
	PyPlot.plot(l[1][ids],l[2][ids],"-k",linewidth=.5)
end

# =#

xy = readdlm("data_utk/utk1_coord.csv",',')
ids = (1:size(xy)[1]).*(xy[:,1] .< -55.)
ids = setdiff(ids,[0,])
PyPlot.plot(xy[ids,1],xy[ids,2],"ok",markersize=5.)
PyPlot.plot(xy[fs1,1],xy[fs1,2],"o",color=cols[2],markersize=10.)
axis([-105.,-63.,24.,48.])
PyPlot.xticks([])
PyPlot.yticks([])


subplot2grid((4,6),(0,3),colspan=3,rowspan=2)

L0red = nL0[[1:fs1-1;fs1+1:size(nL0)[1]],:]
L0max = [maximum(L0red[:,i]) for i in 1:size(L0red)[2]]
L0min = [minimum(L0red[:,i]) for i in 1:size(L0red)[2]]
PyPlot.fill([ks;ks[end:-1:1]]/T,-[L0max;L0min[end:-1:1]],color="C7")
PyPlot.plot(ks/T,-nL0[fs1,:],"-",color=cols[2],linewidth=2.,label="$fs1")
axis([0.,maximum(ks)/T,-.1,1.1])
xlabel("freq")
ylabel("normalized log-likelihodd")
legend()

########################## FORCING TYPES #######################

f = .3
t = LinRange(-1,11,200)

# Multi-sine

subplot2grid((4,6),(3,0),colspan=2,rowspan=1)
PyPlot.plot(t,sin.(2π*f*t) + sin.(4π*f*t) + sin.(6π*f*t),"k",linewidth=3.)
axis([0.,10.,-2.8,2.8])
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
		L[i,j] = readdlm("data/ntw20/ntw20_multisine_l0_$(ls[i]).$(ks[j])_obj.csv",',')[1]
	end
end
nL = (L .- maximum(L))./(maximum(L) - minimum(L))
nLmi,ii = findmin(nL)
iii = ii[1]
nLmax = [maximum(nL[[1:iii-1;iii+1:size(nL)[1]],i]) for i in 1:size(nL)[2]]
nLmin = [minimum(nL[[1:iii-1;iii+1:size(nL)[1]],i]) for i in 1:size(nL)[2]]
freq = round(ks[ii[2]]/100.,digits=3)

subplot2grid((4,6),(2,0),colspan=2,rowspan=1)
PyPlot.fill([ks;ks[end:-1:1]]/100.,-[nLmax;nLmin[end:-1:1]],color="C7")
for i in 1:length(ks)
	k = ks[i]
	PyPlot.plot([1,1]*k/100.,[0.,-nL[iii,i]],color=cols[2])
end
PyPlot.plot(ks/100.,-nL[iii,:],"o",color=cols[2])
xlabel("freq")
ylabel("normalized log-likelihood")
axis([ks[1]/100.,2.,-.1,1.1])


# Saw
subplot2grid((4,6),(3,2),colspan=2,rowspan=1)
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

subplot2grid((4,6),(2,2),colspan=2,rowspan=1)
PyPlot.fill([ks;ks[end:-1:1]]/100.,-[nLmax;nLmin[end:-1:1]],color="C7")
for i in 1:length(ks)
	k = ks[i]
	PyPlot.plot([1,1]*k/100.,[0.,-nL[iii,i]],color=cols[2])
end
PyPlot.plot(ks/100.,-nL[iii,:],"o",color=cols[2])
xlabel("freq")
ylabel("normalized log-likelihood")
axis([ks[1]/100.,2.,-.1,1.1])


# Step
subplot2grid((4,6),(3,4),colspan=2,rowspan=1)
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


subplot2grid((4,6),(2,4),colspan=2,rowspan=1)
PyPlot.fill([ks;ks[end:-1:1]]/100.,-[nLmax;nLmin[end:-1:1]],color="C7")
for i in 1:length(ks)
	k = ks[i]
	PyPlot.plot([1,1]*k/100.,[0.,-nL[iii,i]],color=cols[2])
end
PyPlot.plot(ks/100.,-nL[iii,:],"o",color=cols[2])
xlabel("freq")
ylabel("normalized log-likelihood")
axis([ks[1]/100.,2.,-.1,1.1])






