using PyPlot, DelimitedFiles, RecipesBase, Shapefile, Colors

include("processing_tools.jl")

cmap = get_cmap("plasma")
colshift1 = .5
colshift2 = .5
cols = [cmap(1-(i+colshift1)/(2+colshift1+colshift2)) for i in 0:2]

# #=
ntw = "utk"
n = 99
ls = 1:n
ks = 1:1500
T = 300.
file = "mysterious_forcing"
fs1 = 51
fs2 = 66
fs3 = 85
# =#

L0 = zeros(length(ls),length(ks))
for i in 1:length(ls)
	for j in 1:length(ks)
		L0[i,j] = readdlm("data/"*file*"_l0_$(ls[i]).$(ks[j])_obj.csv",',')[1]
	end
end
nL0 = (L0 .- maximum(L0))./(maximum(L0) - minimum(L0))


##############################################################

figure("UTK (bis)",(14,4.5))

subplot(1,2,1)

s_xy, red_i, name, l_xy, lids, lakes = utk_preprocess_boundaries()

for i in 1:length(s_xy)
	s = s_xy[i]
	for j in 1:length(red_i[i])
		PyPlot.fill(s[1][red_i[i][j]],s[2][red_i[i][j]],color=(.9,.9,.9,1.))
		PyPlot.plot(s[1][red_i[i][j]],s[2][red_i[i][j]],"-k",linewidth=.5)
	end
end

# #=
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
PyPlot.plot(xy[fs1,1],xy[fs1,2],"o",color=cols[1],markersize=10.)
PyPlot.plot(xy[fs2,1],xy[fs2,2],"o",color=cols[3],markersize=10.)
PyPlot.plot(xy[fs3,1],xy[fs3,2],"o",color=cols[2],markersize=10.)
axis([-105.,-63.,24.,48.])
PyPlot.xticks([])
PyPlot.yticks([])

subplot(1,2,2)

L0red = nL0[[1:fs1-1;fs1+1:fs2-1;fs2+1:fs3-1;fs3+1:size(nL0)[1]],:]
L0max = [maximum(L0red[:,i]) for i in 1:size(L0red)[2]]
L0min = [minimum(L0red[:,i]) for i in 1:size(L0red)[2]]
PyPlot.fill([ks;ks[end:-1:1]]/T,-[L0max;L0min[end:-1:1]],color="C7")
PyPlot.plot(ks/T,-nL0[fs3,:],"-",color=cols[2],linewidth=2.,label="$fs3")
PyPlot.plot(ks/T,-nL0[fs1,:],"-",color=cols[1],linewidth=2.,label="$fs1")
PyPlot.plot(ks/T,-nL0[fs2,:],"-",color=cols[3],linewidth=2.,label="$fs2")
axis([0.,maximum(ks)/T,-.1,1.1])
xlabel("freq")
ylabel("normalized log-likelihodd")
legend()



#=
subplot2grid((2,4),(0,2),colspan=1,rowspan=2)
for i in 1:length(ls)
	PyPlot.plot(ks/T,L0[i,:],"-o",color=cmap((i-1)/(length(ls)-1)))
end
xlabel("freq")
ylabel("obj")

subplot2grid((2,4),(0,3),colspan=1,rowspan=1)
PyPlot.plot(ks/T,L1,"-o",color=cmap(.4))
xlabel("freq")
ylabel("obj")

subplot2grid((2,4),(1,3),colspan=1,rowspan=1)
PyPlot.plot(ls,γ1[iii],"-o",color=cmap(0.))
xlabel("node id")
ylabel("amplitude")
=#
