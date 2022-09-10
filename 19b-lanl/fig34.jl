using PyPlot, DelimitedFiles

include("processing_tools.jl")

cmap = get_cmap("plasma")
colshift1 = .5
colshift2 = .5
cols = [cmap(1-(i+colshift1)/(2+colshift1+colshift2)) for i in 0:2]

@info "Plot ebc_x and utk."
@info "x = ? (2, 3, 4, 5, or 8)"
ntw1 = "ebc_"*readline()

n1, ks1, T1, file1, date1, nL01, j1, k1 = ebc_preprocess_data(ntw1)
AAA = nL01

figure("Dan - UTK",(11,8))

subplot2grid((2,11),(0,0),colspan=4,rowspan=1)

ixy = readdlm("data_ebc/coord.dat")
idx = Int64.(ixy[:,1])
x = ixy[:,2]
y = ixy[:,3]
PyPlot.plot(x,y,"ok",markersize=3.)

idate1 = Int64.(vec(readdlm("data_ebc/ebc_"*date1*"_ids.csv",',')))[j1] - 1
a,row1 = findmin(abs.(idx .- idate1))
PyPlot.plot(x[row1],y[row1],"o",color=cols[2],markersize=8.)
xticks([])
yticks([])

subplot2grid((2,11),(1,0),colspan=4,rowspan=1)

L0red = nL01[[1:j1-1;j1+1:size(nL01)[1]],:]
L0max = [maximum(L0red[:,i]) for i in 1:size(L0red)[2]]
L0min = [minimum(L0red[:,i]) for i in 1:size(L0red)[2]]
PyPlot.fill([ks1;ks1[end:-1:1]]/T1,-[L0max;L0min[end:-1:1]],color="C7")
PyPlot.plot(ks1/T1,-nL01[j1,:],"-",color=cols[2],linewidth=2.)
axis([ks1[1]/T1,ks1[end]/T1,-.1,1.1])
xlabel("freq")
ylabel("normalized log-likelihodd")
title(ntw1)


ntw = "utk"
n = 99
ls = 1:n
ks = 1:1500
T = 300.
file = "mysterious_forcing"
fs1 = 85

L0 = zeros(length(ls),length(ks))
for i in 1:length(ls)
	for j in 1:length(ks)
		L0[i,j] = readdlm("data/utk/"*file*"_l0_$(ls[i]).$(ks[j])_obj.csv",',')[1]
	end
end
nL0 = (L0 .- maximum(L0))./(maximum(L0) - minimum(L0))


##############################################################

subplot2grid((2,11),(0,5),colspan=6,rowspan=1)

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

xy = readdlm("data_utk/utk1_coord.csv",',')
ids = (1:size(xy)[1]).*(xy[:,1] .< -55.)
ids = setdiff(ids,[0,])
PyPlot.plot(xy[ids,1],xy[ids,2],"ok",markersize=5.)
PyPlot.plot(xy[fs1,1],xy[fs1,2],"o",color=cols[2],markersize=10.)
axis([-105.,-63.,24.,48.])
PyPlot.xticks([])
PyPlot.yticks([])

subplot2grid((2,11),(1,5),colspan=6,rowspan=1)

L0red = nL0[[1:fs1-1;fs1+1:size(nL0)[1]],:]
L0max = [maximum(L0red[:,i]) for i in 1:size(L0red)[2]]
L0min = [minimum(L0red[:,i]) for i in 1:size(L0red)[2]]
PyPlot.fill([ks;ks[end:-1:1]]/T,-[L0max;L0min[end:-1:1]],color="C7")
PyPlot.plot(ks/T,-nL0[fs1,:],"-",color=cols[2],linewidth=2.)
axis([0.,maximum(ks)/T,-.1,1.1])
xlabel("freq")
ylabel("normalized log-likelihodd")
legend()

