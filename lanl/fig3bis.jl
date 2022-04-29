using PyPlot, DelimitedFiles

include("processing_tools.jl")

cmap = get_cmap("plasma")
colshift = .5
cols = [cmap(1-(i+colshift)/(2+colshift)) for i in 0:2]

@info "Plot pen_x1 and pen_x2."
@info "x1 = ? (2, 3, 4, 5, or 8)"
ntw1 = "pen_"*readline()
@info "x2 = ? (2, 3, 4, 5, or 8)"
ntw2 = "pen_"*readline()

n1, ks1, T1, file1, date1, nL01, j1, k1 = pen_preprocess_data(ntw1)
n2, ks2, T2, file2, date2, nL02, j2, k2 = pen_preprocess_data(ntw2)

figure("Dan",(15,6))

subplot2grid((2,3),(0,0),rowspan=2,colspan=2)

ixy = readdlm("data_pen/deathvalley_coord.dat")
idx = Int64.(ixy[:,1])
x = ixy[:,2]
y = ixy[:,3]
PyPlot.plot(x,y,"ok",markersize=3.)

idate1 = Int64.(vec(readdlm("data_pen/pen_"*date1*"_ids.csv",',')))[j1] - 1
a,row1 = findmin(abs.(idx .- idate1))
PyPlot.plot(x[row1],y[row1],"o",color=cols[1],markersize=6.)

idate2 = Int64.(vec(readdlm("data_pen/pen_"*date2*"_ids.csv",',')))[j2] - 1
a,row2 = findmin(abs.(idx .- idate2))
PyPlot.plot(x[row2],y[row2],"o",color=cols[2],markersize=6.)

subplot2grid((2,3),(0,2),rowspan=1,colspan=2)

L0red = nL01[[1:j1-1;j1+1:size(nL01)[1]],:]
L0max = [maximum(L0red[:,i]) for i in 1:size(L0red)[2]]
L0min = [minimum(L0red[:,i]) for i in 1:size(L0red)[2]]
PyPlot.fill([ks1;ks1[end:-1:1]]/T1,-[L0max;L0min[end:-1:1]],color="C7")
PyPlot.plot(ks1/T1,-nL01[j1,:],"-",color=cols[1],linewidth=2.)
axis([ks1[1]/T1,ks1[end]/T1,-.1,1.1])
xlabel("freq")
ylabel("normalized log-likelihodd")
title(ntw1)

subplot2grid((2,3),(1,2),rowspan=1,colspan=2)

L0red = nL02[[1:j2-1;j2+1:size(nL02)[1]],:]
L0max = [maximum(L0red[:,i]) for i in 1:size(L0red)[2]]
L0min = [minimum(L0red[:,i]) for i in 1:size(L0red)[2]]
PyPlot.fill([ks2;ks2[end:-1:1]]/T2,-[L0max;L0min[end:-1:1]],color="C7")
PyPlot.plot(ks2/T2,-nL02[j2,:],"-",color=cols[2],linewidth=2.)
axis([ks2[1]/T2,ks2[end]/T2,-.1,1.1])
xlabel("freq")
ylabel("normalized log-likelihodd")
title(ntw2)

