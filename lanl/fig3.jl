using PyPlot, DelimitedFiles

include("processing_tools.jl")

cmap = get_cmap("plasma")
colshift = .5
cols = [cmap(1-(i+colshift)/(2+colshift)) for i in 0:2]

@info "Plot pen_x1 and pen_x2."
@info "x1 = ? (2, 3, 4, 5, or 8)"
ntw1 = "pen_"*readline()

n1, ks1, T1, file1, date1, nL01, j1, k1 = pen_preprocess_data(ntw1)

figure("Dan",(14,4.5))

subplot(1,2,1)

ixy = readdlm("data_pen/deathvalley_coord.dat")
idx = Int64.(ixy[:,1])
x = ixy[:,2]
y = ixy[:,3]
PyPlot.plot(x,y,"ok",markersize=3.)

idate1 = Int64.(vec(readdlm("data_pen/pen_"*date1*"_ids.csv",',')))[j1] - 1
a,row1 = findmin(abs.(idx .- idate1))
PyPlot.plot(x[row1],y[row1],"o",color=cols[2],markersize=8.)
xticks([])
yticks([])

subplot(1,2,2)

L0red = nL01[[1:j1-1;j1+1:size(nL01)[1]],:]
L0max = [maximum(L0red[:,i]) for i in 1:size(L0red)[2]]
L0min = [minimum(L0red[:,i]) for i in 1:size(L0red)[2]]
PyPlot.fill([ks1;ks1[end:-1:1]]/T1,-[L0max;L0min[end:-1:1]],color="C7")
PyPlot.plot(ks1/T1,-nL01[j1,:],"-",color=cols[2],linewidth=2.)
axis([ks1[1]/T1,ks1[end]/T1,-.1,1.1])
xlabel("freq")
ylabel("normalized log-likelihodd")
title(ntw1)


