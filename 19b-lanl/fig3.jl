using PyPlot, DelimitedFiles

include("processing_tools.jl")

 #=
cmap = get_cmap("plasma")
colshift1 = .5
colshift2 = .5
cols = [cmap(1-(i+colshift1)/(2+colshift1+colshift2)) for i in 0:2]
# =#
# #=
cols = [(243,111,33)./255,(0,145,194)./255,(161,0,202)./255]
# =#

@info "Plot ebc_x1 and ebc_x2."
@info "x1 = ? (2, 3, 4, 5, or 8)"
ntw1 = "ebc_"*readline()

n1, ks1, T1, file1, date1, nL01, j1, k1 = ebc_preprocess_data(ntw1)

figure("ebc 1",(11,5))

subplot(1,2,1)

#ixy = readdlm("data_ebc/coord.dat")
ixy = readdlm("data_ebc/coord2.csv",',')
idx = Int64.(ixy[:,1])
x = ixy[:,2]
y = ixy[:,3]
PyPlot.plot(x,y,"ok",markersize=3.)

idate1 = Int64.(vec(readdlm("data_ebc/ebc_"*date1*"_ids.csv",',')))[j1] - 1
a,row1 = findmin(abs.(idx .- idate1))
PyPlot.plot(x[row1],y[row1],"o",color=cols[2],markersize=8.)
xticks([])
yticks([])

subplot(1,2,2)

L0red = nL01[[1:j1-1;j1+1:size(nL01)[1]],:]
L0max = [maximum(L0red[:,i]) for i in 1:size(L0red)[2]]
L0min = [minimum(L0red[:,i]) for i in 1:size(L0red)[2]]
PyPlot.fill([ks1;ks1[end:-1:1]]/T1,-[L0max;L0min[end:-1:1]],color="C7")
for i in 1:length(ks1)
	PyPlot.plot([1,1].*(ks1[i]/T1),[0.,-nL01[j1,i]],color=cols[2])
end
PyPlot.plot(ks1/T1,-nL01[j1,:],"o",color=cols[2])
axis([ks1[1]/T1,ks1[end]/T1,-.1,1.1])
xlabel("freq")
ylabel("normalized log-likelihodd")
title(ntw1)


