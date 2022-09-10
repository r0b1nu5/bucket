using PyPlot, DelimitedFiles

include("processing_tools.jl")

cmap = get_cmap("plasma")
colshift1 = .5
colshift2 = .5
cols = [cmap(1-(i+colshift1)/(2+colshift1+colshift2)) for i in 0:2]

colid = idate -> (idate == 69)*2 + (idate == 68)*1 + ((idate != 68)*(idate != 69))*3

@info "Plot four indices out of five."
@info "Which one to exclude? (2, 3, 4, 5, or 8)"
ntw5 = "ebc_"*readline()
ntws = setdiff(["ebc_2","ebc_3","ebc_4","ebc_5","ebc_8"],[ntw5,])
ntw1 = ntws[1]
ntw2 = ntws[2]
ntw3 = ntws[3]
ntw4 = ntws[4]

n1, ks1, T1, file1, date1, nL01, j1, k1 = ebc_preprocess_data(ntw1)
n2, ks2, T2, file2, date2, nL02, j2, k2 = ebc_preprocess_data(ntw2)
n3, ks3, T3, file3, date3, nL03, j3, k3 = ebc_preprocess_data(ntw3)
n4, ks4, T4, file4, date4, nL04, j4, k4 = ebc_preprocess_data(ntw4)

figure("ebc 4",(12,5))

subplot2grid((2,7),(0,0),rowspan=2,colspan=3)

ixy = readdlm("data_ebc/coord2.csv",',')
idx = Int64.(ixy[:,1])
x = ixy[:,2]
y = ixy[:,3]
PyPlot.plot(x,y,"ok",markersize=3.)

idate1 = Int64.(vec(readdlm("data_ebc/ebc_"*date1*"_ids.csv",',')))[j1] - 1
a,row1 = findmin(abs.(idx .- idate1))
PyPlot.plot(x[row1],y[row1],"o",color=cols[colid(idate1)],markersize=5.)

idate2 = Int64.(vec(readdlm("data_ebc/ebc_"*date2*"_ids.csv",',')))[j2] - 1
a,row2 = findmin(abs.(idx .- idate2))
PyPlot.plot(x[row2],y[row2],"o",color=cols[colid(idate2)],markersize=5.)

idate3 = Int64.(vec(readdlm("data_ebc/ebc_"*date3*"_ids.csv",',')))[j3] - 1
a,row3 = findmin(abs.(idx .- idate3))
PyPlot.plot(x[row3],y[row3],"o",color=cols[colid(idate3)],markersize=5.)

idate4 = Int64.(vec(readdlm("data_ebc/ebc_"*date4*"_ids.csv",',')))[j4] - 1
a,row4 = findmin(abs.(idx .- idate4))
PyPlot.plot(x[row4],y[row4],"o",color=cols[colid(idate4)],markersize=5.)

xticks([])
yticks([])

subplot2grid((2,7),(0,3),rowspan=1,colspan=2)

L0red = nL01[[1:j1-1;j1+1:size(nL01)[1]],:]
L0max = [maximum(L0red[:,i]) for i in 1:size(L0red)[2]]
L0min = [minimum(L0red[:,i]) for i in 1:size(L0red)[2]]
PyPlot.fill([ks1;ks1[end:-1:1]]/T1,-[L0max;L0min[end:-1:1]],color="C7")
for i in 1:length(ks1)
	PyPlot.plot([1,1].*ks1[i]/T1,[0.,-nL01[j1,i]],color=cols[colid(idate1)])
end
PyPlot.plot(ks1/T1,-nL01[j1,:],"o",color=cols[colid(idate1)])
axis([ks1[1]/T1,ks1[end]/T1,-.1,1.1])
xlabel("freq")
ylabel("normalized log-likelihood")
title(ntw1)

subplot2grid((2,7),(1,3),rowspan=1,colspan=2)

L0red = nL02[[1:j2-1;j2+1:size(nL02)[1]],:]
L0max = [maximum(L0red[:,i]) for i in 1:size(L0red)[2]]
L0min = [minimum(L0red[:,i]) for i in 1:size(L0red)[2]]
PyPlot.fill([ks2;ks2[end:-1:1]]/T2,-[L0max;L0min[end:-1:1]],color="C7")
for i in 1:length(ks2)
	PyPlot.plot([1,1].*ks2[i]/T2,[0.,-nL02[j2,i]],color=cols[colid(idate2)])
end
PyPlot.plot(ks2/T2,-nL02[j2,:],"o",color=cols[colid(idate2)])
axis([ks2[1]/T2,ks2[end]/T2,-.1,1.1])
xlabel("freq")
ylabel("normalized log-likelihood")
title(ntw2)

subplot2grid((2,7),(0,5),rowspan=1,colspan=2)

L0red = nL03[[1:j3-1;j3+1:size(nL03)[1]],:]
L0max = [maximum(L0red[:,i]) for i in 1:size(L0red)[2]]
L0min = [minimum(L0red[:,i]) for i in 1:size(L0red)[2]]
PyPlot.fill([ks3;ks3[end:-1:1]]/T3,-[L0max;L0min[end:-1:1]],color="C7")
for i in 1:length(ks3)
	PyPlot.plot([1,1].*ks3[i]/T3,[0.,-nL03[j3,i]],color=cols[colid(idate3)])
end
PyPlot.plot(ks3/T3,-nL03[j3,:],"o",color=cols[colid(idate3)])
axis([ks3[1]/T3,ks3[end]/T3,-.1,1.1])
xlabel("freq")
ylabel("normalized log-likelihood")
title(ntw3)

subplot2grid((2,7),(1,5),rowspan=1,colspan=2)

L0red = nL04[[1:j4-1;j4+1:size(nL04)[1]],:]
L0max = [maximum(L0red[:,i]) for i in 1:size(L0red)[2]]
L0min = [minimum(L0red[:,i]) for i in 1:size(L0red)[2]]
PyPlot.fill([ks4;ks4[end:-1:1]]/T4,-[L0max;L0min[end:-1:1]],color="C7")
for i in 1:length(ks4)
	PyPlot.plot([1,1].*ks4[i]/T4,[0.,-nL04[j4,i]],color=cols[colid(idate4)])
end
PyPlot.plot(ks4/T4,-nL04[j4,:],"o",color=cols[colid(idate4)])
axis([ks4[1]/T4,ks4[end]/T4,-.1,1.1])
xlabel("freq")
ylabel("normalized log-likelihood")
title(ntw4)

