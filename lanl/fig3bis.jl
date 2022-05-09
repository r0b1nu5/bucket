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
gra = (.8,.8,.8,1.)

colid = idate -> (idate == 69)*1 + (idate == 68)*2 + ((idate != 68)*(idate != 69))*3

@info "Plot ebc_x1 and ebc_x2."
@info "x1 = ? (2, 3, 4, 5, or 8)"
ntw1 = "ebc_"*readline()
@info "x2 = ? (2, 3, 4, 5, or 8)"
ntw2 = "ebc_"*readline()

n1, ks1, T1, file1, date1, nL01, j1, k1 = ebc_preprocess_data(ntw1)
δk1 = ks1[2]-ks1[1]
n2, ks2, T2, file2, date2, nL02, j2, k2 = ebc_preprocess_data(ntw2)
δk2 = ks2[2]-ks2[1]
col1 = cols[1]

figure("ebc 2",(14,4))

#subplot2grid((2,5),(0,0),rowspan=2,colspan=3)
subplot(1,3,1)

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

xticks([])
yticks([])

#subplot2grid((2,5),(0,3),rowspan=1,colspan=2)
subplot(1,3,2)

L0red = nL01[[1:j1-1;j1+1:size(nL01)[1]],:]
L0max = [maximum(L0red[:,i]) for i in 1:size(L0red)[2]]
L0min = [minimum(L0red[:,i]) for i in 1:size(L0red)[2]]
PyPlot.fill([ks1[1]-δk1;ks1;ks1[end]+δk1;ks1[end]+δk1;ks1[end:-1:1];ks1[1]-δk1]/T1,-[L0max[1];L0max;L0max[end];L0min[end];L0min[end:-1:1];L0min[1]],color=gra)
for i in 1:length(ks1)
	PyPlot.plot([1,1].*ks1[i]/T1,[0.,-nL01[j1,i]],color=cols[colid(idate1)])
end
PyPlot.plot(ks1/T1,-nL01[j1,:],"o",color=cols[colid(idate1)])
axis([(ks1[1]-δk1/2)/T1,(ks1[end]+δk1/2)/T1,-.1,1.1])
xlabel("freq")
ylabel("normalized log-likelihood")
title(date1)

#subplot2grid((2,5),(1,3),rowspan=1,colspan=2)
subplot(1,3,3)

L0red = nL02[[1:j2-1;j2+1:size(nL02)[1]],:]
L0max = [maximum(L0red[:,i]) for i in 1:size(L0red)[2]]
L0min = [minimum(L0red[:,i]) for i in 1:size(L0red)[2]]
PyPlot.fill([ks2[1]-δk2;ks2;ks2[end]+δk2;ks2[end]+δk2;ks2[end:-1:1];ks2[1]-δk2]/T2,-[L0max[1];L0max;L0max[end];L0min[end];L0min[end:-1:1];L0min[1]],color=gra)
for i in 1:length(ks2)
	PyPlot.plot([1,1].*ks2[i]/T2,[0.,-nL02[j2,i]],color=cols[colid(idate2)])
end
PyPlot.plot(ks2/T2,-nL02[j2,:],"o",color=cols[colid(idate2)])
axis([(ks2[1]-δk2/2)/T2,(ks2[end]+δk2/2)/T2,-.1,1.1])
xlabel("freq")
ylabel("normalized log-likelihood")
title(date2)

