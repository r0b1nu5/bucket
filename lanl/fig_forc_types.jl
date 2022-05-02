using PyPlot, DelimitedFiles

figure("fig5",(7.,10.))

cmap = get_cmap("plasma")
colshift1 = .5
colshift2 = .5
cols = [cmap(1-(i+colshift1)/(2+colshift1+colshift2)) for i in 0:2]

f = .3
t = LinRange(-1,11,200)

# Multi-sine

subplot(3,2,1)
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
		L[i,j] = readdlm("data/ntw20_multisine_l0_$(ls[i]).$(ks[j])_obj.csv",',')[1]
	end
end
nL = (L .- maximum(L))./(maximum(L) - minimum(L))
nLmi,ii = findmin(nL)
iii = ii[1]
nLmax = [maximum(nL[[1:iii-1;iii+1:size(nL)[1]],i]) for i in 1:size(nL)[2]]
nLmin = [minimum(nL[[1:iii-1;iii+1:size(nL)[1]],i]) for i in 1:size(nL)[2]]
freq = round(ks[ii[2]]/100.,digits=3)

subplot(3,2,2)
PyPlot.fill([ks;ks[end:-1:1]]/100.,-[nLmax;nLmin[end:-1:1]],color="C7")
PyPlot.plot(ks/100.,-nL[iii,:],color=cols[2])
xlabel("freq")
ylabel("normalized log-likelihood")
axis([ks[1]/100.,2.,-.1,1.1])


# Saw
subplot(3,2,3)
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
		L[i,j] = readdlm("data/ntw20_saw_l0_$(ls[i]).$(ks[j])_obj.csv",',')[1]
	end
end
nL = (L .- maximum(L))./(maximum(L) - minimum(L))
nLmi,ii = findmin(nL)
iii = ii[1]
nLmax = [maximum(nL[[1:iii-1;iii+1:size(nL)[1]],i]) for i in 1:size(nL)[2]]
nLmin = [minimum(nL[[1:iii-1;iii+1:size(nL)[1]],i]) for i in 1:size(nL)[2]]
freq = round(ks[ii[2]]/100.,digits=3)

subplot(3,2,4)
PyPlot.fill([ks;ks[end:-1:1]]/100.,-[nLmax;nLmin[end:-1:1]],color="C7")
PyPlot.plot(ks/100.,-nL[iii,:],color=cols[2])
xlabel("freq")
ylabel("normalized log-likelihood")
axis([ks[1]/100.,2.,-.1,1.1])


# Step
subplot(3,2,5)
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
		L[i,j] = readdlm("data/ntw20_step_l0_$(ls[i]).$(ks[j])_obj.csv",',')[1]
	end
end
nL = (L .- maximum(L))./(maximum(L) - minimum(L))
nLmi,ii = findmin(nL)
iii = ii[1]
nLmax = [maximum(nL[[1:iii-1;iii+1:size(nL)[1]],i]) for i in 1:size(nL)[2]]
nLmin = [minimum(nL[[1:iii-1;iii+1:size(nL)[1]],i]) for i in 1:size(nL)[2]]
freq = round(ks[ii[2]]/100.,digits=3)

subplot(3,2,6)
PyPlot.fill([ks;ks[end:-1:1]]/100.,-[nLmax;nLmin[end:-1:1]],color="C7")
PyPlot.plot(ks/100.,-nL[iii,:],color=cols[2])
xlabel("freq")
ylabel("normalized log-likelihood")
axis([ks[1]/100.,2.,-.1,1.1])




