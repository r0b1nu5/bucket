using PyPlot, DelimitedFiles, FFTW, LinearAlgebra

ntw = "ieee57"
ex = 2
ks = 1:160
ls = 1:57
ksF = 1:160
τ = .1
fs = 17
ff = 1.85/2π

Xs = readdlm("data_melvyn/"*ntw*"/"*ntw*"_ex$(ex)_Xs.csv",',')

nn,N = size(Xs)
n = Int64(nn/2)
T = (N-1)*τ


L0 = zeros(n,length(ks))
for i in ls
	for j in 1:length(ks)
		L0[i,j] = readdlm("data_melvyn/"*ntw*"/"*ntw*"_ex$(ex)_l0_$(i).$(ks[j])_obj.csv",',')[1]
	end
end
nL0 = (L0[ls,:] .- maximum(L0[ls,:]))./(maximum(L0[ls,:]) - minimum(L0[ls,:]))

L1 = Vector{Float64}()
γ1 = Vector{Vector{Float64}}()
for j in 1:length(ks)
	push!(L1,readdlm("data_melvyn/"*ntw*"/"*ntw*"_ex$(ex)_l1_$(ks[j])_obj.csv",',')[1])
	push!(γ1,vec(readdlm("data_melvyn/"*ntw*"/"*ntw*"_ex$(ex)_l1_$(ks[j])_gamma.csv",',')))
end
nL1 = (L1 .- maximum(L1))./(maximum(L1) - minimum(L1))

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

FX = Matrix{Complex{Float64}}(undef,nn,N-1)
nFX = Matrix{Float64}(undef,nn,N-1)

for i in 1:nn
	FX[i,:] = fft(Xs[i,:])[2:end]
	nFX[i,:] = norm.(FX[i,:])
end

nFX[1:n,:] ./= maximum(nFX[1:n,:])
nFX[n+1:nn,:] ./= maximum(nFX[n+1:nn,:])

xxx = findmax(nFX[1:n,:])
Fsolx = xxx[2][1]
Fxred = nFX[1:n,:]
Fxmax = [maximum(Fxred[:,i]) for i in 1:size(Fxred)[2]]
Fxmin = [minimum(Fxred[:,i]) for i in 1:size(Fxred)[2]]

xxx = findmax(nFX[n+1:nn,:])
Fsolp = xxx[2][1]
Fpred = nFX[n+1:nn,:]
Fpmax = [maximum(Fpred[:,i]) for i in 1:size(Fpred)[2]]
Fpmin = [minimum(Fpred[:,i]) for i in 1:size(Fpred)[2]]

xxx = findmin(nL0)
sol0 = xxx[2][1]
nL0red = nL0[[1:sol-1;sol+1:size(nL0)[1]],:]
nL0max = [maximum(nL0red[:,i]) for i in 1:size(nL0red)[2]]
nL0min = [minimum(nL0red[:,i]) for i in 1:size(nL0red)[2]]

xxx = findmin(nL1)
solk1 = xxx[2]
xxx = findmax(γ1[solk1])
soll1 = xxx[2]
#nγ1 = (γ1[solk1] .- minimum(γ1[solk1]))./(maximum(γ1[solk1]) - minimum(γ1[solk1]))
nγ1 = γ1[solk1]./maximum(γ1[solk1])


figure(ntw*"FT vs. l0/l1 (sup)",(10,9))

#subplot(5,1,1)
subplot2grid((5,20),(0,0),rowspan=1,colspan=20)

PyPlot.plot([ff,ff],[-.1,1.1],"--",color="C7")
PyPlot.fill([ksF;ksF[end:-1:1]]/(N*τ),[Fxmax[ksF];Fxmin[ksF[end:-1:1]]],color=gra)
PyPlot.plot(ksF/(N*τ),nFX[Fsolx,ksF],color=cols[2])
PyPlot.plot(ksF/(N*τ),nFX[sol,ksF],color=cols[1])
axis([ksF[1]/(N*τ),ksF[end]/(N*τ),-.1,1.1])
#axis([ksF[1],ksF[end],-.1,1.1])
#xlabel("freq")
#xlable("k")
ylabel("normalized FT(x)")

#subplot(5,1,2)
subplot2grid((5,20),(1,0),rowspan=1,colspan=20)

PyPlot.plot([ff,ff],[-.1,1.1],"--",color="C7")
PyPlot.fill([ksF;ksF[end:-1:1]]/(N*τ),[Fpmax[ksF];Fpmin[ksF[end:-1:1]]],color=gra)
PyPlot.plot(ksF/(N*τ),nFX[n+Fsolp,ksF],color=cols[3])
PyPlot.plot(ksF/(N*τ),nFX[n+sol,ksF],color=cols[1])
axis([ksF[1]/(N*τ),ksF[end]/(N*τ),-.1,1.1])
#axis([ksF[1],ksF[end],-.1,1.1])
#xlabel("freq")
#xlable("k")
ylabel("normalized FT(p)")

#subplot(5,1,3)
subplot2grid((5,20),(2,0),rowspan=1,colspan=20)

PyPlot.plot([ff,ff],[-.1,1.1],"--",color="C7")
PyPlot.fill([ks;ks[end:-1:1]]/(N*τ),-[nL0max;nL0min[end:-1:1]],color="C7")
PyPlot.plot(ks/(N*τ),-nL0[sol,:],color=cols[1])
#PyPlot.fill([1;ks;ks[end:-1:1];1]/(N*τ),-[nLmax[1];nLmax;nLmin[end:-1:1];nLmin[1]],color=gra)
#PyPlot.plot([1;ks]/(N*τ),-[nL0[sol0,1];nL0[sol0,:]],color=cols[1])
axis([ksF[1]/(N*τ),ksF[end]/(N*τ),-.1,1.1])
xlabel("freq")
ylabel("normalized \n log-likelihood")

#subplot(5,1,4)
subplot2grid((5,20),(3,0),rowspan=1,colspan=20)

PyPlot.plot([ff,ff],[-.1,1.1],"--",color="C7")
PyPlot.plot(ks/(N*τ),-nL1,color=cols[1])
axis([ksF[1]/(N*τ),ksF[end]/(N*τ),-.1,1.1])
xlabel("freq")
ylabel("normalized \n log-likelihood")

#subplot(5,1,5)
subplot2grid((5,20),(4,0),rowspan=1,colspan=12)

PyPlot.plot([fs,fs],[-.1,1.1],"--",color="C7")
PyPlot.plot(ls,nγ1,"ok",markersize=4.)
PyPlot.plot(soll1,nγ1[soll1],"o",color=cols[1],markersize=6.)
PyPlot.plot(Fsolx,nγ1[Fsolx],"o",color=cols[2],markersize=6.)
PyPlot.plot(Fsolp,nγ1[Fsolp],"o",color=cols[3],markersize=6.)
xlabel("node index")
ylabel("normalized \n amplitude")
axis([0.,58.,-.1,1.1])

figure(ntw*" (sup): ntw")

xy = readdlm("data_melvyn/"*ntw*"/"*ntw*"_xy.csv",',')
adj = Int64.(readdlm("data_melvyn/"*ntw*"/"*ntw*"_adj.csv",','))

for i in 1:2:size(adj)[1]
	PyPlot.plot(xy[adj[i,1:2],1],xy[adj[i,1:2],2],"k",linewidth=1.)
end
PyPlot.plot(xy[sol0,1],xy[sol0,2],"o",color=cols[1],markersize=8.)
PyPlot.plot(xy[Fsolx,1],xy[Fsolx,2],"o",color=cols[2],markersize=8.)
PyPlot.plot(xy[Fsolp,1],xy[Fsolp,2],"o",color=cols[3],markersize=8.)
for l in ls
	α = nγ1[l]
	PyPlot.plot(xy[l,1],xy[l,2],"o",markersize=5.,color=1 .- α.*(1 .- cols[1]),markeredgecolor="k")
end
xticks([])
yticks([])







