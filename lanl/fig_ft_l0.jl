using PyPlot, DelimitedFiles, FFTW, LinearAlgebra

# #=
ntw = "ieee57"
ex = 1
ks = [1:15;150:205;260:275]
ksF = 1:400
τ = .1
fs = 17
ff = 1.85/2π
# =#

 #=
ntw = "uk_ex1"
ks = ...
ksF = ...
τ = ...
ff = ...
# =#

Xs = readdlm("data_melvyn/"*ntw*"_ex$(ex)_Xs.csv",',')


L0 = zeros(n,length(ks))
for i in 1:n
	for j in 1:length(ks)
#		L0[i,j] = readdlm("data_melvyn/temp/28_test__l0_$(i).$(ks[j])_obj.csv",',')[1]
		L0[i,j] = readdlm("data_melvyn/"*ntw*"_ex$(ex)_l0_$(i).$(ks[j])_obj.csv",',')[1]
	end
end
nL0 = (L0 .- maximum(L0))./(maximum(L0) - minimum(L0))

cmap = get_cmap("plasma")
colshift = .5
cols = [cmap(1-(i+colshift)/(2+colshift)) for i in 0:2]

nn,N = size(Xs)
n = Int64(nn/2)
T = (N-1)*τ

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
Fxred = nFX[[1:Fsolx-1;Fsolx+1:n],:]
Fxmax = [maximum(Fxred[:,i]) for i in 1:size(Fxred)[2]]
Fxmin = [minimum(Fxred[:,i]) for i in 1:size(Fxred)[2]]

xxx = findmax(nFX[n+1:nn,:])
Fsolp = xxx[2][1]
Fpred = nFX[[n+1:n+Fsolp-1;n+Fsolp+1:nn],:]
Fpmax = [maximum(Fpred[:,i]) for i in 1:size(Fpred)[2]]
Fpmin = [minimum(Fpred[:,i]) for i in 1:size(Fpred)[2]]

xxx = findmin(nL0)
sol = xxx[2][1]
nLred = nL0[[1:sol-1;sol+1:n],:]
nLmax = [maximum(nLred[:,i]) for i in 1:size(nLred)[2]]
nLmin = [minimum(nLred[:,i]) for i in 1:size(nLred)[2]]


figure(ntw*"FT vs. l0",(6,8))

subplot(3,1,1)
PyPlot.plot([ff,ff],[-.1,1.1],"--",color="C7")
PyPlot.fill([ksF;ksF[end:-1:1]]/(N*τ),[Fxmax[ksF];Fxmin[ksF[end:-1:1]]],color="C7")
PyPlot.plot(ksF/(N*τ),nFX[Fsolx,ksF],color=cols[2])
axis([ksF[1]/(N*τ),ksF[end]/(N*τ),-.1,1.1])
#axis([ksF[1],ksF[end],-.1,1.1])
xlabel("freq")
#xlable("k")
ylabel("normalized FT(x)")

subplot(3,1,2)
PyPlot.plot([ff,ff],[-.1,1.1],"--",color="C7")
PyPlot.fill([ksF;ksF[end:-1:1]]/(N*τ),[Fpmax[ksF];Fpmin[ksF[end:-1:1]]],color="C7")
PyPlot.plot(ksF/(N*τ),nFX[n+Fsolp,ksF],color=cols[2])
axis([ksF[1]/(N*τ),ksF[end]/(N*τ),-.1,1.1])
#axis([ksF[1],ksF[end],-.1,1.1])
xlabel("freq")
#xlable("k")
ylabel("normalized FT(p)")

subplot(3,1,3)
PyPlot.plot([ff,ff],[-.1,1.1],"--",color="C7")
PyPlot.fill([ks;ks[end:-1:1]]/(N*τ),-[nLmax;nLmin[end:-1:1]],color="C7")
PyPlot.plot(ks/(N*τ),-nL0[sol,:],color=cols[2])
axis([ksF[1]/(N*τ),ksF[end]/(N*τ),-.1,1.1])
xlabel("freq")
ylabel("normalized log-likelihood")


figure(ntw*": ntw")

xy = readdlm("data_melvyn/"*ntw*"_xy.csv",',')
adj = Int64.(readdlm("data_melvyn/"*ntw*"_adj.csv",','))

for i in 1:2:size(adj)[1]
	PyPlot.plot(xy[adj[i,1:2],1],-xy[adj[i,1:2],2],"k",linewidth=1.)
end
PyPlot.plot(xy[:,1],-xy[:,2],"ok",markersize=5.)
PyPlot.plot(xy[fs,1],-xy[fs,2],"o",color=cols[2],markersize=8.)
xticks([])
yticks([])
