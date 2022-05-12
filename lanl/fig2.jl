using PyPlot, DelimitedFiles, FFTW, LinearAlgebra

ntw = "ntw20"
ex = 1
show_graph = false
ks = Int64[]
ls = Int64[]
ksF = Int64[]
τ = .1
fs = 1
ff = 1/2π

if ntw == "uk" && ex == 1
	ks = [80:10:220;222:2:238;240:10:330]
	ls = 1:120
	ksF = 1:205
	fs = 14
	ff = 2.41/2π
elseif ntw == "uk" && ex == 2
	ks = [5:10:105;110:120;125:10:205]
	ls = 1:120
	ksF = 1:205
	fs = 14
	ff = 2.41/2π
elseif ntw == "ieee57" && ex == 1
	ks = [1:205;260:275;275:10:330]
	ls = 1:57
	ksF = 1:160
	fs = 17
	ff = 1.85/2π
elseif ntw == "ieee57" && ex == 2
	ks = 10:160
	ls = 1:57
	ksF = 1:160
	fs = 17
	ff = 1.85/2π
elseif ntw == "ntw20" && ex == 1
	ks = 1:40
	ls = 1:20
	ksF = 1:40
	fs = 20
	ff = .08
end

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
sol = xxx[2][1]
nLred = nL0[[1:sol-1;sol+1:size(nL0)[1]],:]
nLmax = [maximum(nLred[:,i]) for i in 1:size(nLred)[2]]
nLmin = [minimum(nLred[:,i]) for i in 1:size(nLred)[2]]


figure(ntw*"FT vs. l0",(10,5))

subplot(3,1,1)
PyPlot.plot([ff,ff],[-.1,1.1],"--",color="C7")
PyPlot.fill([ksF;ksF[end:-1:1]]/(N*τ),[Fxmax[ksF];Fxmin[ksF[end:-1:1]]],color=gra)
PyPlot.plot(ksF/(N*τ),nFX[Fsolx,ksF],color=cols[2])
PyPlot.plot(ksF/(N*τ),nFX[sol,ksF],color=cols[1])
axis([ksF[1]/(N*τ),ksF[end]/(N*τ),-.1,1.1])
#axis([ksF[1],ksF[end],-.1,1.1])
#xlabel("freq")
#xlable("k")
ylabel("normalized FT(x)")

subplot(3,1,2)
PyPlot.plot([ff,ff],[-.1,1.1],"--",color="C7")
PyPlot.fill([ksF;ksF[end:-1:1]]/(N*τ),[Fpmax[ksF];Fpmin[ksF[end:-1:1]]],color=gra)
PyPlot.plot(ksF/(N*τ),nFX[n+Fsolp,ksF],color=cols[3])
PyPlot.plot(ksF/(N*τ),nFX[n+sol,ksF],color=cols[1])
axis([ksF[1]/(N*τ),ksF[end]/(N*τ),-.1,1.1])
#axis([ksF[1],ksF[end],-.1,1.1])
#xlabel("freq")
#xlable("k")
ylabel("normalized FT(p)")

subplot(3,1,3)
PyPlot.plot([ff,ff],[-.1,1.1],"--",color="C7")
#PyPlot.fill([ks;ks[end:-1:1]]/(N*τ),-[nLmax;nLmin[end:-1:1]],color="C7")
#PyPlot.plot(ks/(N*τ),-nL0[sol,:],color=cols[2])
PyPlot.fill([1;ks;ks[end:-1:1];1]/(N*τ),-[nLmax[1];nLmax;nLmin[end:-1:1];nLmin[1]],color=gra)
PyPlot.plot([1;ks]/(N*τ),-[nL0[sol,1];nL0[sol,:]],color=cols[1])
axis([ksF[1]/(N*τ),ksF[end]/(N*τ),-.1,1.1])
xlabel("freq")
ylabel("normalized \n log-likelihood")


if show_graph
	figure(ntw*": ntw")

	xy = readdlm("data_melvyn/"*ntw*"/"*ntw*"_xy.csv",',')
	adj = Int64.(readdlm("data_melvyn/"*ntw*"/"*ntw*"_adj.csv",','))
	
	for i in 1:2:size(adj)[1]
		PyPlot.plot(xy[adj[i,1:2],1],xy[adj[i,1:2],2],"k",linewidth=1.)
	end
	PyPlot.plot(xy[:,1],xy[:,2],"ok",markersize=5.)
	PyPlot.plot(xy[sol,1],xy[sol,2],"o",color=cols[1],markersize=8.)
	PyPlot.plot(xy[Fsolx,1],xy[Fsolx,2],"o",color=cols[2],markersize=8.)
	PyPlot.plot(xy[Fsolp,1],xy[Fsolp,2],"o",color=cols[3],markersize=8.)
	xticks([])
	yticks([])
end







