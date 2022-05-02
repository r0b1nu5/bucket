using PyPlot, DelimitedFiles, FFTW, LinearAlgebra

run = 1

L = readdlm("data_marc/L.csv",',')
us = eigvecs(L)
λs = round.(eigvals(L))

Xs = readdlm("data_marc/Xs$(run).csv",',')
nn,N = size(Xs)
n = Int64(nn/2)
τ = .001
T = (N-1)*τ
ff = 1/2π

FX = Matrix{Complex{Float64}}(undef,nn,N)
nFX = Matrix{Float64}(undef,nn,N)

for i in 1:nn
	FX[i,:] = fft(Xs[i,:])
	nFX[i,:] = norm.(FX[i,:])
end

nFX[1:n,:] ./= maximum(nFX[1:n,:])
nFX[n+1:nn,:] ./= maximum(nFX[n+1:nn,:])

ks = 1:50
L0 = zeros(n,length(ks))
for i in 1:n
	for j in 1:length(ks)
		L0[i,j] = readdlm("data/ntw3_$(run)_l0_$(i).$(ks[j])_obj.csv",',')[1]
	end
end
nL0 = (L0 .- maximum(L0))./(maximum(L0) - minimum(L0))

 #=
L1 = Array{Float64,1}()
γ1 = Array{Array{Float64,1},1}()
for j in 1:length(ks)
	push!(L1,readdlm("data/ntw3_$(run)_l1_$(ks[j])_obj.csv",',')[1])
	push!(γ1,vec(readdlm("data/ntw3_$(run)_l1_$(ks[j])_g.csv",',')))
end

Lmi = minimum(L1)
Lma = maximum(L1)
# =#

cmap = get_cmap("plasma")
colshift1 = .5
colshift2 = .5
cols = [cmap(1-(i+colshift1)/(2+colshift1+colshift2)) for i in 0:2]

 #=
figure("fig1",(15,4.5))

subplot(1,3,1)
x = LinRange(0,2π,4)
lws = [4.,2.,4.]
lts = ["-k","-k","--k"]
bss = [40,20,10]
for i in 1:3
	PyPlot.plot(sin.(x[[i,i+1]]),cos.(x[[i,i+1]]),lts[i],linewidth=lws[i])
end
for i in 1:3
	PyPlot.plot(sin(x[i]),cos(x[i]),"o",color=cols[i],markersize=bss[i])
end
PyPlot.plot(LinRange(-1.,-.5,200),.1*sin.(LinRange(0.,4π,200)) .+ 1.,color="C7")
axis([-1.,1.,-.7,1.3])
#axis("off")

subplot(1,3,2)
PyPlot.plot([ff,ff],[-.1,1.1],"--",color="C7")
for i in 1:3
	PyPlot.plot((0:50)/(50001*τ),nFX[i,1:51],color=cols[i])
	PyPlot.plot((0:50)/(50001*τ),nFX[i+3,1:51],"--",color=cols[i])
end
axis([0.,50/(50001*τ),-.1,1.1])
xlabel("freq")
ylabel("normalized FT")

subplot(1,3,3)
PyPlot.plot([ff,ff],[-1.1,.1],"--",color="C7")
for i in 1:3
	PyPlot.plot(ks/T,nL0[i,:],color=cols[i])
end
axis([0.,50/(50001*τ),-1.1,.1])
xlabel("freq")
ylabel("normalized inverse log-likelihood")

# =#

figure("test",(18,7))

subplot2grid((2,3),(0,0),colspan=1,rowspan=1)
for i in 1:3
	PyPlot.plot((0:N-1)*τ,Xs[i,:],color=cols[i])
end
xlabel("t[s]")
ylabel("x(t)")

subplot2grid((2,3),(0,1),colspan=1,rowspan=1)
for i in 1:3
	PyPlot.plot((0:N-1)*τ,Xs[i+3,:],color=cols[i])
end
xlabel("t[s]")
ylabel("p(t)")

subplot2grid((2,3),(0,2),colspan=1,rowspan=1)
PyPlot.plot([ff,ff],[-.1,1.1],"--",color="C7")
for i in 1:3
	PyPlot.plot((0:50)/(N*τ),nFX[i,1:51],color=cols[i])
	PyPlot.plot((0:50)/(N*τ),nFX[i+3,1:51],"--",color=cols[i])
end
xlabel("freq")
ylabel("normalized FT")
axis([0.,50/(50001*τ),-.1,1.1])

subplot2grid((2,3),(1,0),colspan=2,rowspan=1)
θ = LinRange(0,2π,4)
α = 2.5
#ms = (2.).^(log.(abs.(us)) .+ 5)
ms = 5 .+ 15*(abs.(us) .- minimum(abs.(us)))./(maximum(abs.(us)) - minimum(abs.(us)))
for v in 1:3
	PyPlot.plot(sin.(θ) .+ α*v,cos.(θ),"k")
	PyPlot.text(α*v,0.,"λ = $(λs[v])")
	for i in 1:3
		PyPlot.plot(sin(θ[i]) + α*v,cos(θ[i]),"o",color=cols[i],markersize=ms[i,v])
	end
end
axis([1.,9.,-.9,1.4])
xticks([])
yticks([])

subplot2grid((2,3),(1,2),colspan=1,rowspan=1)
PyPlot.plot([ff,ff],[-.1,1.1],"--",color="C7")
for i in 1:3
	PyPlot.plot(ks/T,-nL0[i,:],color=cols[i])
end
xlabel("freq")
ylabel("normalized log-likelihood")
axis([0.,50/(50001*τ),-.1,1.1])

#=

subplot2grid((2,44),(0,22),colspan=10,rowspan=2)
PyPlot.plot(ks/T,L1,color=cmap(.4))
PyPlot.text(.95,-14.57,"(d)")
xlabel("freq")
ylabel("obj")

subplot2grid((2,44),(0,33),colspan=10,rowspan=2)
for i in 1:length(γ1)
	α = cmap(1 - (L1[i] - Lma)/(Lmi - Lma))
	PyPlot.plot(1:length(γ1[i]),γ1[i],"o",color=α)
end
PyPlot.text(2.9,.85,"(e)")
xlabel("node id")
ylabel("amplitude")

subplot2grid((2,44),(0,43),colspan=1,rowspan=2)
yticks([])
twinx()
for ll in LinRange(Lmi,Lma,400)
	α = cmap(1 - (ll - Lma)/(Lmi - Lma))
	PyPlot.plot([0,1],[ll,ll],color=α)
end
axis([0,1,Lmi,Lma])
xticks([])
ylabel("obj")

=#





