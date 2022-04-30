using PyPlot, DelimitedFiles, FFTW, LinearAlgebra

Xs = readdlm("data_melvyn/temp/ieee57_Xs.csv",',')
nn,N = size(Xs)
n = Int64(nn/2)
τ = .1
T = (N-1)*τ
ff = 1.85/2π

FX = Matrix{Complex{Float64}}(undef,nn,N-1)
nFX = Matrix{Float64}(undef,nn,N-1)

for i in 1:nn
	FX[i,:] = fft(Xs[i,:])[2:end]
	nFX[i,:] = norm.(FX[i,:])
end

nFX[1:n,:] ./= maximum(nFX[1:n,:])
nFX[n+1:nn,:] ./= maximum(nFX[n+1:nn,:])

ks = [1:15;150:205;260:275]
ksF = 1:400
L0 = zeros(n,length(ks))
for i in 1:n
	for j in 1:length(ks)
		L0[i,j] = readdlm("data_melvyn/temp/28_test__l0_$(i).$(ks[j])_obj.csv",',')[1]
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
colshift = .5
cols = [cmap(1-(i+colshift)/(2+colshift)) for i in 0:2]

figure("fig1",(7,10))

subplot(2,1,1)
PyPlot.plot([ff,ff],[-.1,1.1],"--",color="C7")
for i in 1:n
	PyPlot.plot(ksF/(N*τ),nFX[i,ksF],color=cmap(i/n))
	PyPlot.plot(ksF/(N*τ),nFX[i+n,ksF],"--",color=cmap(i/n))
#	PyPlot.plot(ksF,nFX[i,ksF],color=cmap(i/n))
#	PyPlot.plot(ksF,nFX[i+n,ksF],"--",color=cmap(i/n))
end
#axis([0.,50/(N*τ),-.1,1.1])
axis([0.,ks[end]/(N*τ),-.1,1.1])
xlabel("freq")
ylabel("normalized FT")

subplot(2,1,2)
PyPlot.plot([ff,ff],[-.1,1.1],"--",color="C7")
for i in 1:n
	PyPlot.plot(ks/T,-nL0[i,:],"o-",color=cmap(i/n))
end
#axis([35/(N*τ),45/(N*τ),-1.1,.1])
axis([0.,ks[end]/(N*τ),-.1,1.1])
xlabel("freq")
ylabel("normalized inverse log-likelihood")


