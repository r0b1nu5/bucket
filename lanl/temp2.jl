using PyPlot, DelimitedFiles, FFTW, LinearAlgebra

figure(k,(10,7))

#for k in 1:56
#Xs = readdlm("data_melvyn/temp/$(k)_1.85_data",',')
Xs = readdlm("data_melvyn/temp/$(k)_2.0_data",',')
nn,N = size(Xs)
n = Int64(nn/2)
τ = .1
T = (N-1)*τ
ff = 2.0/2π

ks = 180:200

FX = Matrix{Complex{Float64}}(undef,nn,N-1)
nFX = Matrix{Float64}(undef,nn,N-1)

for i in 1:nn
	FX[i,:] = fft(Xs[i,:])[2:end]
	nFX[i,:] = norm.(FX[i,:])
end

nFX[1:n,:] ./= maximum(nFX[1:n,:])
nFX[n+1:nn,:] ./= maximum(nFX[n+1:nn,:])




#=
ks = 6:10
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
=#

cmap = get_cmap("plasma")
colshift = .5
cols = [cmap(1-(i+colshift)/(2+colshift)) for i in 0:2]


#subplot(1,2,1)
PyPlot.plot([ff,ff],[-.1,1.1],"--",color="C7")
for i in 1:n
	PyPlot.plot(ks/(N*τ),nFX[i,ks],color=cmap(i/n))
	PyPlot.plot(ks/(N*τ),nFX[i+n,ks],"--",color=cmap(i/n))
end
#axis([0.,50/(N*τ),-.1,1.1])
xlabel("freq")
ylabel("normalized FT")

PyPlot.plot(ks/(N*τ),nFX[k,ks],"k")
PyPlot.plot(ks/(N*τ),nFX[k+n,ks],"--k")

#end
#=
subplot(1,2,2)
PyPlot.plot([ff,ff],[-1.1,.1],"--",color="C7")
for i in 1:n
	PyPlot.plot(ks/T,nL0[i,:],"o-",color=cmap(i/n))
end
#axis([35/(N*τ),45/(N*τ),-1.1,.1])
xlabel("freq")
ylabel("normalized inverse log-likelihood")

=#
