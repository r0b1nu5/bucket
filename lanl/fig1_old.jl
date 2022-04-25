using PyPlot, DelimitedFiles, FFTW, LinearAlgebra

run = 1

Xs = readdlm("data_marc/Xs$(run).csv",',')
nn,N = size(Xs)
n = Int64(nn/2)
τ = .001
T = (N-1)*τ

FX = Array{Complex{Float64},2}(undef,nn,N)

for i in 1:nn
	FX[i,:] = fft(Xs[i,:])
end

ks = 1:50
L0 = zeros(n,length(ks))
for i in 1:n
	for j in 1:length(ks)
		L0[i,j] = readdlm("data/ntw3_$(run)_l0_$(i).$(ks[j])_obj.csv",',')[1]
	end
end

L1 = Array{Float64,1}()
γ1 = Array{Array{Float64,1},1}()
for j in 1:length(ks)
	push!(L1,readdlm("data/ntw3_$(run)_l1_$(ks[j])_obj.csv",',')[1])
	push!(γ1,vec(readdlm("data/ntw3_$(run)_l1_$(ks[j])_g.csv",',')))
end

Lmi = minimum(L1)
Lma = maximum(L1)

cmap = get_cmap("plasma")

figure("fig1_$(run)",(19,4.5))

subplot2grid((2,44),(0,0),colspan=10,rowspan=1)
for i in 1:3
	PyPlot.plot((0:50)/(50001*τ),norm.(FX[i,1:51]),color=cmap((i-1)/2))
end
PyPlot.text(.95,5500.,"(a)")

subplot2grid((2,44),(1,0),colspan=10,rowspan=1)
for i in 1:3
	PyPlot.plot((0:50)/(50001*τ),norm.(FX[i+3,1:51]),color=cmap((i-1)/2))
end
PyPlot.text(.95,19000.,"(b)")

subplot2grid((2,44),(0,11),colspan=10,rowspan=2)
for i in 1:3
	PyPlot.plot(ks/T,L0[i,:],"-",color=cmap((i-1)/2))
end
PyPlot.text(.95,-14.57,"(c)")
xlabel("freq")
ylabel("obj")

subplot2grid((2,44),(0,22),colspan=10,rowspan=2)
PyPlot.plot(ks/T,L1,color=cmap(.4))
#= 
for i in 1:length(L1)
	α = cmap(1 - (L1[i] - Lma)/(Lmi - Lma))
	PyPlot.plot(ks[i]/T,L1[i],".",color=α)
end
=#
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





