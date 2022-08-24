using PyPlot, DelimitedFiles, FFTW, LinearAlgebra

τ = .1
N = 3000
T = N*τ
n1 = 19
n2 = 20

ls1 = 1:n1
ls2 = 1:n2
ks = 5:205

ff1 = 2.4
ff2 = ff1*.4

ntw = "ntw20"
show_graph = true
cols = [(243,111,33)./255,(0,145,194)./255,(161,0,202)./255,"C2"]
gra = (.8,.8,.8,1.)

# #=
# l0 ===================================
L0 = zeros(n1,length(ks))
for i in ls1
	for j in 1:length(ks)
		L0[i,j] = readdlm("data_melvyn/"*ntw*"/l0_2_missing_hidden_20_test__l0_$(i).$(ks[j])_obj.csv",',')[1]
	end
end
nL0 = (L0[ls1,:] .- maximum(L0[ls1,:]))./(maximum(L0[ls1,:]) - minimum(L0[ls1,:]))

sol11 = 13
sol12 = 14
sol13 = 15
sol14 = 17
nL0red = nL0[[1:12;16;18:19],:]
nL0max = [maximum(nL0red[:,i]) for i in 1:size(nL0red)[2]]
nL0min = [minimum(nL0red[:,i]) for i in 1:size(nL0red)[2]]

figure("ntw20sup l0")
PyPlot.plot([ff1,ff1],[-.1,1.1],"--",color="C7")
PyPlot.plot([ff2,ff2],[-.1,1.1],"--",color="C7")
PyPlot.fill([1;ks;ks[end:-1:1];1]*2π/(N*τ),-[nL0max[1];nL0max;nL0min[end:-1:1];nL0min[1]],color=gra)
PyPlot.plot([1;ks]*2π/(N*τ),-[nL0[sol11,1];nL0[sol11,:]],color=cols[1])
PyPlot.plot([1;ks]*2π/(N*τ),-[nL0[sol12,1];nL0[sol12,:]],color=cols[2])
PyPlot.plot([1;ks]*2π/(N*τ),-[nL0[sol13,1];nL0[sol13,:]],color=cols[3])
PyPlot.plot([1;ks]*2π/(N*τ),-[nL0[sol14,1];nL0[sol14,:]],color=cols[4])
axis([ks[1]*2π/(N*τ),ks[end]*2π/(N*τ),-.1,1.1])
xlabel("freq")
ylabel("rescaled \n log-likelihood")

# =#



# #=
figure("ntw20sup l1")

L1 = Float64[]
γ1 = Vector{Vector{Float64}}()
for k in ks
	push!(L1,readdlm("data_melvyn/"*ntw*"/l1_missing_hidden_20_test__l1_$(k)_obj.csv",',')[1])
	push!(γ1,vec(readdlm("data_melvyn/"*ntw*"/l1_missing_hidden_20_test__l1_$(k)_g.csv",',')))
end
nL1 = (L1 .- maximum(L1))./(maximum(L1) - minimum(L1))

subplot(3,2,1)
PyPlot.plot([ff1,ff1],[-.1,1.1],"--",color="C7")
PyPlot.plot(ks*2π/(N*τ),-nL1,"k")
axis([ks[1]*2π/(N*τ),ks[end]*2π/(N*τ),-.1,1.1])
xlabel("freq")
ylabel("rescaled log-likelihood")

subplot(3,2,2)
PyPlot.plot([14,14],[-.1,1.1],"--",color="C7")
γ = γ1[111]
PyPlot.plot([1:12;17:20],γ[[1:12;16:19]]./maximum(γ),"ok")
PyPlot.plot(13,γ[13]/maximum(γ),"o",color=cols[1])
PyPlot.plot(15,γ[14]/maximum(γ),"o",color=cols[3])
PyPlot.plot(16,γ[15]/maximum(γ),"o",color=cols[2])
axis([.5,20.5,-.1,1.1])
xlabel("node index")
ylabel("rescaled amplitude")


L2 = Float64[]
γ2 = Vector{Vector{Float64}}()
for k in ks
	push!(L2,readdlm("data_melvyn/"*ntw*"/l1_2_hidden_20_test__l1_$(k)_obj.csv",',')[1])
	push!(γ2,vec(readdlm("data_melvyn/"*ntw*"/l1_2_hidden_20_test__l1_$(k)_g.csv",',')))
end
nL2 = (L2 .- maximum(L2))./(maximum(L2) - minimum(L2))

subplot(3,2,3)
PyPlot.plot([ff1,ff1],[-.1,1.1],"--",color=cols[1])
PyPlot.plot([ff2,ff2],[-.1,1.1],"--",color=cols[2])
PyPlot.plot(ks*2π/(N*τ),-nL2,"k")
axis([ks[1]*2π/(N*τ),ks[end]*2π/(N*τ),-.1,1.1])
xlabel("freq")
ylabel("rescaled log-likelihood")

subplot(3,2,4)
PyPlot.plot([14,14],[-.1,1.1],"--",color="C7")
PyPlot.plot([18,18],[-.1,1.1],"--",color="C7")
ids = [111,42]
for i in 1:length(ids)
	γ = γ2[ids[i]]
	PyPlot.plot(1:20,γ./maximum(γ),"o",color=cols[i])
end
axis([.5,20.5,-.1,1.1])
xlabel("node index")
ylabel("rescaled amplitude")


L3 = Float64[]
γ3 = Vector{Vector{Float64}}()
for k in ks
	push!(L3,readdlm("data_melvyn/"*ntw*"/l1_2_missing_hidden_20_test__l1_$(k)_obj.csv",',')[1])
	push!(γ3,vec(readdlm("data_melvyn/"*ntw*"/l1_2_missing_hidden_20_test__l1_$(k)_g.csv",',')))
end
nL3 = (L3 .- maximum(L3))./(maximum(L3) - minimum(L3))

subplot(3,2,5)
PyPlot.plot([ff1,ff1],[-.1,1.1],"--k")
PyPlot.plot([ff2,ff2],[-.1,1.1],"--",color=cols[4])
PyPlot.plot(ks*2π/(N*τ),-nL3,"k")
axis([ks[1]*2π/(N*τ),ks[end]*2π/(N*τ),-.1,1.1])
xlabel("freq")
ylabel("rescaled log-likelihood")

subplot(3,2,6)
PyPlot.plot([14,14],[-.1,1.1],"--",color="C7")
PyPlot.plot([18,18],[-.1,1.1],"--",color="C7")
γ = γ3[111]
PyPlot.plot([1:12;17:20],γ[[1:12;16:19]]./maximum(γ),"ok")
PyPlot.plot(13,γ[13]/maximum(γ),"o",color=cols[1])
PyPlot.plot(15,γ[14]/maximum(γ),"o",color=cols[3])
PyPlot.plot(16,γ[15]/maximum(γ),"o",color=cols[2])
γ = γ3[42]
PyPlot.plot([1:13;15:20],γ./maximum(γ),"o",color=cols[4])
axis([.5,20.5,-.1,1.1])
xlabel("node index")
ylabel("rescaled amplitude")

# =#


