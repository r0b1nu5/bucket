using PyPlot, DelimitedFiles

cmap = get_cmap("plasma")

ls = Array(1:99)
ls2 = Array([1:84;86:99])
ks = Array(1:1500)

L = zeros(99,1500)
for i in 1:99
	for j in 1:1500
		L[i,j] = readdlm("data/mysterious_forcing_l0_$(i).$(j)_obj.csv",',')[1]
	end
end
L2 = L[ls2,:]

Lmi,iii = findmin(L)
Lmi2,iii2 = findmin(L2)

node_id = iii[1]
node_id2 = iii2[1]

freq = round(ks[iii[2]]/300.,digits=3)
freq2 = round(ks[iii2[2]]/300.,digits=3)

sort_nodes = sortslices([L[:,iii[2]] ls Array(1:99)],dims=1)
sort_nodes2 = sortslices([L2[:,iii2[2]] ls2 Array(1:98)],dims=1)

figure("...",(18.,10.))

subplot2grid((11,14),(0,0),colspan=5,rowspan=5)
for i in 1:length(ls)
	j = Int(sort_nodes[i,3])
	PyPlot.plot(ks/300.,L[j,:],"-o",color=cmap((i-1)/(length(ls)-1)))
end
xlabel("freq")
ylabel("obj")
twiny()
PyPlot.plot([ks[1],ks[end]],[Lmi,Lmi],"--k")
xlabel("k")
title("Node id: $(node_id) \n freq: $freq [Hz]")

subplot2grid((11,14),(0,5),colspan=1,rowspan=5)
xticks([])
yticks([])
twinx()
yticks(-ls,sort_nodes[:,2])
for i in 1:length(ls)
	PyPlot.plot([0,1],[-i,-i],color=cmap((i-1)/(length(ls)-1)),linewidth=2)
end
axis([-.1,1.1,-maximum(ls)-1,-minimum(ls)+1])

subplot2grid((11,14),(0,7),colspan=5,rowspan=5)
for i in 1:length(ls2)
	j = Int(sort_nodes2[i,3])
	PyPlot.plot(ks/300.,L2[j,:],"-o",color=cmap((i-1)/(length(ls2)-1)))
end
xlabel("freq")
ylabel("obj")
twiny()
PyPlot.plot([ks[1],ks[end]],[Lmi2,Lmi2],"--k")
xlabel("k")
title("Node id: $(node_id2) \n freq: $freq2 [Hz]")

subplot2grid((11,14),(0,12),colspan=1,rowspan=5)
xticks([])
yticks([])
twinx()
yticks(-Array(1:length(ls2)),sort_nodes2[:,2])
for i in 1:length(ls2)
	PyPlot.plot([0,1],[-i,-i],color=cmap((i-1)/(length(ls2)-1)),linewidth=2)
end
axis([-.1,1.1,-maximum(ls2)-1,-minimum(ls2)+1])


### l1

L = Array{Float64,1}()
γ = Array{Array{Float64,1},1}()
for j in 1:1500
	push!(L,readdlm("data/mysterious_forcing_l1_$(j)_obj.csv",',')[1])
	push!(γ,vec(readdlm("data/mysterious_forcing_l1_$(j)_g.csv",',')))
end

Lmi,iii = findmin(L)
Lma = maximum(L)
nma,jjj = findmax(γ[iii])

node_id = jjj
freq = round(ks[iii]/300,digits=3)

subplot2grid((11,15),(6,0),colspan=5,rowspan=5)
PyPlot.plot(ks/300,L,"-o",color=cmap(.4))
xlabel("freq")
ylabel("obj")
twiny()
PyPlot.plot([ks[1],ks[end]],[Lmi,Lmi],"--k")
xlabel("k")

subplot2grid((11,14),(6,7),colspan=5,rowspan=5)
PyPlot.plot(1:length(γ[iii]),γ[iii],"-o",color=cmap(0.),label="Min obj, k = $(ks[iii])")
legend()
xlabel("node idx")
ylabel("amplitude")
title("Node id: $(node_id) \n freq: $freq [Hz]")

