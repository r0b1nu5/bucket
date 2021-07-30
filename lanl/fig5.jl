using PyPlot, DelimitedFiles

figure("fig5")
cmap = get_cmap("plasma")

f = .3
t = LinRange(0,10,200)

# Multi-sine
subplot2grid((3,38),(0,0),colspan=9)
PyPlot.plot(t,sin.(2π*f*t) + sin.(4π*f*t) + sin.(6π*f*t))

# l0
ls = Array(1:20)
ks = Array(10:10:200)
L = zeros(20,20)
for i in 1:20
	for j in 1:20
		L[i,j] = readdlm("data/ntw20_multisine_l0_$(ls[i]).$(ks[j])_obj.csv",',')[1]
	end
end
Lmi,iii = findmin(L)
freq = round(ks[iii[2]]/100.,digits=3)
sort_nodes = sortslices([L[:,iii[2]] Array(1:20) Array(1:length(ls))],dims=1)

subplot2grid((3,38),(0,9),colspan=9)
for i in 1:length(ls)
	j = Int.(sort_nodes[i,3])
	PyPlot.plot(ks/100.,L[j,:],"-",color=cmap((i-1)/(length(ls)-1)))
end
xlabel("freq")
ylabel("obj")
twiny()
PyPlot.plot([ks[1],ks[end]],[Lmi,Lmi],"--k")
xlabel("k")

#l1
L = Array{Float64,1}()
γ = Array{Array{Float64,1},1}()
for j in 1:20
	push!(L,readdlm("data/ntw20_multisine_l1_$(ks[j])_obj.csv",',')[1])
	push!(γ,vec(readdlm("data/ntw20_multisine_l1_$(ks[j])_g.csv",',')))
end
Lmi,iii = findmin(L)
Lma = maximum(L)
nma,jjj = findmax(γ[iii])
freq = round(ks[iii]/100.,digits=3)

subplot2grid((3,38),(0,18),colspan=9)
PyPlot.plot(ks/100.,L,"-",color=cmap(.4))
xlabel("freq")
ylabel("obj")
twiny()
PyPlot.plot([ks[1],ks[end]],[Lmi,Lmi],"--k")
xlabel("k")

subplot2grid((3,38),(0,27),colspan=9)
#for j in 1:length(ks)
#	α = cmap(1 - (L[j] - Lma)/(Lmi - Lma))
#	PyPlot.plot(1:length(γ[i]),γ[i],color=α)
#end
#PyPlot.plot(jjj,γ[iii][jjj],"or",makersize=10)
#PyPlot.plot(jjj,γ[iii][jjj],"ow",makersize=7)
PyPlot.plot(1:length(γ[iii]),γ[iii],"-",color=cmap(0.))
xlabel("idx")
ylabel("ampl")

#=
subplot2grid((3,38),(0,36),colspan=1)
yticks([])
twinx()
for ll in LinRange(Lmi,Lma,1000)
	α = cmap(1 - (ll - Lma)/(Lmi - Lma))
	PyPlot.plot([0,1],[ll,ll],color=α)
end
axis([0,1,Lmi,Lma])
xticks([])
ylabel("obj")
=#




# Saw
subplot2grid((3,38),(1,0),colspan=9)
PyPlot.plot(t,2*mod.(f*t,1.) .- 1)

# l0
ls = Array(1:20)
ks = Array(10:10:200)
L = zeros(20,20)
for i in 1:20
	for j in 1:20
		L[i,j] = readdlm("data/ntw20_saw_l0_$(ls[i]).$(ks[j])_obj.csv",',')[1]
	end
end
Lmi,iii = findmin(L)
freq = round(ks[iii[2]]/100.,digits=3)
sort_nodes = sortslices([L[:,iii[2]] Array(1:20) Array(1:length(ls))],dims=1)

subplot2grid((3,38),(1,9),colspan=9)
for i in 1:length(ls)
	j = Int.(sort_nodes[i,3])
	PyPlot.plot(ks/100.,L[j,:],"-",color=cmap((i-1)/(length(ls)-1)))
end
xlabel("freq")
ylabel("obj")
twiny()
PyPlot.plot([ks[1],ks[end]],[Lmi,Lmi],"--k")
xlabel("k")

#l1
L = Array{Float64,1}()
γ = Array{Array{Float64,1},1}()
for j in 1:20
	push!(L,readdlm("data/ntw20_saw_l1_$(ks[j])_obj.csv",',')[1])
	push!(γ,vec(readdlm("data/ntw20_saw_l1_$(ks[j])_g.csv",',')))
end
Lmi,iii = findmin(L)
Lma = maximum(L)
nma,jjj = findmax(γ[iii])
freq = round(ks[iii]/100.,digits=3)

subplot2grid((3,38),(1,18),colspan=9)
PyPlot.plot(ks/100.,L,"-",color=cmap(.4))
xlabel("freq")
ylabel("obj")
twiny()
PyPlot.plot([ks[1],ks[end]],[Lmi,Lmi],"--k")
xlabel("k")

subplot2grid((3,38),(1,27),colspan=9)
#for j in 1:length(ks)
#	α = cmap(1 - (L[j] - Lma)/(Lmi - Lma))
#	PyPlot.plot(1:length(γ[i]),γ[i],color=α)
#end
#PyPlot.plot(jjj,γ[iii][jjj],"or",makersize=10)
#PyPlot.plot(jjj,γ[iii][jjj],"ow",makersize=7)
PyPlot.plot(1:length(γ[iii]),γ[iii],"-",color=cmap(0.))
xlabel("idx")
ylabel("ampl")

#=
subplot2grid((3,38),(1,36),colspan=1)
yticks([])
twinx()
for ll in LinRange(Lmi,Lma,1000)
	α = cmap(1 - (ll - Lma)/(Lmi - Lma))
	PyPlot.plot([0,1],[ll,ll],color=α)
end
axis([0,1,Lmi,Lma])
xticks([])
ylabel("obj")
=#


# Step
subplot2grid((3,38),(2,0),colspan=9)
PyPlot.plot(t,2*mod.(floor.(f*t*2),2).-1)

# l0
ls = Array(1:20)
ks = Array(10:10:200)
L = zeros(20,20)
for i in 1:20
	for j in 1:20
		L[i,j] = readdlm("data/ntw20_step_l0_$(ls[i]).$(ks[j])_obj.csv",',')[1]
	end
end
Lmi,iii = findmin(L)
freq = round(ks[iii[2]]/100.,digits=3)
sort_nodes = sortslices([L[:,iii[2]] Array(1:20) Array(1:length(ls))],dims=1)

subplot2grid((3,38),(2,9),colspan=9)
for i in 1:length(ls)
	j = Int.(sort_nodes[i,3])
	PyPlot.plot(ks/100.,L[j,:],"-",color=cmap((i-1)/(length(ls)-1)))
end
xlabel("freq")
ylabel("obj")
twiny()
PyPlot.plot([ks[1],ks[end]],[Lmi,Lmi],"--k")
xlabel("k")

#l1
L = Array{Float64,1}()
γ = Array{Array{Float64,1},1}()
for j in 1:20
	push!(L,readdlm("data/ntw20_step_l1_$(ks[j])_obj.csv",',')[1])
	push!(γ,vec(readdlm("data/ntw20_step_l1_$(ks[j])_g.csv",',')))
end
Lmi,iii = findmin(L)
Lma = maximum(L)
nma,jjj = findmax(γ[iii])
freq = round(ks[iii]/100.,digits=3)

subplot2grid((3,38),(2,18),colspan=9)
PyPlot.plot(ks/100.,L,"-",color=cmap(.4))
xlabel("freq")
ylabel("obj")
twiny()
PyPlot.plot([ks[1],ks[end]],[Lmi,Lmi],"--k")
xlabel("k")

subplot2grid((3,38),(2,27),colspan=9)
#for j in 1:length(ks)
#	α = cmap(1 - (L[j] - Lma)/(Lmi - Lma))
#	PyPlot.plot(1:length(γ[i]),γ[i],color=α)
#end
#PyPlot.plot(jjj,γ[iii][jjj],"or",makersize=10)
#PyPlot.plot(jjj,γ[iii][jjj],"ow",makersize=7)
PyPlot.plot(1:length(γ[iii]),γ[iii],"-",color=cmap(0.))
xlabel("idx")
ylabel("ampl")

#=
subplot2grid((3,38),(2,36),colspan=1)
yticks([])
twinx()
for ll in LinRange(Lmi,Lma,1000)
	α = cmap(1 - (ll - Lma)/(Lmi - Lma))
	PyPlot.plot([0,1],[ll,ll],color=α)
end
axis([0,1,Lmi,Lma])
xticks([])
ylabel("obj")
=#





