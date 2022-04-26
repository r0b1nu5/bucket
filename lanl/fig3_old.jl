using PyPlot, DelimitedFiles

cmap = get_cmap("plasma")

inset = true

# #=
ntw = "ieee57"
n = 57
ls = 1:n
ks0 = 1:50
ks1 = 1:50
T = 200000*2e-3
file = "mysterious_forcing_57"
# =#

L0 = zeros(length(ls),length(ks0))
for i in 1:length(ls)
	for j in 1:length(ks0)
		L0[i,j] = readdlm("data/"*file*"_l0_$(ls[i]).$(ks0[j])_obj.csv",',')[1]
	end
end
L0mi = minimum(L0)

L1 = Array{Float64,1}()
γ1 = Array{Array{Float64,1},1}()

for j in 1:length(ks1)
	push!(L1,readdlm("data/"*file*"_l1_$(ks1[j])_obj.csv",',')[1])
	push!(γ1,vec(readdlm("data/"*file*"_l1_$(ks1[j])_g.csv",',')))
end
L1mi,iii = findmin(L1)

figure("ieee57_uk",(19,4.5))

subplot2grid((2,4),(0,0),colspan=1,rowspan=2)
for i in 1:length(ls)
	PyPlot.plot(ks0/T,L0[i,:],"-o",color=cmap((i-1)/(length(ls)-1)))
end
PyPlot.text(.12,-.00003401,"(a)")
xlabel("freq")
ylabel("obj")

if inset
	xy0 = readdlm("data_melvyn/ieee57_xy.csv",',')
	adj = Int.(readdlm("data_melvyn/ieee57_adj.csv",','))
	xmi = .038
	xma = .1145
	ymi = -3.43e-5
	yma = -3.407e-5
	xy = [(xy0[:,1]*(xma-xmi) .+ xmi) (xy0[:,2]*(yma-ymi) .+ ymi)]
	for i in 1:2:size(adj)[1]
		PyPlot.plot(xy[adj[i,1:2],1],xy[adj[i,1:2],2],"k",linewidth=1.)
	end
	PyPlot.plot(xy[:,1],xy[:,2],"ok",markersize=5.)
end	

subplot2grid((2,4),(0,1),colspan=1,rowspan=1)
PyPlot.plot(ks1/T,L1,"-o",color=cmap(.4))
PyPlot.text(.12,-.0000340,"(b)")
xlabel("freq")
ylabel("obj")

subplot2grid((2,4),(1,1),colspan=1,rowspan=1)
PyPlot.plot(ls,γ1[iii],"-o",color=cmap(0.))
PyPlot.text(55.,.0017,"(c)")
xlabel("node id")
ylabel("amplitude")


# #=
ntw = "uk"
n = 120
ls = 1:n
ks0 = 1:50
ks1 = 1:50
T = 50000*.01
file = "mysterious_forcing_UK"
# =#

L0 = zeros(length(ls),length(ks0))
for i in 1:length(ls)
	for j in 1:length(ks0)
		L0[i,j] = readdlm("data/"*file*"_l0_$(ls[i]).$(ks0[j])_obj.csv",',')[1]
	end
end
L0mi = minimum(L0)

L1 = Array{Float64,1}()
γ1 = Array{Array{Float64,1},1}()

for j in 1:length(ks1)
	push!(L1,readdlm("data/"*file*"_l1_$(ks1[j])_obj.csv",',')[1])
	push!(γ1,vec(readdlm("data/"*file*"_l1_$(ks1[j])_g.csv",',')))
end
L1mi,iii = findmin(L1)

subplot2grid((2,4),(0,2),colspan=1,rowspan=2)
for i in 1:length(ls)
	PyPlot.plot(ks0/T,L0[i,:],"-o",color=cmap((i-1)/(length(ls)-1)))
end
PyPlot.text(.095,-2.1091,"(d)")
xlabel("freq")
ylabel("obj")

if inset
	xy0 = readdlm("data_melvyn/uk_xy.csv",',')
	xyb = readdlm("data_melvyn/uk_bord.csv",',')
	adj = Int.(readdlm("data_melvyn/uk_adj.csv",','))
	xmi = .01
	xma = .1
	ymi = -2.1106
	yma = -2.1093
	xy = [(xyb[:,1]*(xma-xmi) .+ xmi) (xyb[:,2]*(yma-ymi) .+ ymi)]
	PyPlot.plot(xy[:,1],xy[:,2],"k",linewidth=.5)
	xy = [(xy0[:,1]*(xma-xmi) .+ xmi) (xy0[:,2]*(yma-ymi) .+ ymi)]
	for i in 1:2:size(adj)[1]
		PyPlot.plot(xy[adj[i,1:2],1],xy[adj[i,1:2],2],"k",linewidth=1.)
	end
	PyPlot.plot(xy[:,1],xy[:,2],"ok",markersize=4.)
end

subplot2grid((2,4),(0,3),colspan=1,rowspan=1)
PyPlot.plot(ks1/T,L1,"-o",color=cmap(.4))
PyPlot.text(.095,-.3091,"(e)")
xlabel("freq")
ylabel("obj")

subplot2grid((2,4),(1,3),colspan=1,rowspan=1)
PyPlot.plot(ls,γ1[iii],"-o",color=cmap(0.))
PyPlot.text(118.,.09,"(f)")
xlabel("node id")
ylabel("amplitude")

