using PyPlot, DelimitedFiles

cmap = get_cmap("plasma")

ntw1 = "pen_3"
ntw2 = "pen_4"

ntws = Dict{String,Any}("pen_2" => (129,Array(5000:50:6000),1260.,1000.),
			"pen_3" => (130,Array(8500:50:10000),660.,-.024),
			"pen_4" => (129,Array(6000:50:7000),600.,-.03),
			"pen_5" => (130,Array(5000:50:6500),1260.,-.015),
			"pen_8" => (134,Array(9000:50:10000),1260.,1000.)
			)

n1,ks1,T1,omax1 = ntws[ntw1]
n2,ks2,T2,omax2 = ntws[ntw2]
ls1 = 1:n1
ls2 = 1:n2


L0 = zeros(length(ls1),length(ks1))
for i in 1:length(ls1)
	for j in 1:length(ks1)
		L0[i,j] = readdlm("data/"*ntw1*"_l0_$(ls1[i]).$(ks1[j])_obj.csv",',')[1]
	end
end
L0mi = minimum(L0)

L1 = Array{Float64,1}()
γ1 = Array{Array{Float64,1},1}()

for j in 1:length(ks1)
	push!(L1,readdlm("data/"*ntw1*"_l1_$(ks1[j])_obj.csv",',')[1])
	push!(γ1,vec(readdlm("data/"*ntw1*"_l1_$(ks1[j])_g.csv",',')))
end
L1mi,iii = findmin(L1)

figure(ntw1*" - "*ntw2,(19,4.5))

subplot2grid((2,4),(0,0),colspan=1,rowspan=2)
for i in 1:length(ls1)
	if minimum(L0[i,:]) < omax1
		PyPlot.plot(ks1/T1,L0[i,:],"-o",color=cmap((i-1)/(length(ls1)-1)))
	end
end
xlabel("freq")
ylabel("obj")

subplot2grid((2,4),(0,1),colspan=1,rowspan=1)
PyPlot.plot(ks1/T1,L1,"-o",color=cmap(.4))
xlabel("freq")
ylabel("obj")

subplot2grid((2,4),(1,1),colspan=1,rowspan=1)
PyPlot.plot(ls1,γ1[iii],"-o",color=cmap(0.))
xlabel("node id")
ylabel("amplitude")


L0 = zeros(length(ls2),length(ks2))
for i in 1:length(ls2)
	for j in 1:length(ks2)
		L0[i,j] = readdlm("data/"*ntw2*"_l0_$(ls2[i]).$(ks2[j])_obj.csv",',')[1]
	end
end
L0mi = minimum(L0)

L1 = Array{Float64,1}()
γ1 = Array{Array{Float64,1},1}()

for j in 1:length(ks2)
	push!(L1,readdlm("data/"*ntw2*"_l1_$(ks2[j])_obj.csv",',')[1])
	push!(γ1,vec(readdlm("data/"*ntw2*"_l1_$(ks2[j])_g.csv",',')))
end
L1mi,iii = findmin(L1)

subplot2grid((2,4),(0,2),colspan=1,rowspan=2)
for i in 1:length(ls2)
	if minimum(L0[i,:]) < omax2
		PyPlot.plot(ks2/T2,L0[i,:],"-o",color=cmap((i-1)/(length(ls2)-1)))
	end
end
xlabel("freq")
ylabel("obj")

subplot2grid((2,4),(0,3),colspan=1,rowspan=1)
PyPlot.plot(ks2/T2,L1,"-o",color=cmap(.4))
xlabel("freq")
ylabel("obj")

subplot2grid((2,4),(1,3),colspan=1,rowspan=1)
PyPlot.plot(ls2,γ1[iii],"-o",color=cmap(0.))
xlabel("node id")
ylabel("amplitude")

