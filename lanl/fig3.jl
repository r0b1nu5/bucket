using PyPlot, DelimitedFiles

# #=
ntw = "uk"
n = 120
ls = 1:n
ks0 = 1:10
ks1 = 31:120
T = 50000*.01
file = "mysterious_forcing_UK"
# =#

 #=
ntw = "ieee57"
n = 57
ls = 1:n
ks0 = 1:10
ks1 = 1:57
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

figure(ntw,(9.5,4.5))

subplot2grid((2,2),(0,0),colspan=1,rowspan=2)
for i in 1:length(ls)
	PyPlot.plot(ks0/T,L0[i,:],"-o")
end
xlabel("freq")
ylabel("obj")

subplot2grid((2,2),(0,1),colspan=1,rowspan=1)
PyPlot.plot(ks1/T,L1,"-o")
xlabel("freq")
ylabel("obj")

subplot2grid((2,2),(1,1),colspan=1,rowspan=1)
PyPlot.plot(ls,γ1[iii],"-o")
xlabel("node id")
ylabel("amplitude")




