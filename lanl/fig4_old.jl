using PyPlot, DelimitedFiles

cmap = get_cmap("plasma")

@info "Plot pen_x."
@info "x = ? (2, 3, 4, 5, or 8)"

ntw = "pen_"*readline()

if ntw == "pen_2"
	n = 129
	ks = Array(5000:50:6000)
	T = 1260. 
	file = "pen_2"
elseif ntw == "pen_3"
	n = 130
	ks = Array(8500:50:10000)
	T = 660.
	file = "pen_3"
elseif ntw == "pen_4"
	n = 129
	ks = Array(6000:50:7000)
	T = 600.
	file = "pen_4"
elseif ntw == "pen_5"
	n = 130
	ks = Array(5000:50:6500)
	T = 1260.
	file = "pen_5"
elseif ntw == "pen_8"
	n = 134
	ks = Array(9000:50:10000)
	T = 1260.
	file = "pen_8"
else
	@info "Invalid entry."
end

ls = 1:n

L0 = zeros(length(ls),length(ks))
for i in 1:length(ls)
	for j in 1:length(ks)
		L0[i,j] = readdlm("data/"*file*"_l0_$(ls[i]).$(ks[j])_obj.csv",',')[1]
	end
end
L0mi = minimum(L0)

L1 = Array{Float64,1}()
γ1 = Array{Array{Float64,1},1}()

for j in 1:length(ks)
	push!(L1,readdlm("data/"*file*"_l1_$(ks[j])_obj.csv",',')[1])
	push!(γ1,vec(readdlm("data/"*file*"_l1_$(ks[j])_g.csv",',')))
end
L1mi,iii = findmin(L1)

figure("pen - utk",(19,4.5))

subplot2grid((2,4),(0,0),colspan=1,rowspan=2)
for i in 1:length(ls)
	PyPlot.plot(ks/T,L0[i,:],"-o",color=cmap((i-1)/(length(ls)-1)))
end
xlabel("freq")
ylabel("obj")

subplot2grid((2,4),(0,1),colspan=1,rowspan=1)
PyPlot.plot(ks/T,L1,"-o",color=cmap(.4))
xlabel("freq")
ylabel("obj")

subplot2grid((2,4),(1,1),colspan=1,rowspan=1)
PyPlot.plot(ls,γ1[iii],"-o",color=cmap(0.))
xlabel("node id")
ylabel("amplitude")


# #=
ntw = "utk"
n = 99
ls = 1:n
ks = 1:1500
T = 300.
file = "mysterious_forcing"
# =#

L0 = zeros(length(ls),length(ks))
for i in 1:length(ls)
	for j in 1:length(ks)
		L0[i,j] = readdlm("data/"*file*"_l0_$(ls[i]).$(ks[j])_obj.csv",',')[1]
	end
end
L0mi = minimum(L0)

L1 = Array{Float64,1}()
γ1 = Array{Array{Float64,1},1}()

for j in 1:length(ks)
	push!(L1,readdlm("data/"*file*"_l1_$(ks[j])_obj.csv",',')[1])
	push!(γ1,vec(readdlm("data/"*file*"_l1_$(ks[j])_g.csv",',')))
end
L1mi,iii = findmin(L1)

subplot2grid((2,4),(0,2),colspan=1,rowspan=2)
for i in 1:length(ls)
	PyPlot.plot(ks/T,L0[i,:],"-o",color=cmap((i-1)/(length(ls)-1)))
end
xlabel("freq")
ylabel("obj")

subplot2grid((2,4),(0,3),colspan=1,rowspan=1)
PyPlot.plot(ks/T,L1,"-o",color=cmap(.4))
xlabel("freq")
ylabel("obj")

subplot2grid((2,4),(1,3),colspan=1,rowspan=1)
PyPlot.plot(ls,γ1[iii],"-o",color=cmap(0.))
xlabel("node id")
ylabel("amplitude")

