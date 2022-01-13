using PyPlot, PowerModels, DelimitedFiles, SparseArrays

include("plot_voltages.jl")
include("NR.jl")

name_cmap = "viridis"
cmap = get_cmap(name_cmap)
# #=
col1 = cmap(1/6)
col2 = cmap(2.5/6)
col3 = cmap(4/6)
# =#

cycle = [13,39,40,43,44,47,66,65,64,67,68,71,73,21,15,16,19,20,23,13]
c = length(cycle)

Y = readdlm("ntw_data/rts96_G2.csv",',') + im*readdlm("ntw_data/rts96_B2.csv",',')

n = size(Y)[1]

θ0 = zeros(n)
θ1 = vec(readdlm("ntw_data/rts96_th1init.csv"))
v0 = ones(n)

P0 = vec(readdlm("ntw_data/rts96_P2.csv",','))
Q0 = vec(readdlm("ntw_data/rts96_Q2.csv",','))

bt = Int64.(vec(readdlm("ntw_data/rts96_bt.csv",',')))

# #=
# Adapting parameters...
l = [47,66]
Y[l,l] *= 1.3

Q0[73] -= 5.

l = [21,73]
Y[l,l] *= .9

l = [71,73]
Y[l,l] *= 1.1

l = [13,39]
Y[l,l] *= .7
l = [23,41]
Y[l,l] *= .8
l = [7,8]
Y[l,l] *= .8

P0[13] += 10.2
P0[39] -= 10.2
# =#

v1,t1,i1 = NR_solver(Y,v0,θ0,P0,Q0,bt)
v2,t2,i2 = NR_solver(Y,v0,θ1,P0,Q0,bt)

figure("rts 96",(5,10))

subplot(2,1,2)
plot_voltages(v1.*exp.(im*t1),cycle,col1)
plot_voltages(v2.*exp.(im*(t2 .- π/4)),cycle,col2)

# #=
subplot(2,1,1)
xy = readdlm("ntw_data/rts96_xy.csv",',')
x = xy[:,1]
y = xy[:,2]
I,J,V = findnz(calc_admittance_matrix(nd).matrix)
m = length(I)

for i in 1:n-1
	for j in i+1:n
		if norm(Y[i,j]) > 1e-2
			PyPlot.plot(x[[i,j]],y[[i,j]],"k")
		end
	end
end
for k in 2:c
	i = cycle[k-1]
	j = cycle[k]
	PyPlot.plot(x[[i,j]],y[[i,j]],color=col3,linewidth=2)
end

sh = Dict{Int64,String}(1 => "ok", 2 => "sk", 3 => "or")
for i in 1:n
	PyPlot.plot(x[i],y[i],sh[ty[i]])
	PyPlot.text(x[i],y[i],"$i")
end

# =#



