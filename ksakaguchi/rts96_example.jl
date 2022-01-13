using PyPlot, PowerModels, DelimitedFiles, SparseArrays

include("plot_voltages.jl")
include("NR.jl")
include("../../DFNSolvers/cyclic_iterations.jl")

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



figure("rts 96",(5,10))

# #=
subplot(2,1,1)
xy = readdlm("ntw_data/rts96_xy.csv",',')
x = xy[:,1]
y = xy[:,2]

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
	PyPlot.plot(x[i],y[i],sh[bt[i]])
	PyPlot.text(x[i],y[i],"$i")
end

# =#

# Solving with our iteration scheme for KSakaguchi
h,hi,γ,B,as,ϕs = load_ksakaguchi(Y)
n,m = size(B)
Δ0 = [((γ[i][2]-γ[i][1])*rand() + γ[i][1]) for i in 1:m]
include("ntw_data/rts96_cycles.jl")
δ = .01
s = 1.
max_iter = 1000
tol = 1e-5

u1 = Int64.(vec(readdlm("ntw_data/rts96_w_u1.csv",',')))
Δ1,Δs1 = iterations(Δ0,B,C,u1,P0,h,γ,δ,s,max_iter,tol)

u2 = Int64.(vec(readdlm("ntw_data/rts96_w_u2.csv",',')))
Δ2,Δs2 = iterations(Δ0,B,C,u2,P0,h,γ,δ,s,max_iter,tol)

Lst = readdlm("ntw_data/rts96_spantree_L.csv",',')
Bst,w,Bstt = L2B(Lst)
id_st = Int64.(vec(readdlm("ntw_data/rts96_ids_spantree.csv",',')))

θ1 = Bstd'*Δ1[id_st]
θ2 = Bstd'*Δ2[id_st]


# Solving by Newton-Raphson, with appropriate initial conditions
t00 = zeros(n)
t01 = vec(readdlm("ntw_data/rts96_th1init.csv"))
v0 = ones(n)

v1,t1,i1 = NR_solver(Y,v0,t00,P0,Q0,bt)
v2,t2,i2 = NR_solver(Y,v0,t01,P0,Q0,bt)


#subplot(2,1,2)
figure()
shift = -2.7
PyPlot.plot([-1,1],[0,0],"k")
PyPlot.plot([0,0],[-1,1],"k")
plot_voltages(.6*exp.(im*(θ1 .- θ1[73] .+ t1[73])),cycle,"C1",.6,.6)
plot_voltages(.6*exp.(im*(θ2 .- θ2[39] .+ t2[39] .+ shift)),cycle,"C3",-1.)
plot_voltages(v1.*exp.(im*t1),cycle,col1)
plot_voltages(v2.*exp.(im*(t2 .+ shift)),cycle,col2,-1.)




