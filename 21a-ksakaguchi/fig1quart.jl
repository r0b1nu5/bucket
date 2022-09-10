using PyPlot, DelimitedFiles

include("plot_ntws.jl")
include("runs_dict.jl")

name_cmap = "viridis"
cmap = get_cmap(name_cmap)
# #=
col1 = cmap(1/6)
col2 = cmap(3/6)
col3 = cmap(5/6)
# =#

run = 6179
run = 6786

lw = 2.

αs = vec(readdlm("temp_data/alphas_$(run).csv",','))
β0 = readdlm("temp_data/beta0_$(run).csv",',')
β1 = readdlm("temp_data/beta1_$(run).csv",',')
β2 = readdlm("temp_data/beta2_$(run).csv",',')

T = length(αs)
T = 70
id = Array(1:T)
di = Array(T:-1:1)

figure("1cyc")
subplot(1,3,1)
PyPlot.plot(cos.(2π*(1:7)/18),sin.(2π*(1:7)/18),color=col1,linewidth=lw)
PyPlot.plot(cos.(2π*(7:19)/18),sin.(2π*(7:19)/18),color=col2,linewidth=lw)
PyPlot.plot(cos.(2π*(1:18)/18),sin.(2π*(1:18)/18),"ok",markersize=8)
xticks([])
yticks([])

subplot(1,3,2)
θ1 = vec(readdlm("temp_data/th0f_6179.csv",','))
PyPlot.plot(cos.(LinRange(0,2π,200)),sin.(LinRange(0,2π,200)),"k")
for i in 1:18
	PyPlot.plot([.9,1.1]*cos.(θ1[i]),[.9,1.1]*sin.(θ1[i]),"-k",linewidth=5)
end
PyPlot.plot(.9*cos.(θ1[1:7]),.9*sin.(θ1[1:7]),"-o",color=col1,linewidth=lw)
PyPlot.plot(1.1*cos.([θ1[7:18];θ1[1]]),1.1*sin.([θ1[7:18];θ1[1]]),"-o",color=col2,linewidth=lw)
title("q = 0")

subplot(1,3,3)
θ2 = vec(readdlm("temp_data/th2f_6179.csv",','))
PyPlot.plot(cos.(LinRange(0,2π,200)),sin.(LinRange(0,2π,200)),"k")
for i in 1:18
	PyPlot.plot([.9,1.1]*cos.(θ2[i]),[.9,1.1]*sin.(θ2[i]),"-k",linewidth=5)
end
PyPlot.plot(1.1*cos.(θ2[1:7]),1.1*sin.(θ2[1:7]),"-o",color=col1,linewidth=lw)
PyPlot.plot([cos.(θ2[7:10]);.9*cos.(θ2[[11:18;1]])],[sin.(θ2[7:10]);.9*sin.(θ2[[11:18;1]])],"-o",color=col2,linewidth=lw)
title("q = -1")

figure("2cyc")
subplot(1,3,1)
PyPlot.plot(cos.(2π*(1:10)/18),sin.(2π*(1:10)/18),color=col1,linewidth=lw)
PyPlot.plot(cos.(2π*(10:19)/18),sin.(2π*(10:19)/18),color=col2,linewidth=lw)
PyPlot.plot(cos.(2π*[1,10]/18),sin.(2π*[1,10]/18),color=col3,linewidth=lw)
PyPlot.plot(cos.(2π*(1:18)/18),sin.(2π*(1:18)/18),"ok",markersize=8)
xticks([])
yticks([])

subplot(1,3,2)
θ3 = vec(readdlm("temp_data/th0f_6786.csv",','))
PyPlot.plot(cos.(LinRange(0,2π,200)),sin.(LinRange(0,2π,200)),"k")
for i in 1:18
	PyPlot.plot([.9,1.1]*cos.(θ3[i]),[.9,1.1]*sin.(θ3[i]),"-k",linewidth=5)
end
PyPlot.plot(1.1*cos.(θ3[1:10]),1.1*sin.(θ3[1:10]),"-o",color=col1,linewidth=lw)
PyPlot.plot(cos.([θ3[10:18];θ3[1]]),sin.([θ3[10:18];θ3[1]]),"-o",color=col2,linewidth=lw)
PyPlot.plot(.9*cos.(θ3[[1,10]]),.9*sin.(θ3[[1,10]]),"-o",color=col3,linewidth=lw)
title("q = [0,0]")

subplot(1,3,3)
θ4 = vec(readdlm("temp_data/th1f_6786.csv",','))
PyPlot.plot(cos.(LinRange(0,2π,200)),sin.(LinRange(0,2π,200)),"k")
for i in 1:18
	PyPlot.plot([.9,1.1]*cos.(θ4[i]),[.9,1.1]*sin.(θ4[i]),"-k",linewidth=5)
end
PyPlot.plot(1.1*cos.(θ4[1:10]),1.1*sin.(θ4[1:10]),"-o",color=col1,linewidth=lw)
PyPlot.plot(cos.([θ4[10:18];θ4[1]]),sin.([θ4[10:18];θ4[1]]),"-o",color=col2,linewidth=lw)
PyPlot.plot(.9*cos.(θ4[[1,10]]),.9*sin.(θ4[[1,10]]),"-o",color=col3,linewidth=lw)
title("q = [-1,1]")



