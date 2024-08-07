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
 #=
col1 = "C0"
col2 = "C1"
col3 = "C2"
# =#

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

figure("$run",(11.,5.))
subplot(1,2,1)
PyPlot.fill([αs[id];αs[di]],[β1[1,id];β1[2,di]]/sqrt(2),"k",alpha=.1)
PyPlot.plot(αs[id],β1[1,id]/sqrt(2),"-",color=col1,linewidth=lw,label=runs[run]["labs"][1])
PyPlot.plot(αs[id],β1[2,id]/sqrt(2),"--",color=col1,linewidth=lw)
PyPlot.fill([αs[id];αs[di]],[β0[1,id];-β0[1,di]]/sqrt(2),"k",alpha=.1)
PyPlot.plot(αs[id],β0[1,id]/sqrt(2),"-",color=col2,linewidth=lw,label=runs[run]["labs"][0])
PyPlot.plot(αs[id],-β0[1,id]/sqrt(2),"--",color=col2,linewidth=lw)
PyPlot.fill([αs[id];αs[di]],[-β1[2,id];-β1[1,di]]/sqrt(2),"k",alpha=.1)
PyPlot.plot(αs[id],-β1[2,id]/sqrt(2),"-",color=col3,linewidth=lw,label=runs[run]["labs"][2])
PyPlot.plot(αs[id],-β1[1,id]/sqrt(2),"--",color=col3,linewidth=lw)
title(runs[run]["tit"]*" (run: $run)")
xlabel("ϕ")
ylabel("β")
legend()

subplot(1,2,2)
θ = vec(readdlm("temp_data/th0f_6786.csv",','))
#θ = vec(readdlm("temp_data/th1f_6786.csv",',')) .+ .5
#plot_ntws(runs[run]["ntw"],w,10.,2.,name_cmap)
plot_ntws(runs[run]["ntw"],w,θ,10.,2.,name_cmap)
axis([-1.3,1.3,-1.3,1.3])

#flows0 = readdlm("temp_data/f0_$(run).csv",',')
#flows1 = readdlm("temp_data/f1_$(run).csv",',')
#flows2 = readdlm("temp_data/f2_$(run).csv2,',')


