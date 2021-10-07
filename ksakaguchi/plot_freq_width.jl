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

run = 6179 #run = 3695
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
#PyPlot.fill(αs[[id;di]],[β2[1,id];β2[3,di]],color=cmap(0.),alpha=.3,label=runs[run]["labs"][2])
PyPlot.plot(αs[id],β2[1,id],"-",color=cmap(0.),linewidth=lw,label=runs[run]["labs"][2])
#PyPlot.plot(αs[id],β2[2,id],"--",color=cmap(0.))
PyPlot.plot(αs[id],max.(β2[size(β2)[1],id],zeros(T)),"--",color=cmap(0.))
#PyPlot.fill(αs[[id;di]],[β0[1,id];β0[3,di]],color=cmap(1/3),alpha=.3,label=runs[run]["labs"][0])
PyPlot.plot(αs[id],β0[1,id],"-",color=cmap(1/3),linewidth=lw,label=runs[run]["labs"][0])
#PyPlot.plot(αs[id],β0[2,id],"--",color=cmap(1/3))
PyPlot.plot(αs[id],max.(β0[size(β0)[1],id],zeros(T)),"--",color=cmap(1/3))
#PyPlot.fill(αs[[id;di]],[β1[1,id];β1[3,di]],color=cmap(2/3),alpha=.3,label=runs[run]["labs"][1])
PyPlot.plot(αs[id],β1[1,id],"-",color=cmap(2/3),linewidth=lw,label=runs[run]["labs"][1])
#PyPlot.plot(αs[id],β1[2,id],"--",color=cmap(2/3))
PyPlot.plot(αs[id],max.(β1[size(β1)[1],id],zeros(T)),"--",color=cmap(2/3))
title(runs[run]["tit"]*" (run: $run)")
xlabel("ϕ")
ylabel("β")
legend()

subplot(1,2,2)
w = vec(readdlm("temp_data/w_$(run).csv",','))
plot_ntws(runs[run]["ntw"],w,10.,2.,name_cmap)
axis([-1.3,1.3,-1.3,1.3])
