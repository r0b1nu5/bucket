using PyPlot, DelimitedFiles

include("plot_ntws.jl")
include("runs_dict.jl")

cmap = get_cmap("plasma")

run = 3538

αs = vec(readdlm("temp_data/alphas_$(run).csv",','))
β0 = readdlm("temp_data/beta0_$(run).csv",',')
β1 = readdlm("temp_data/beta1_$(run).csv",',')
β2 = readdlm("temp_data/beta2_$(run).csv",',')

T = length(αs)
T = 70
id = Array(1:T)
di = Array(T:-1:1)

figure()
subplot(1,2,1)
PyPlot.fill(αs[[id;di]],[β2[1,id];β2[3,di]],color=cmap(0.),alpha=.3,label=runs[run]["labs"][2])
PyPlot.plot(αs[id],β2[1,id],"-o",color=cmap(0.))
PyPlot.plot(αs[id],β2[2,id],"-",color=cmap(0.))
PyPlot.fill(αs[[id;di]],[β0[1,id];β0[3,di]],color=cmap(1/3),alpha=.3,label=runs[run]["labs"][0])
PyPlot.plot(αs[id],β0[1,id],"-o",color=cmap(1/3))
PyPlot.plot(αs[id],β0[2,id],"-",color=cmap(1/3))
PyPlot.fill(αs[[id;di]],[β1[1,id];β1[3,di]],color=cmap(2/3),alpha=.3,label=runs[run]["labs"][1])
PyPlot.plot(αs[id],β1[1,id],"-o",color=cmap(2/3))
PyPlot.plot(αs[id],β1[2,id],"-",color=cmap(2/3))
title(runs[run]["tit"]*" (run: $run)")
xlabel("α")
ylabel("width")
legend()

subplot(1,2,2)
w = vec(readdlm("temp_data/w_$(run).csv",','))
plot_ntws(runs[run]["ntw"],w)

