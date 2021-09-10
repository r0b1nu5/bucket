using PyPlot, DelimitedFiles

include("plot_ntws.jl")

run = 3695; ntw = "cyc1_18"; tit = "1-cycle, n=18"; labs = Dict{Int64,String}(0 => "q=0", 1 => "q=1", 2 => "q=-1")
#run = 6786; ntw = "cyc2_18"; tit = "2-cycle, n=18"; labs = Dict{Int64,String}(0 => "q=[0,0]", 1 => "q=[1,-1]", 2 => "q=[-1,1]")
#run = 6822; ntw = "cyc2_12"; tit = "2-cycle, n=12"; labs = Dict{Int64,String}(0 => "q=[0,0]", 1 => "q=[1,-1]", 2 => "q=[-1,1]")

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
PyPlot.fill(αs[[id;di]],[β2[1,id];β2[2,di]],color="C2",alpha=.5,label=labs[2])
PyPlot.fill(αs[[id;di]],[β0[1,id];β0[2,di]],color="C0",alpha=.5,label=labs[0])
PyPlot.fill(αs[[id;di]],[β1[1,id];β1[2,di]],color="C1",alpha=.5,label=labs[1])
title(tit*" (run: $run)")
xlabel("α")
ylabel("width")
legend()

subplot(1,2,2)
plot_ntws(ntw)

