using PyPlot, DelimitedFiles

include("plot_ntws.jl")
include("runs_dict.jl")

name_cmap = "viridis"
cmap = get_cmap(name_cmap)
# #=
col1 = cmap(1/6)
col2 = cmap(2.5/6)
col3 = cmap(4/6)
# =#
 #=
col1 = "C0"
col2 = "C1"
col3 = "C2"
# =#
 #=
base_col = (0,54,96)./255
base_col = (66,0,192)./255
ts = LinRange(0.,.5,3)
col1 = ts[1] .+ (1-ts[1]).*base_col
col2 = ts[2] .+ (1-ts[2]).*base_col
col3 = ts[3] .+ (1-ts[3]).*base_col
# =#

run = 6179 #run = 3695

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
PyPlot.fill([αs[id];αs[di]],[β2[1,id];zeros(T)]/sqrt(2),"k",alpha=.1)
PyPlot.plot(αs[id],β2[1,id]/sqrt(2),"-",color=col1,linewidth=lw,label=runs[run]["labs"][2])
#PyPlot.plot(αs[id],-β1[1,id],"--",color=col1,linewidth=lw)
PyPlot.fill([αs[id];αs[di]],[β0[1,id];zeros(T)]/sqrt(2),"k",alpha=.1)
PyPlot.plot(αs[id],β0[1,id]/sqrt(2),"-",color=col2,linewidth=lw,label=runs[run]["labs"][0])
#PyPlot.plot(αs[id],-β0[1,id],"--",color=col2,linewidth=lw)
PyPlot.fill([αs[id];αs[di]],[β1[1,id];zeros(T)]/sqrt(2),"k",alpha=.1)
PyPlot.plot(αs[id],β1[1,id]/sqrt(2),"-",color=col3,linewidth=lw,label=runs[run]["labs"][1])
#PyPlot.plot(αs[id],-β2[1,id],"--",color=col3,linewidth=lw)
title(runs[run]["tit"]*" (run: $run)")
xlabel("ϕ")
ylabel("β")
legend()

subplot(1,2,2)
w = vec(readdlm("temp_data/w_$(run).csv",','))
#θ = vec(readdlm("temp_data/th0f_$(run).csv",',')) .+ .4
θ = vec(readdlm("temp_data/th2f_$(run).csv",',')) .+ .5
#plot_ntws(runs[run]["ntw"],w,10.,2.,name_cmap)
plot_ntws(runs[run]["ntw"],w,θ,10.,2.,name_cmap)
axis([-1.3,1.3,-1.3,1.3])

figure("test",(8,8))
PyPlot.plot(cos.(LinRange(0,2π,200)),sin.(LinRange(0,2π,200)),"k")
PyPlot.plot(cos.([θ;θ[1]]),sin.([θ;θ[1]]),"ok")
PyPlot.plot(1.1*cos.(θ[1:7]),1.1*sin.(θ[1:7]))
PyPlot.plot(.9*cos.([θ[7:18];θ[1]]),.9*sin.([θ[7:18];θ[1]]))


