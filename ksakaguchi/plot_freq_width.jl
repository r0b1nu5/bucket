using PyPlot, DelimitedFiles

include("plot_ntws.jl")
cmap = get_cmap("plasma")

run = 3695

runs = Dict{Int64,Any}(3695 => Dict{String,Any}("ntw" => "cyc1_18", 
						"tit" => "1-cycle, n=18", 
						"labs" => Dict{Int64,String}(0 => "q=0", 
									     1 => "q=1", 
									     2 => "q=-1"
									     ),
						),
		      6786 => Dict{String,Any}("ntw" => "cyc2_18", 
					       "tit" => "2-cycle, n=18", 
					       "labs" => Dict{Int64,String}(0 => "q=[0,0]", 
									    1 => "q=[1,-1]", 
									    2 => "q=[-1,1]"
									    ),
					       ),
		      6822 => Dict{String,Any}("ntw" => "cyc2_12", 
					       "tit" => "2-cycle, n=12", 
					       "labs" = Dict{Int64,String}(0 => "q=[0,0]", 
									   1 => "q=[1,-1]", 
									   2 => "q=[-1,1]"
									   ),
					       ),
		      )

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
PyPlot.plot(αs[id],β2[2,id],"-",color=cmpa(0.))
PyPlot.fill(αs[[id;di]],[β0[1,id];β0[3,di]],color=cmap(1/3),alpha=.3,label=runs[run]["labs"][0])
PyPlot.plot(αs[id],β0[1,id],"-o",color=cmap(1/3))
PyPlot.plot(αs[id],β0[2,id],"-",color=cmpa(1/3))
PyPlot.fill(αs[[id;di]],[β1[1,id];β1[3,di]],color=cmap(2/3),alpha=.3,label=runs[run]["labs"][1])
PyPlot.plot(αs[id],β1[1,id],"-o",color=cmap(2/3))
PyPlot.plot(αs[id],β1[2,id],"-",color=cmpa(2/3))
title(runs[run]["tit"]*" (run: $run)")
xlabel("α")
ylabel("width")
legend()

subplot(1,2,2)
plot_ntws(runs[run]["ntw"])

