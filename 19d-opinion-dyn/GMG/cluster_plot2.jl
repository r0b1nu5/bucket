using PyPlot, DelimitedFiles, Statistics

include("scripts.jl")

# Generates plots, based on the data obtained with "cluster_run2.jl". 
# (1) effort vs. inital outcome
# (2) effort vs. sum of comps in final state

nx = 100
n0 = 299
emi = .15
ema = .35
ne = 3
epss = Array(LinRange(emi,ema,ne))
n_run = 20
d = 1.
sig = .2
n1 = 1000
n2 = 1001
n = n1 + n2
n_modes = [2,3,4]
symbols = ["o","x","s","^","*","D","+"]

for j in 1:ne
	os = Array{Float64,1}()
	sxs = Array{Float64,1}()
	efr = Array{Float64,1}()
	eff = Array{Float64,2}(undef,0,length(n_modes))
	efm = Array{Float64,1}()
	
	for i in n0+1:n0+nx
		push!(os,readdlm("data/o_x$(i)_$(epss[j]).csv",',')[1])
		push!(sxs,sum(readdlm("data/x_x$(i)_$(epss[j]).csv",',')))
		xxx = Array{Float64,1}()
		for k in 1:n_run
			push!(xxx,readdlm("data/eff_x$(i)_rand$(k)_$(epss[j]).csv",',')[1])
		end
		push!(efr,mean(xxx))
		eff = [eff;vec(readdlm("data/eff_x$(i)_fiedler_$(epss[j]).csv",','))']
		push!(efm,readdlm("data/eff_x$(i)_mini_$(epss[j]).csv",',')[1])

	end

	figure(222)

	subplot(2,3,1)
	PyPlot.plot(os,efr,"o",color="C$(j-1)",label="r = $(round(cor(os,efr),digits=3)), ε = $(epss[j])")

	subplot(2,3,2)
	for m in 1:length(n_modes)
		PyPlot.plot(os,eff[:,m],symbols[m],color="C$(j-1)",label="mode $(n_modes[m]): r = $(round(cor(os,eff[:,m]),digits=3)), ε = $(epss[j])")
	end

	subplot(2,3,3)
	PyPlot.plot(os,efm,"o",color="C$(j-1)",label="r = $(round(cor(os,efm),digits=3)), ε = $(epss[j])")

	subplot(2,3,4)
	PyPlot.plot(sxs,efr,"o",color="C$(j-1)",label="r = $(round(cor(sxs,efr),digits=3)), ε = $(epss[j])")

	subplot(2,3,5)
	for m in 1:length(n_modes)
		PyPlot.plot(sxs,eff[:,m],symbols[m],color="C$(j-1)",label="mode $(n_modes[m]): r = $(round(cor(sxs,eff[:,m]),digits=3)), ε = $(epss[j])")
	end

	subplot(2,3,6)
	PyPlot.plot(sxs,efm,"o",color="C$(j-1)",label="r = $(round(cor(sxs,efm),digits=3)), ε = $(epss[j])")
end

figure(222)

subplot(2,3,1)
xlabel("initial outcome")
ylabel("ξ")
legend()
title("Random")

subplot(2,3,2)
xlabel("initial outcome")
ylabel("ξ")
legend()
title("Fiedler")

subplot(2,3,3)
xlabel("initial outcome")
ylabel("ξ")
legend()
title("Minimum")

subplot(2,3,4)
xlabel("Σx")
ylabel("ξ")
legend()

subplot(2,3,5)
xlabel("Σx")
ylabel("ξ")
legend()

subplot(2,3,6)
xlabel("Σx")
ylabel("ξ")
legend()






