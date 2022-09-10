using PyPlot,PyCall,DelimitedFiles,SparseArrays,Statistics
import PyPlot.plt
@pyimport matplotlib as mpl

include("NR.jl")

ntw = "uk"
ranking_type = "init"
ranking_measure = "Omega"
P0 = .6

function load_tracking(ntw::String,ranking_type::String,ranking_measure::String,P0::Float64)
	global mpl
	
	rmvd = Int.(vec(readdlm("data/rmvd_1b1_"*ntw*"_$(P0)_"*ranking_type*"_"*ranking_measure*".csv",',')))
	
	Bsp = readdlm(ntw*"_data/"*ntw*"_inc_mat.csv",',')
	B = sparse(Bsp[:,1],Bsp[:,2],Bsp[:,3])
	n,m = size(B)
	w = Bsp[2*(1:m),4]
	
	P = P0*vec(readdlm(ntw*"_data/P_"*ntw*".csv",','))
	P .-= mean(P)
	
	th,n_iter = NR_kuramoto(B,w,P,zeros(n))
	
	loads = Array{Float64,2}(undef,m,0)
	p00 = Array{Float64,1}()
	p10 = Array{Float64,1}()
	p25 = Array{Float64,1}()
	p50 = Array{Float64,1}()
	p75 = Array{Float64,1}()
	p90 = Array{Float64,1}()
	p100 = Array{Float64,1}()
	
	load = vec(abs.(sin.(transpose(B)*th)))
	loads = [loads load]
	push!(p00,minimum(load))
	push!(p10,quantile(load,.1))
	push!(p25,quantile(load,.25))
	push!(p50,quantile(load,.5))
	push!(p75,quantile(load,.75))
	push!(p90,quantile(load,.9))
	push!(p100,maximum(load))
	
	for i in 1:length(rmvd)
#		global m,rmvd,B,w,P,th,loads,p00,p25,p50,p75,p100
		
		nrmvd = setdiff(1:m,rmvd[1:i])
		Br = B[:,nrmvd]
		wr = w[nrmvd]
		wtest = ones(m)
		wtest[rmvd[1:i]] = zeros(i)
		
		th,n_iter = NR_kuramoto(Br,wr,P,th)
		
		load = wtest.*vec(abs.(sin.(transpose(B)*th)))
		loads = [loads load]
		
		push!(p00,minimum(load[nrmvd]))
		push!(p10,quantile(load[nrmvd],.1))
		push!(p25,quantile(load[nrmvd],.25))
		push!(p50,quantile(load[nrmvd],.5))
		push!(p75,quantile(load[nrmvd],.75))
		push!(p90,quantile(load[nrmvd],.9))
		push!(p100,maximum(load[nrmvd]))
	end
		
	figs = get_fignums()
	if length(figs) > 0
		num = maximum(get_fignums()) + 1
	else
		num = 1
	end
#=
	for i in 1:size(loads)[1]
		figure(num)
		PyPlot.plot(loads[i,:],"-o")
	end
=#
	decs = [0:length(p10)-1 p10;length(p90)-1:-1:0 p90[end:-1:1]]
	quarts = [0:length(p25)-1 p25;length(p75)-1:-1:0 p75[end:-1:1]]
	
	patches = mpl.pymember("patches")
	
	fig = plt.figure(num,figsize=(8,5))
	ax = fig[:add_axes]([.1,.1,.8,.8])
	
	p = patches[:Polygon](decs,closed=true,edgecolor="none",facecolor="blue",alpha=.2,rasterized=true)
	ax[:add_patch](p)
	
	p = patches[:Polygon](quarts,closed=true,edgecolor="none",facecolor="blue",alpha=.5,rasterized=true)
	ax[:add_patch](p)
	
	PyPlot.plot(p50,"-ok")
	PyPlot.plot(p00,"xk")
	PyPlot.plot(p100,"xk")
	
	PyPlot.plot([0,length(p25)-1],[1.,1.],":k")
	
	title(ntw*", "*ranking_measure*", "*ranking_type*", $(P0): Load distribution")
	xlabel("Removed lines")
	xticks(1:length(rmvd),labels=rmvd)
	ylabel("Relative load")
		
	return loads
end






