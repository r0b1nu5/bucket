using PyPlot, DelimitedFiles, LinearAlgebra

ntws = ["cyc5","ntw20"]
ntw = ntws[1]

cols = [(161,207,239)./255,(76,163,224)./255,(31,119,180)./255]

c = 500
	
	cohes = Vector{Float64}()
	avcoh = Vector{Float64}()
	μs = Vector{Float64}()
	
	for i in 1:c
		global cohes,avcoh,μs
	
		x = readdlm("temp_data/"*ntw*"_$i.csv",',')
		cohes = [cohes;x[:,1]]
		avcoh = [avcoh;x[:,2]]
		μs = [μs;x[:,3]]
	end
	
	L = readdlm("data/"*ntw*"_L.csv",',')
	n = size(L)[1]
	λ2 = eigvals(L)[2]
	δ = maximum(diag(L))
	ϕ = .3
	γbar = atan(λ2/(δ*tan(ϕ)))
	
	cm = Dict{String,String}("cyc5" => "Blues", "ntw20" => "Greens")
	cmap = get_cmap(cm[ntw])
	α = .7*(avcoh .- minimum(avcoh))./(maximum(avcoh) - minimum(avcoh)) .+ .3
	
	PyPlot.plot([0,maximum(cohes)],[0,0],"k")
	PyPlot.plot([0,0],[min(0.,minimum(μs)),max(0.,maximum(μs))],"k")
	PyPlot.fill([0,γbar,γbar,0,0],[0,0,maximum(μs),maximum(μs),0],color=(.8,.8,.8,1.))
	PyPlot.plot(cohes[1:1000],μs[1:1000],".")
	#PyPlot.plot(cohes,μs,".",alpha=.005)
	#=
	for i in 1:1000
		PyPlot.plot(cohes[i],μs[i],".",color=cmap(α[i]))
	end
	=#

PyPlot.text(γbar,0,"γbar")
ylabel("μ_Π(Dh)")
xlabel("cohesiveness")

res = 100
cs0 = LinRange(minimum(cohes),maximum(cohes),res)
dc = cs0[2]-cs0[1]
cs = [cs0[1] - .5*dc;cs0 .+ .5*dc]
μmax = Float64[]
μmin = Float64[]
for i in 1:res
	ids = setdiff((1:length(cohes)).*(cs[i] .<= cohes .< cs[i+1]),[0,])
	push!(μmax,maximum(μs[ids]))
	push!(μmin,minimum(μs[ids]))
end

figure()

PyPlot.plot([0,maximum(cohes)],[0,0],"k")
PyPlot.plot([0,0],[min(0.,minimum(μs)),max(0.,maximum(μs))],"k")
PyPlot.fill([0,γbar,γbar,0,0],[0,0,maximum(μs),maximum(μs),0],color=(.8,.8,.8,1.))
PyPlot.fill([cs0[1:res];cs0[res:-1:1]],[μmax;μmin[res:-1:1]],color=cols[1])
PyPlot.plot(cohes[1:2000],μs[1:2000],".",markersize=2,color=cols[3])
PyPlot.text(γbar,0,"γbar")
ylabel("μ_Π(Dh)")
xlabel("cohesiveness")

