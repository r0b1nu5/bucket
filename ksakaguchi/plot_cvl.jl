using PyPlot, DelimitedFiles

ntws = ["cyc5","ntw20"]
ntw = ntws[1]

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
	#=
	for i in 1:1000
		PyPlot.plot(cohes[i],μs[i],".",color=cmap(α[i]))
	end
	=#

PyPlot.text(γbar,0,"γbar")
ylabel("μ_Π(Dh)")
xlabel("cohesiveness")



