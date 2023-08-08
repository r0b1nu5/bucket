using Statistics, DelimitedFiles

include("tools.jl")

ns = [83,]

type = "undir" # "dir", "undir"
n_intra = 1000

Mf_slope = Float64[]
Mi_slope = Float64[]
Σf = Float64[]
Σi = Float64[]

Qs = Dict{Any,Any}()

for n in ns
	ids = Int64.(vec(readdlm("data/ids-"*type*"-$n-$(n_intra).csv",',')))
	global Qs = load_Qs("data/",ids,type,n,n_intra)

	xxx = get_normal_par(Qs)
	push!(Mf_slope,xxx[1])	
	push!(Mi_slope,xxx[2])	
	push!(Σf,xxx[3])	
	push!(Σi,xxx[4])
end

if length(ns) > 1
	figure("n-scaling")
	PyPlot.plot(ns,Mf_slope,color="C0")
	PyPlot.plot(ns,Σf,"--",color="C0")
	PyPlot.plot(ns,Mi_slope,color="C1")
	PyPlot.plot(ns,Σi,"--",color="C1")
else
	plot_Qs(Qs)
	@info "μf_slope = $(Mf_slope[1])"
	@info "μi_slope = $(Mi_slope[1])"
	@info "σf = $(Σf[1])"
	@info "σi = $(Σi[1])"
end



