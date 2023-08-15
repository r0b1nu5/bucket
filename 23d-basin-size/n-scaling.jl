using Statistics, DelimitedFiles

include("tools.jl")

ns = [83,]
ns = [23,43,83,123,163]

type = "undir" # "dir", "undir"
n_intra = 1000

Mf_slope = Float64[]
Mi_slope = Float64[]
Σf = Float64[]
Σi = Float64[]
Λf = Float64[]
Λi = Float64[]

Qs = Dict{Any,Any}()

for n in ns
	qmax = floor(Int64,n/4)
	ids = Int64.(vec(readdlm("data/ids-"*type*"-$n-$(n_intra).csv",',')))
	global Qs = load_Qs("data/",ids,type,n,n_intra)

	xxx = get_normal_par(Qs)
	push!(Mf_slope,xxx[1])	
	push!(Mi_slope,xxx[2])	
	push!(Σf,xxx[3])	
	push!(Σi,xxx[4])
	yyy = get_exp_par(Qs,qmax)
	push!(Λf,yyy[3])
	push!(Λi,yyy[4])
end

if length(ns) > 1
	figure("n-scaling")
	PyPlot.plot(ns,Mf_slope,"-o",color="C0",label="μf slope")
	PyPlot.plot(ns,Σf./sqrt.(ns),"--o",color="C0",label="σf/√n")
	PyPlot.plot(ns,Λf.*sqrt.(ns),":o",color="C0",label="λf*√n")
	PyPlot.plot(ns,Mi_slope,"-o",color="C1",label="μi slope")
	PyPlot.plot(ns,Σi./sqrt.(ns),"--o",color="C1",label="σi/√n")
	PyPlot.plot(ns,Λi.*sqrt.(ns),":o",color="C1",label="λi*√n")
	legend()
else
	plot_Qs(Qs)
	@info "μf_slope = $(Mf_slope[1])"
	@info "μi_slope = $(Mi_slope[1])"
	@info "σf = $(Σf[1]) (if normal dist.)"
	@info "σi = $(Σi[1]) (if normal dist.)"
	@info "λf = $(Λf[1]) (if expo. dist.)"
	@info "λi = $(Λi[1]) (if expo. dist.)"
end



