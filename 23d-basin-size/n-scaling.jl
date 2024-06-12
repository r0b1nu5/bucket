using Statistics, DelimitedFiles

include("tools.jl")
include("volumes.jl")

#ns = [83,]
#ns = [23,43,83,123,163]
ns = [43,83,123,163]

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
	PyPlot.plot(ns,Λf./sqrt.(ns),":o",color="C0",label="λf*√n")
	PyPlot.plot(ns,Mi_slope,"-o",color="C1",label="μi slope")
	PyPlot.plot(ns,Σi./sqrt.(ns),"--o",color="C1",label="σi/√n")
	PyPlot.plot(ns,Λi./sqrt.(ns),":o",color="C1",label="λi*√n")
	legend()
	xlabel("n")
	ylabel("parameter")

	mμf = mean(Mf_slope)
	dμf = norm(Mf_slope .- mμf,Inf)
	eμf = round(dμf/mμf*100,sigdigits=4)
	mμi = mean(Mi_slope)
	dμi = norm(Mi_slope .- mμi,Inf)
	eμi = round(dμi/mμi*100,sigdigits=4)
	mσf = mean(Σf./sqrt.(ns))
	dσf = norm(Σf./sqrt.(ns) .- mσf,Inf)
	eσf = round(dσf/mσf*100,sigdigits=4)
	mσi = mean(Σi./sqrt.(ns))
	dσi = norm(Σi./sqrt.(ns) .- mσi,Inf)
	eσi = round(dσi/mσi*100,sigdigits=4)
	mλf = mean(Λf./sqrt.(ns))
	dλf = norm(Λf./sqrt.(ns) .- mλf,Inf)
	eλf = round(dλf/mλf*100,sigdigits=4)
	mλi = mean(Λi./sqrt.(ns))
	dλi = norm(Λi./sqrt.(ns) .- mλi,Inf)
	eλi = round(dλi/mλi*100,sigdigits=4)

	@info "Av. μf_slope = $(round(mμf,digits=2)) [±$(eμf)%]"
	@info "Av. μi_slope = $(round(mμi,digits=2)) [±$(eμi)%]"
	@info "Av. σf/√n = $(round(mσf,digits=2)) [±$(eσf)%] (if normal dist.)"
	@info "Av. σi/√n = $(round(mσi,digits=2)) [±$(eσi)%] (if normal dist.)"
	@info "Av. λf/√n = $(round(mλf,digits=2)) [±$(eλf)%] (if expo dist.)"
	@info "Av. λi/√n = $(round(mλi,digits=2)) [±$(eλi)%] (if expo dist.)"
else
	plot_Qs(Qs)
	@info "μf_slope = $(Mf_slope[1])"
	@info "μi_slope = $(Mi_slope[1])"
	@info "σf = $(Σf[1]) (if normal dist.)"
	@info "σi = $(Σi[1]) (if normal dist.)"
	@info "λf = $(Λf[1]) (if expo. dist.)"
	@info "λi = $(Λi[1]) (if expo. dist.)"
end



