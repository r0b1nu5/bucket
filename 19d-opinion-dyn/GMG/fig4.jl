using PyPlot, NPZ

syst = ["SR","WTAR","PR","HOR"]
dist = [1,2,3]
years = [2012,2014,2016,2018,2020]
n_iter = 100

N1 = 43960
N2 = 501*435
Ns = Dict{String,Int64}("SR" => N1, "WTAR" => N1, "PR" => N1, "HOR" => N2)

ϵs = Vector(LinRange(.05,.8,25))

ξ = Dict{Tuple{String,Int64},Vector{Float64}}((s,d) => zeros(25) for s in syst for d in dist)

for d in dist
	for y in years
		x = npzread("./data_gmg/US_STATE/RealDout_$(d)_$(y).npz")
		y = npzread("./data_gmg/US_DISTRICT/usrepub_$(d)_$(y).npz")
		
		ξ[("PR",d)] += vec(sum(x["E_PR"],dims=1))./(n_iter*length(years))
		ξ[("WTAR",d)] += vec(sum(x["E_WTAR"],dims=1))./(n_iter*length(years))
		ξ[("SR",d)] += vec(sum(x["E_SR"],dims=1))./(n_iter*length(years))
		ξ[("HOR",d)] += vec(sum(y["Effort_MM"],dims=1))./(n_iter*length(years))	
	end
end

shp = ["-","--",":"]
col = Dict{String,String}("SR" => "C0", "WTAR" => "C1", "PR" => "C2", "HOR" => "C3")

p = [.7,.5,1.]
μ = [0.,.3,.3]
Δ = [1.,1.,0.]
σ = [.3,.3,.3]

for d in dist
	figure(666)
	subplot2grid((3,3),(0,d-1),rowspan=2,colspan=1)
	mi = minimum([minimum(ξ[(s,d)]/Ns[s]) for s in syst])
	ma = maximum([maximum(ξ[(s,d)]/Ns[s]) for s in syst])
#	mi = minimum([minimum([minimum(ξ[(s,dd)]/Ns[s]) for s in syst]) for dd in dist])
#	ma = maximum([maximum([maximum(ξ[(s,dd)]/Ns[s]) for s in syst]) for dd in dist])
	
	for s in syst
		PyPlot.semilogy(ϵs/.2,ξ[(s,d)]/Ns[s],label=s)

#		figure(667)
#		PyPlot.semilogy(ϵs/.2,ξ[(s,d)]/Ns[s],shp[d],color=col[s],label=s)
	end
	
	xlabel("ϵ/σ")
	ylabel("ξ/N")
#	legend()
	axis([ϵs[1]/.2,ϵs[end]/.2,.95*mi,10*ma])

	
	x = LinRange(-2,2,200)
	y = p[d]*exp.(-(x .- μ[d] .- Δ[d]).^2/σ[d]^2) + (1-p[d])*exp.(-(x .- μ[d] .+ Δ[d]).^2/σ[d]^2)
	
	figure(666)
	subplot2grid((3,3),(2,d-1))
	PyPlot.plot([-2,2],[0,0],"k")
	PyPlot.plot([0,0],[-.05,max(p[d],1-p[d])+.05],"k")
	PyPlot.plot(x,y,linewidth=3)
	axis([-2,2,-.05,max(p[d],1-p[d])+.05])

end

#figure(667)
#xlabel("ϵ/σ")
#ylabel("ξ/N")

				

