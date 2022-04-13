using PyPlot, DelimitedFiles

include("journals.jl")

ref_year = "00"

#journal = "science"
journal = "prl"
#journal = "lancet"

yi = 99
yf = 8

Nmin = 100
Nmax = 200

xi = readdlm("data/time_split/"*journal*"_$(ref_year)_$yi.txt",'\t')[2:end-2,1:2]
auti = xi[:,1]
numi = Int64.(xi[:,2])

xf = readdlm("data/time_split/"*journal*"_$(ref_year)_$yf.txt",'\t')[2:end-2,1:2]
autf = xf[:,1]
numf = Int64.(xf[:,2])

as = Vector{String}()
ni = Vector{Int64}()
nf = Vector{Int64}()
for k in 1:length(auti)
	if Nmin .<= numi[k] .<= Nmax
		idx = findmax(autf .== auti[k])[2]
		if numi[k] < numf[idx]
			push!(as,auti[k])
			push!(ni,numi[k])
			push!(nf,numf[idx])
		end
	end
end

ρi = Vector{Float64}()
ρf = Vector{Float64}()
for i in 1:length(as)-1
	@info "$i/$(length(as))"
	for j in i+1:length(as)
		push!(ρi,ni[i]/ni[j])
		push!(ρf,nf[i]/nf[j])
	end
end

PyPlot.plot(ρi,ρf,"o",color=journals_colors[journal][1])


