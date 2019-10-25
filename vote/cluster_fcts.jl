using Distributed, DelimitedFiles, Distributions

include("scripts.jl")
include("big_rand.jl")

function eff(xes::Tuple{Array{Float64,1},Float64,String,Int64})
	x0 = xes[1]
	eps = xes[2]
	strat = xes[3]
	thr_id = xes[4]
	if strat == "fiedler"
		ef = influence_effort_fiedler(x0,eps)
	elseif strat == "mini"
		ef = influence_effort_mini(x0,eps)
	else
		ef = influence_effort_rand(x0,eps)
	end
	
	writedlm("data/eff_x$(thr_id)_"*strat*"_$eps.csv",ef,',')
end
				    







