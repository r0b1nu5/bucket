using Distributed, DelimitedFiles, Distributions

include("scripts.jl")
include("big_rand.jl")

# Function computing the effort (eff), initial outcome (o), and initial state (x).

function eff(xes::Tuple{Array{Float64,1},Float64,String,Int64})
	x0 = xes[1]
	eps = xes[2]
	strat = xes[3]
	thr_id = xes[4]

	w0 = .1

	if strat == "fiedler"
		ef = Array{Float64,1}()
		for m in [2,3,4]
			eeff,o,x = influence_effort_fiedler(x0,eps,w0,m)
			push!(ef,eeff)
		end
	elseif strat == "mini"
		ef,o,x = influence_effort_mini(x0,eps,w0)
	else
		ef,o,x = influence_effort_rand(x0,eps,w0)
	end
	
	writedlm("data/eff_x$(thr_id)_"*strat*"_$eps.csv",ef,',')
	writedlm("data/o_x$(thr_id)_$eps.csv",o,',')
	writedlm("data/x_x$(thr_id)_$eps.csv",x,',')
end
				    







