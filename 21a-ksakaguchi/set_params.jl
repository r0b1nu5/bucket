using PowerModels

function set_angles(nd0::Dict{String,Any}, i::String, θ::Float64)
	nd = nd0

	nd["bus"][i]["va"] = θ

	return nd
end

function set_angles(nd0::Dict{String,Any}, ids::Vector{String}, θ::Vector{Float64})
	nd = nd0

	for i in 1:length(ids)
		nd["bus"][ids[i]]["va"] = θ[i]
	end

	return nd
end

function scale_apower(nd0::Dict{String,Any}, ρ::Float64)
	nd = nd0

	for k in keys(nd["gen"])
		nd["gen"][k]["pg"] *= ρ
	end

	for k in keys(nd["load"])
		nd["load"][k]["pd"] *= ρ
	end

	return nd
end

function scale_rpower(nd0::Dict{String,Any}, ρ::Float64)
	nd = nd0

	for k in keys(nd["load"])
		nd["load"][k]["qd"] *= ρ
	end
	
	return nd
end



