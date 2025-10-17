using DelimitedFiles, PowerModels, PyPlot

function load_testcase_1()
	bus = readdlm("rv/data-marc/bus.csv",',')
	load = readdlm("rv/data-marc/load.csv",',')
	line = readdlm("rv/data-marc/line.csv",',')
	
	network_data = Dict{String,Any}("name" => "tapis-testcase-1",
					"bus" => Dict{String,Any}(),
					"gen" => Dict{String,Any}(),
					"branch" => Dict{String,Any}(),
					"load" => Dict{String,Any}(),
					"dcline" => Dict{String,Any}(),
					"switch" => Dict{String,Any}(),
					"storage" => Dict{String,Any}(),
				        "shunt" => Dict{String,Any}(),
					"per_unit" => true
					)
	
	for i in 2:size(bus)[1]
		network_data["bus"]["$(bus[i,1])"] = Dict{String,Any}("number" => bus[i,1],
								      "bus_i" => bus[i,1],
								      "bus_type" => 1,
								      "vmin" => .9,
								      "vmax" => 1.1,
								      "index" => bus[i,1],
								      "va" => 0.,
								      "vm" => 1.,
								      "base_kv" => bus[i,3]
								      )
	end
	network_data["bus"]["1"]["bus_type"] = 3 # Slack bus
	network_data["gen"]["1"]["pg"] = 0.
	network_data["gen"]["1"]["qg"] = 0.
	network_data["gen"]["1"]["gen_bus"] = 1
	network_data["gen"]["1"]["pg"] = 0.
	network_data["gen"]["1"]["pmax"] = Inf
	network_data["gen"]["1"]["pmin"] = -Inf
	network_data["gen"]["1"]["qmax"] = Inf
	network_data["gen"]["1"]["qmin"] = -Inf
	network_data["gen"]["1"]["vg"] = 1.
	network_data["gen"]["1"]["gen_bus"] = 1
	network_data["gen"]["1"]["index"] = 1
	
	for i in 2:size(load)[1]
		network_data["load"]["$(load[i,1]+1)"] = Dict{String,Any}("load_bus" => load[i,3],
									  "status" => 1,
									  "pd" => load[i,4],
									  "qd" => load[i,5],
									  "index" => load[i,1]+1
									  )
	end
	
	for i in 2:size(line)[1]
		network_data["branch"]["$(line[i,1]+1)"] = Dict{String,Any}("name" => line[i,2],
									    "f_bus" => line[i,4],
									    "t_bus" => line[i,5],
									    "br_r" => line[i,6]*line[i,7],
									    "br_x" => line[i,6]*line[i,8],
									    "angmin" => -π/6,
									    "angmax" => π/6,
									    "transformer" => false,
									    "tap" => 1.,
									    "index" => line[i,1]+1,
									    "br_status" => 1,
									    "shift" => 0.,
	# See computation of Y_{0,1,2} at https://pandapower.readthedocs.io/en/develop/elements/line.html
									    "g_fr" => .5*line[i,10]*1e-6*line[i,6],
									    "g_to" => .5*line[i,10]*1e-6*line[i,6],
									    "b_fr" => .5*2*π*50*line[i,9]*1e-9*line[i,6],
									    "b_to" => .5*2*π*50*line[i,9]*1e-9*line[i,6]
									    )
	end

	return network_data
end


function set_pq_and_solve(p::Dict{String,Float64}, q::Dict{String,Float64}, nd::Dict{String,Any})
	n = length(p)

	for k in keys(p)
		nd["load"][k]["pd"] = p[k]*1e-6
		nd["load"][k]["qd"] = q[k]*1e-6
	end

	sol = PowerModels.compute_ac_pf(nd,show_trace=true)

	v = Float64[]
	for k in keys(p)
		push!(v,sol["solution"]["bus"]["$(nd["load"][k]["load_bus"])"]["vm"])
	end

	return v
end



