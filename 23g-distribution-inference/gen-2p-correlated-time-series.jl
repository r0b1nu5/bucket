using DelimitedFiles, PowerModels, PyPlot

bus = readdlm("rv/data-marc/bus.csv",',')
load = readdlm("rv/data-marc/load.csv",',')
line = readdlm("rv/data-marc/line.csv",',')

network_data = Dict{String,Any}("name" => "tapis-testcase-1",
				"bus" => Dict{String,Any}(),
				"gen" => Dict{String,Any}(),
				"branch" => Dict{String,Any}(),
				"load" => Dict{String,Any}(),
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
								    "br_status" => 1
								    )
end








