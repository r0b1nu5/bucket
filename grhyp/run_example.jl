include("load_rate.jl")

ntw_data = PowerModels.parse_file("data/IEEE_RTS96_bus.raw"; import_all=true)

rescale_gen = [("24",1.),]
for rg in rescale_gen
	ntw_data["gen"][rg[1]]["pg"] *= rg[2]
end

rescale_load = [("32",1.),]
for rl in rescale_load
	ntw_data["load"][rl[1]]["pd"] *= rl[2]
	ntw_data["load"][rl[1]]["qd"] *= rl[2]
end

rescale_admittance = [("29",1.),]
for ra in rescale_admittance
	ntw_data["branch"][ra[1]]["br_r"] *= ra[2]
	ntw_data["branch"][ra[1]]["br_x"] *= ra[2]
end

compute_ac_pf!(ntw_data)

load_rate(ntw_data)



