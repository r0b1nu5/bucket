using DelimitedFiles

include("load_rate.jl")

ntw = "ieee300"
ntw_data = PowerModels.parse_file("data/IEEE_300_bus.raw"; import_all=true)
ntw_data["dcline"] = Dict{String,Any}()

 #=
case = "case01"
rescale_gen = [("24",1.),]
rescale_load = [("32",1.),]
rescale_admittance = [("29",1.),]
# =#

 #=
case = "case02"
rescale_gen = []
rescale_load = [(k,1.033) for k in keys(ntw_data["load"])]
rescale_admittance = []
# =#

# #=
case = "case03"
rescale_gen = []
rescale_load = [(k,1.03738) for k in keys(ntw_data["load"])]
rescale_admittance = []
# =#

for rg in rescale_gen
	ntw_data["gen"][rg[1]]["pg"] *= rg[2]
end

for rl in rescale_load
	ntw_data["load"][rl[1]]["pd"] *= rl[2]
	ntw_data["load"][rl[1]]["qd"] *= rl[2]
end

for ra in rescale_admittance
	ntw_data["branch"][ra[1]]["br_r"] *= ra[2]
	ntw_data["branch"][ra[1]]["br_x"] *= ra[2]
end

compute_ac_pf!(ntw_data)

θ = angle_vec(ntw_data)
v = volt_vec(ntw_data)
ad_ratio = findnz(anglediff_rate(ntw_data))
fl_ratio = findnz(flow_rate(ntw_data))
active_fl = findnz(active_flow(ntw_data))

writedlm("cases/"*ntw*"_"*case*"_angles.csv",θ,',')
writedlm("cases/"*ntw*"_"*case*"_volt.csv",v,',')
writedlm("cases/"*ntw*"_"*case*"_angleratio.csv",[ad_ratio[1] ad_ratio[2] ad_ratio[3]],',')
writedlm("cases/"*ntw*"_"*case*"_flowratio.csv",[fl_ratio[1] fl_ratio[2] fl_ratio[3]],',')
writedlm("cases/"*ntw*"_"*case*"_flows.csv",[active_fl[1] active_fl[2] active_fl[3]],',')


