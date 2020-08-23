using DelimitedFiles

include("journals.jl")

js = journals_short
js = ["nature_1950","energy_2005","ieee_trans_autom_control_2000","j_acs_1930","j_acs","lancet_1910","neng_j_med_1950","plant_cell_2000","pnas_1950","science_1940"]
for j in js
	num = Array{Float64,1}(vec(readdlm("./data/"*j*".txt",'\t')[2:end-2,2]))
	
	dat = sort(union(num))
	
	val_num = Array{Int64,2}(undef,2,0)
	
	for n in dat
		val_num = [val_num [n,Int(sum(num .== n))]]
	end
	
	writedlm("data/"*j*"_parsed.csv",val_num,',')
end









