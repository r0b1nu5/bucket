using DelimitedFiles

id = 100012
num = Int64.(vec(readdlm("./synth_data/d_$(id).csv",',')))

dat = sort(union(num))
	
val_num = Array{Int64,2}(undef,2,0)
	
for n in dat
	global val_num = [val_num [n,Int(sum(num .== n))]]
end
	
writedlm("data/d_$(id)_parsed.csv",val_num,',')










