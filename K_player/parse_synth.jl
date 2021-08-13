using DelimitedFiles

file = "synth01_2000_1000_70"
num = Int64.(vec(readdlm("./data/"*file*".csv",',')))

dat = sort(union(num))
	
val_num = Array{Int64,2}(undef,2,0)
	
for n in dat
	global val_num = [val_num [n,Int(sum(num .== n))]]
end
	
writedlm("data/"*file*"_parsed.csv",val_num,',')










