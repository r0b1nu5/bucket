using DelimitedFiles

include("journals.jl")

js = journals_short

for j in js
	num = Array{Float64,1}(vec(readdlm("./data/"*j*".txt",'\t')[2:end-2,2]))
	
	dat = sort(union(num))
	
	val_num = Array{Int64,2}(undef,2,0)
	
	for n in dat
		val_num = [val_num [n,Int(sum(num .== n))]]
	end
	
	writedlm("data/"*j*"_parsed.csv",val_num,',')
end









