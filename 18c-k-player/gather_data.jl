using DelimitedFiles

include("journals.jl")

journals_short = ["prd","chaos-18"]

for j in journals_short
	num = Array{Int64,1}(vec(readdlm("./data/"*j*".txt",'\t')[2:end-2,2]))
	
	vals = intersect(num)
	
	distr = Array{Float64,1}()
	for v in vals
		push!(distr,sum(num .== v)/length(num))
	end
	
	writedlm("./data/"*j*"_vals.txt",vals,',')
	writedlm("./data/"*j*"_distr.txt",distr,',')
end
