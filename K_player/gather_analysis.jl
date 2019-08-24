using DelimitedFiles

#include("journals.jl")
journals_shors = ["prd",]

distri = ["pl","plc","yule","expo","poisson"]

n_threads = 25
n_samp = 100

for j in journals_short
	for dist in distri
		params = vec(readdlm("analysis/"*j*"_"*dist*"_params_$(n_samp)_1.csv",','))
		
		p_gofs = Array{Float64,1}()
		for i in 1:n_threads
			push!(p_gofs,readdlm("analysis/"*j*"_"*dist*"_p-gof_$(n_samp)_$i.csv",',')[1])

			rm("analysis/"*j*"_"*dist*"_params_$(n_samp)_$i.csv")
			rm("analysis/"*j*"_"*dist*"_p-gof_$(n_samp)_$i.csv")
		end
		
		p_gof = sum(p_gofs)/n_threads

		writedlm("analysis/"*j*"_"*dist*"_params_$(n_samp*n_threads).csv",params,',')
		writedlm("analysis/"*j*"_"*dist*"_p-gof_$(n_samp*n_threads).csv",p_gof,',')
	end
end



