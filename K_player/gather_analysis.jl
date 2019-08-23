using DelimitedFiles

#include("journals.jl")
journals_shors = ["prd",]

for j in journals_short
	params = vec(readdlm("analysis/"*j*"_params_100_1.csv",','))
	
	p_gofs = zeros(3,25)
	for i in 1:25
		p_gofs[:,i] = vec(readdlm("analysis/"*j*"_p-gof_100_$i.csv",','))
	end
	
	ps = sum(p_gofs,dims=2)./25
	
#	writedlm("analysis/"*j*"_pl_params_2500.csv",params[1],',')
#	writedlm("analysis/"*j*"_pl_p-gof_2500.csv",ps[1],',')
	writedlm("analysis/"*j*"_plco_params_2500.csv",[params[2],params[3]],',')
	writedlm("analysis/"*j*"_plco_p-gof_2500.csv",ps[2],',')
#	writedlm("analysis/"*j*"_yule_params_2500.csv",params[4],',')
#	writedlm("analysis/"*j*"_yule_p-gof_2500.csv",ps[3],',')
end




