using DelimitedFiles

include("journals.jl")

for j in journals_short
	params = readdlm("analysis/"*j*"_params_100_1.csv",',')
	
	p_gofs = zeros(3,25)
	for i in 1:25
		p_gofs[:,i] = vec(readdlm("analysis/"*j*"_p-gof_100_$i.csv",','))
	end
	
	ps = sum(p_gofs,dims=2)./25
	
	writedlm("analysis/"*j*"_params.csv",',')
	writedlm("analysis/"*j*"_p-gof.csv",',')
end





