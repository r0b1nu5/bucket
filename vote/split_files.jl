using DelimitedFiles

for i in 1:4
	for eps in LinRange(0.,1.,40)
		eo = vec(readdlm("data/eff_x$(i)_fiedler_$eps.csv",','))
		writedlm("data/eff_x$(i)_fiedler_$eps.csv",eo[1],',')
		writedlm("data/o_x$(i)_$eps.csv",eo[2],',')

		eo = vec(readdlm("data/eff_x$(i)_mini_$eps.csv",','))
		writedlm("data/eff_x$(i)_mini_$eps.csv",eo[1],',')
		writedlm("data/o_x$(i)_$eps.csv",eo[2],',')

		for j in 1:100
			eo = vec(readdlm("data/eff_x$(i)_rand$(j)_$eps.csv",','))
			writedlm("data/eff_x$(i)_rand$(j)_$eps.csv",eo[1],',')
			writedlm("data/o_x$(i)_$eps.csv",eo[2],',')
		end
	end
end







