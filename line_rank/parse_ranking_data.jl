using DelimitedFiles

P0 = 1.7

x = readdlm("./data/ranks_init_rank_$(P0)_rev.csv",',')
ranks1 = Array{Array{Int64,1},1}()
for i in 1:size(x)[1]
	push!(ranks1,Array{Int64,1}(vec(x[i,1:size(x)[2]-i+1])))
end

x = readdlm("./data/cuts_init_rank_$(P0)_rev.csv",',')
cuts1 = Array{Array{Int64,1},1}()
for i in 1:size(x)[1]
	push!(cuts1,Array{Int64,1}(vec(x[i,1:size(x)[2]-i+1])))
end

rmvd1 = Array{Int64,1}(vec(readdlm("./data/rmvd_init_rank_$(P0)_rev.csv",',')))

x = readdlm("./data/ranks_updated_rank_$(P0)_rev.csv",',')
ranks2 = Array{Array{Int64,1},1}()
for i in 1:size(x)[1]
	push!(ranks2,Array{Int64,1}(vec(x[i,1:size(x)[2]-i+1])))
end

x = readdlm("./data/cuts_updated_rank_$(P0)_rev.csv",',')
cuts2 = Array{Array{Int64,1},1}()
for i in 1:size(x)[1]
	push!(cuts2,Array{Int64,1}(vec(x[i,1:size(x)[2]-i+1])))
end

rmvd2 = Array{Int64,1}(vec(readdlm("./data/rmvd_updated_rank_$(P0)_rev.csv",',')))

ranks3 = Array{Array{Array{Int64,1},1},1}()
cuts3 = Array{Array{Array{Int64,1},1},1}()
rmvds3 = Array{Array{Int64,1},1}()
for k in 1:20
	x = readdlm("./data/ranks_random_$(P0)_$k.csv",',')
	ranks = Array{Array{Int64,1},1}()
	for i in 1:size(x)[1]
		push!(ranks,Array{Int64,1}(vec(x[i,1:size(x)[2]-i+1])))
	end
	push!(ranks3,ranks)
	
	x = readdlm("./data/cuts_random_$(P0)_$k.csv",',')
	cuts = Array{Array{Int64,1},1}()
	for i in 1:size(x)[1]
		push!(cuts,Array{Int64,1}(vec(x[i,1:size(x)[2]-i+1])))
	end
	push!(cuts3,cuts)
	
	push!(rmvds3,Array{Int64,1}(vec(readdlm("./data/rmvd_random_$(P0)_$k.csv",','))))
end


	




