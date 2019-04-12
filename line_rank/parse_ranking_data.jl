using DelimitedFiles,Statistics

function parse_rank_data(P0::Float64,n_rand::Int)
# =================== INITIAL RANKING ====================
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
	
	m1 = length(ranks1)
	
# =================== UPDATED RANKING ========================
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
	
	m2 = length(ranks2)

# ================= RANDOM RANKING ============================
	ranks3 = Array{Array{Array{Int64,1},1},1}()
	cuts3 = Array{Array{Array{Int64,1},1},1}()
	rmvds3 = Array{Array{Int64,1},1}()
	ms = Array{Float64,1}()
	for k in 1:n_rand
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
		
		push!(ms,length(ranks))
	end
	
	m3 = mean(ms)
	min_3 = minimum(ms)
	p25_3 = quantile(ms,.25)
	p50_3 = quantile(ms,.5)
	p75_3 = quantile(ms,.75)
	max_3 = maximum(ms)
	
	return m1,m2,m3,[min_3,p25_3,p50_3,p75_3,max_3],ranks1,ranks2,ranks3,cuts1,cuts2,cuts3,rmvd1,rmvd2,rmvds3
end

function plot_data(Ps,n_rand)
	m1 = Array{Float64,1}()
	m2 = Array{Float64,1}()
	m3 = Array{Float64,1}()
	p3 = Array{Float64,2}(undef,5,0)
	for P0 in Ps
		X = parse_rank_data(P0,n_rand)
		push!(m1,X[1])
		push!(m2,X[2])
		push!(m3,X[3])
		p3 = [p3 X[4]]
	end
	
	figure()
	PyPlot.plot(Ps,m1,"-o",label="Inital ranking")
	PyPlot.plot(Ps,m2,"--o",label="Updated ranking")
	PyPlot.plot(Ps,p3[3,:],"-o",label="Random ranking (median)")
	PyPlot.plot(Ps,p3[2,:],"--k",label="Quartiles")
	PyPlot.plot(Ps,p3[4,:],"--k")
	xlabel("P0")
	ylabel("Number of lines")
	title("Number of lines to cut before no sync state")
	legend()
end





