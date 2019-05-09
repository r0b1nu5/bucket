using DelimitedFiles,Statistics,PyPlot

function parse_ranking_data(ntw::String,ranking_measure::String,P0::Float64,n_rand::Int=0)
# ========================== INITIAL RANKING ==========================
	x = readdlm("data/ranks_1b1_"*ntw*"_$(P0)_init_"*ranking_measure*".csv",',')
	ranks1 = Array{Array{Int64,1},1}()
	for i in 1:size(x)[1]
		push!(ranks1,Array{Int64,1}(vec(x[i,1:size(x)[2]-i+1])))
	end
	
#	x = readdlm("data/cuts_1b1_"*ntw*"_$(P0)_init_"*ranking_measure*".csv",',')
	cuts1 = Array{Array{Int64,1},1}()
	#for i in 1:size(x)[1]
		#push!(cuts1,Array{Int64,1}(vec(x[i,1:size(x)[2]-i+1])))
	#end
	
	rmvd1 = Array{Int64,1}(vec(readdlm("data/rmvd_1b1_"*ntw*"_$(P0)_init_"*ranking_measure*".csv",',')))
	
	m1 = length(ranks1)

# ========================== UPDATED RANKING ==========================
	x = readdlm("data/ranks_1b1_"*ntw*"_$(P0)_updated_"*ranking_measure*".csv",',')
	ranks2 = Array{Array{Int64,1},1}()
	for i in 1:size(x)[1]
		push!(ranks2,Array{Int64,1}(vec(x[i,1:size(x)[2]-i+1])))
	end
	
	#x = readdlm("data/cuts_1b1_"*ntw*"_$(P0)_updated_"*ranking_measure*".csv",',')
	cuts2 = Array{Array{Int64,1},1}()
	#for i in 1:size(x)[1]
		#push!(cuts2,Array{Int64,1}(vec(x[i,1:size(x)[2]-i+1])))
	#end
	
	rmvd2 = Array{Int64,1}(vec(readdlm("data/rmvd_1b1_"*ntw*"_$(P0)_updated_"*ranking_measure*".csv",',')))
	
	m2 = length(ranks2)
	
# ========================== RANDOM RANKING ===========================
if n_rand > 0
	ranks3 = Array{Array{Array{Int64,1},1},1}()
	cuts3 = Array{Array{Array{Int64,1},1},1}()
	rmvds3 = Array{Array{Int64,1},1}()
	ms = Array{Float64,1}()
	for k in 1:n_rand
		x = readdlm("data/random/ranks_1b1_"*ntw*"_$(P0)_random_$k.csv",',')
		ranks = Array{Array{Int64,1},1}()
		for i in 1:size(x)[1]
			push!(ranks,Array{Int64,1}(vec(x[i,1:size(x)[2]-i+1])))
		end
		push!(ranks3,ranks)
		
		#x = readdlm("data/random/cuts_1b1"*ntw*"_$(P0)_random_$k.csv",',')
		cuts = Array{Array{Int64,1},1}()
		#for i in 1:size(x)[1]
		#	push!(cuts,Array{Int64,1}(vec(x[i,1:size(x)[2]-i+1])))
		#end
		push!(cuts3,cuts)
		
		push!(rmvds3,Array{Int64,1}(vec(readdlm("data/random/rmvd_1b1_"*ntw*"_$(P0)_random_$k.csv",','))))
		
		push!(ms,length(ranks))
	end
	
	m3 = mean(ms)
	min_3 = minimum(ms)
	p25_3 = quantile(ms,.25)
	p50_3 = quantile(ms,.5)
	p75_3 = quantile(ms,.75)
	max_3 = maximum(ms)
else	
	m3 = 0.
	min_3 = 0.
	p25_3 = 0.
	p50_3 = 0.
	p75_3 = 0.
	max_3 = 0.
end

# ========================== LAITINI RANKING ==========================
	x = readdlm("data/ranks_1b1_"*ntw*"_$(P0)_tini_"*ranking_measure*".csv",',')
	ranks4 = Array{Array{Int64,1},1}()
	for i in 1:size(x)[1]
		push!(ranks4,Array{Int64,1}(vec(x[i,1:size(x)[2]-i+1])))
	end
	
	#x = readdlm("data/cuts_1b1_"*ntw*"_$(P0)_tini_"*ranking_measure*".csv",',')
	cuts4 = Array{Array{Int64,1},1}()
	#for i in 1:size(x)[1]
	#	push!(cuts4,Array{Int64,1}(vec(x[i,1:size(x)[2]-i+1])))
	#end
	
	rmvd4 = Array{Int64,1}(vec(readdlm("data/rmvd_1b1_"*ntw*"_$(P0)_tini_"*ranking_measure*".csv",',')))
	
	m4 = length(ranks4)
	
# ========================== DETADPU RANKING ==========================
	x = readdlm("data/ranks_1b1_"*ntw*"_$(P0)_detadpu_"*ranking_measure*".csv",',')
	ranks5 = Array{Array{Int64,1},1}()
	for i in 1:size(x)[1]
		push!(ranks5,Array{Int64,1}(vec(x[i,1:size(x)[2]-i+1])))
	end
	
	#x = readdlm("data/cuts_1b1_"*ntw*"_$(P0)_detadpu_"*ranking_measure*".csv",',')
	cuts5 = Array{Array{Int64,1},1}()
	#for i in 1:size(x)[1]
	#	push!(cuts5,Array{Int64,1}(vec(x[i,1:size(x)[2]-i+1])))
	#end
	
	rmvd5 = Array{Int64,1}(vec(readdlm("data/rmvd_1b1_"*ntw*"_$(P0)_detadpu_"*ranking_measure*".csv",',')))
	
	m5 = length(ranks5)
	
	return m1,m2,m3,[min_3,p25_3,p50_3,p75_3,max_3],m4,m5
end

function plot_ranking_data(ntw::String,ranking_measure::String,Ps::Array{Float64,1},n_rand::Int64)
	m1 = Array{Float64,1}()
	m2 = Array{Float64,1}()
	m3 = Array{Float64,1}()
	p3 = Array{Float64,2}(undef,5,0)
	m4 = Array{Float64,1}()
	m5 = Array{Float64,1}()
	for P0 in Ps
		@info "P0 = $P0"
		
		X = parse_ranking_data(ntw,ranking_measure,P0,n_rand)
		push!(m1,X[1])
		push!(m2,X[2])
		push!(m3,X[3])
		p3 = [p3 X[4]]
		push!(m4,X[5])
		push!(m5,X[6])
	end
	
	n = 0
	m = 0
	if ntw == "uk"
		n = 120
		m = 165
	elseif ntw == "ieee57"
		n = 57
		m = 78
	elseif ntw == "ieee118"
		n = 118
		m = 179
	elseif ntw == "ieee300"
		n = 300
		m = 408
	else
		n = 1
		m = 0
	end
	
	rm_legend = ranking_measure
	if ranking_measure == "Omega"
		rm_legend = "b*Ω"
	elseif ranking_measure == "load"
		rm_legend = "Relative load [%]"
	elseif ranking_measure == "Omega+load"
		rm_legend = "Ω*(b*Δθ)^2"
	end
	
	figure()
	PyPlot.plot([minimum(Ps),maximum(Ps)],[m-n+1,m-n+1],":k",linewidth=1.5)
	PyPlot.plot(Ps,m1,"-o",label="Initial ranking")
	PyPlot.plot(Ps,m2,"--",label="Updated ranking")
	PyPlot.plot(Ps,p3[3,:],"-o",label="Random ranking (median)")
	PyPlot.plot(Ps,p3[2,:],"--k",label="Quartiles")
	PyPlot.plot(Ps,p3[4,:],"--k")
	PyPlot.plot(Ps,m4,"-o",label="Laitini ranking")
	PyPlot.plot(Ps,m5,"--",label="Detadpu ranking")
	xlabel("P0")
	ylabel("Number of lines")
	title(ntw*", "*rm_legend*": Number of lines to cut before no sync state")
	legend()
end




