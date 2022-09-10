using DelimitedFiles,Statistics,LinearAlgebra,SparseArrays,Dates

include("rm_line.jl")
include("isconnected.jl")
include("kuramoto.jl")
include("res_dist.jl")
include("sync.jl")
include("NR.jl")

# Removes the lines of a network one by one until there is no fixed point anymore. Lines are removed such that the network is not splitted.

## INPUT
# ntw: network considered ("uk", "ieee57", "ieee118", "ieee300", or "pegase" (long to converge).)
# ranking_type: ranking used to remove the lines ("init": at each step, removes the most central line according to the measure in the initial network. "updaed": same but updates the ranking after removing each line. "random": removes lines randomly. "tini": same as "init" but removes the least central line. "detadpu": same as "updated" but removes the least central line.)
# ranking_measure: measure used to make the ranking ("Omega": resistance distance, "Womega": resistance distance of the Laplacian weighted with the cosines of the angles differences, "dKf1": based on the variation of Kf1, "dKf2": based on the variation of Kf2, "load": uses the load on the line, "Omega+load": uses the resistance distance rescaled by the load on the line.)

function rmv_1b1(ntw::String,ranking_type::String,ranking_measure::String,P0::Float64,M::Array{Float64,1},D::Array{Float64,1},max_iter::Int64=50,eps::Float64=1e-6,h::Float64=.1)
	
	Asp = readdlm(ntw*"_data/"*ntw*"_adj_mat.csv",',')
	Lsp = readdlm(ntw*"_data/"*ntw*"_lap_mat.csv",',')
	Bsp = readdlm(ntw*"_data/"*ntw*"_inc_mat.csv",',')
	
	A = sparse(Array{Int,1}(Asp[:,1]),Array{Int,1}(Asp[:,2]),vec(Asp[:,3]))
	L = sparse(Array{Int,1}(Lsp[:,1]),Array{Int,1}(Lsp[:,2]),vec(Lsp[:,3]))
	B = sparse(Array{Int,1}(Bsp[:,1]),Array{Int,1}(Bsp[:,2]),vec(Bsp[:,3]))
	w = Bsp[2*(1:Int(size(Bsp)[1]/2)),4]

	n,m = size(B)
		
	P = P0*vec(readdlm(ntw*"_data/P_"*ntw*".csv",','))
	P .-= mean(P)
	
	list_sync = readdir("sync_states")
	if length(intersect(list_sync,[ntw*"_sync_$P0.csv",])) == 0
		x1 = sync(ntw,P0,M,D)
	else
		x1 = vec(readdlm("sync_states/"*ntw*"_sync_$P0.csv",','))
	end
	
	if x1[1] == "nope"
		global cuts,ranks,rmvd,run
		@info "$(now()) -- "*ntw*": No sync possible!"
		cuts = Array{Array{Int64,1},1}()
		ranks = Array{Array{Int64,1},1}()
		run = false
	else
		ll = Array{Tuple{Int64,Int64},1}()
		mes = Array{Float64,1}()
		
		if ranking_measure == "Omega"
			Om = res_dist(L)
			for i in 1:m
				l = (Int(Bsp[2*i-1,1]),Int(Bsp[2*i,1]))
				push!(ll,l)
				push!(mes,-L[l[1],l[2]]*Om[l[1],l[2]])
			end
		elseif ranking_measure == "Womega"
			cdx = cos.(x1[1:n]*ones(1,n)-ones(n)*transpose(x1[1:n])).*(1 .- diagm(0 => ones(n)))
			Lw = L.*cdx
			Lw = diagm(0 => vec(sum(Lw,dims=2)))
			Om = res_dist(Lw)
			for i in 1:m
				l = (Int(Bsp[2*i-1,1]),Int(Bsp[2*i,1]))
				push!(ll,l)
				push!(mes,-Lw[l[1],l[2]]*Om[l[1],l[2]])
			end
		elseif ranking_measure == "dKf1"
			Om = res_dist(L)
			Om2 = res_dist(L,2)
			for i in 1:m
				l = (Int(Bsp[2*i-1,1]),Int(Bsp[2*i,1]))
				push!(ll,l)
				b = -L[l[1],l[2]]
				push!(mes,(b*Om2[l[1],l[2]])/(1-b*Om[l[1],l[2]]))
			end
		elseif ranking_measure == "dKf2"
			Om1 = res_dist(L)
			Om2 = res_dist(L,2)
			Om3 = res_dist(L,3)
			for i in 1:m
				l = (Int(Bsp[2*i-1,1]),Int(Bsp[2*i,1]))
				push!(ll,l)
				b = -L[l[1],l[2]]
				O1 = Om1[l[1],l[2]]
				O2 = Om2[l[1],l[2]]
				O3 = Om3[l[1],l[2]]
				push!(mes,(2*b*O3)/(1-b*O1)+(b*O2)^2/(1-b*O1)^2)
			end
		elseif ranking_measure == "load"
			dx1 = x1[1:n]*ones(1,n) - ones(n)*transpose(x1[1:n])
			load = L.*dx1
			load_rate = abs.(dx1./(pi/2))
			for i in 1:m
				l = (Int(Bsp[2*i-1,1]),Int(Bsp[2*i,1]))
				push!(ll,l)
				push!(mes,load_rate[l[1],l[2]])
			end
		elseif ranking_measure == "Omega+load"
			Om = res_dist(L)
			dx1 = x1[1:n]*ones(1,n) - ones(n)*transpose(x1[1:n])
			load = L.*dx1
			for i in 1:m
				l = (Int(Bsp[2*i-1,1]),Int(Bsp[2*i,1]))
				push!(ll,l)
				push!(mes,Om[l[1],l[2]]*load[l[1],l[2]]^2)
			end
		else
			for i in 1:m
				l = (Int(Bsp[2*i-1,1]),Int(Bsp[2*i,1]))
				push!(ll,l)
				@info "$(now()) -- "*ntw*": No clear ranking strategy..."
			end
		end
		
		rank0 = Array{Int64,1}(sortslices([mes 1:m],dims=1,rev=true)[:,2])
		rank = copy(rank0)
		ranks = Array{Array{Int64,1},1}([rank,])
		run = true
		rmvd = Array{Int64,1}()
		cut = Array{Int64,1}()
		for i in 1:m
			Lt = copy(L)
			l = ll[i]
			Lt = rm_line(L,ll[i])
			if !isconnected(Lt)
				push!(cut,i)
			end
		end
		cuts = Array{Array{Int64,1},1}([cut,])
		Lr = copy(L)
		Br = copy(B)
		wr = copy(w)
	end
		
	count = 0
		
	while run && count < m-n+1
#		global count,run,cuts,rmvd,ranks,Lr,rank0,ranks
		count += 1
		@info "$(now()) -- "*ntw*": round $count"
		
		used_rank = Array{Int64,1}()
		if ranking_type == "init"
			used_rank = rank0
		elseif ranking_type == "updated"
			used_rank = rank
		elseif ranking_type == "random"
			used_rank = Array{Int64,1}(sortslices([rand(length(rank)) 1:m],dims=1)[:,2])
		elseif ranking_type == "tini"
			used_rank = rank0[end:-1:1]
		elseif ranking_type == "detadpu"
			used_rank = rank[end:-1:1]
		end
		k = 1
		kmax = length(setdiff(used_rank,cuts[end],rmvd))
		i = setdiff(used_rank,cuts[end],rmvd)[k]
		Lt = rm_line(Lr,ll[i])
		cut = copy(cuts[end])
		while (!isconnected(Lt)) && (k < kmax)
			push!(cut,i)
			k += 1
			i = setdiff(used_rank,cuts[end],rmvd)[k]
			Lt = rm_line(Lr,ll[i])
		end
		
		if k >= kmax
			@info "$(now()) -- "*ntw*": Network splits => EXIT"
			run = false
		else
			Lr = copy(Lt)
			Br,wr = rm_line_inc(Br,wr,ll[i])
			push!(rmvd,i)
			push!(cuts,cut)
			l = ll[i]
			
			xs,n_iter = NR_kuramoto(Br,wr,P)
#			xs,dxs,n_iter = kuramoto2(Lr,M,D,P,x1[1:n],x1[(n+1):(2*n)])
			stable = isstable(Lr,xs)
			
			
			mes = Array{Float64,1}()
			if ranking_measure == "Omega"
				Om = res_dist(Lr)
				for i in setdiff(1:m,rmvd)
					l = ll[i]
					push!(mes,-Lr[l[1],l[2]]*Om[l[1],l[2]])
				end
			elseif ranking_measure == "Womega"
				cdx = cos.(xs[1:n]*ones(1,n)-ones(n)*transpose(xs[1:n])).*(1 .- diagm(0 => ones(n)))
				Lw = Lr.*cdx
				Lw = diagm(0 => vec(sum(Lw,dims=2)))
				Om = res_dist(Lw)
				for i in 1:m
					l = ll[i]
					push!(mes,-Lw[l[1],l[2]]*Om[l[1],l[2]])
				end
			elseif ranking_measure == "dKf1"
				Om = res_dist(Lr)
				Om2 = res_dist(Lr,2)
				for i in setdiff(1:m,rmvd)
					l = ll[i]
					b = -Lr[l[1],l[2]]
					push!(mes,(b*Om2[l[1],l[2]])/(1-b*Om[l[1],l[2]]))
				end
			elseif ranking_measure == "dKf2"
				Om1 = res_dist(Lr)
				Om2 = res_dist(Lr,2)
				Om3 = res_dist(Lr,3)
				for i in 1:m
					l = ll[i]
					b = -Lr[l[1],l[2]]
					O1 = Om1[l[1],l[2]]
					O2 = Om2[l[1],l[2]]
					O3 = Om3[l[1],l[2]]
					push!(mes,(2*b*O3)/(1-b*O1)+(b*O2)^2/(1-b*O1)^2)
				end
			elseif ranking_measure == "load"
				ddx = xs[1:n]*ones(1,n) - ones(n)*transpose(xs[1:n])	
				load = Lr.*ddx
				load_rate = abs.(ddx./(pi/2))
				for i in setdiff(1:m,rmvd)
					l = ll[i]
					push!(mes,load_rate[l[1],l[2]])
				end
			elseif ranking_measure == "Omega+load"
				Om = res_dist(Lr)
				ddx = xs[1:n]*ones(1,n) - ones(n)*transpose(xs[1:n])	
				load = Lr.*ddx
				for i in setdiff(1:m,rmvd)
					l = ll[i]
					push!(mes,Om[l[1],l[2]]*load[l[1],l[2]]^2)
				end
			else 
				for i in setdiff(1:m,rmvd)
					l = ll[i]
					@info "$(now()) -- "*ntw*": No clear ranking strategy..."
				end
			end
			rank = Array{Int64,1}(sortslices([mes setdiff(1:m,rmvd)],dims=1,rev=true)[:,2])
			push!(ranks,rank)
			
			if n_iter >= max_iter
				run = false
				@info "$(now()) -- "*ntw*": No sync anymore => EXIT"
			elseif !stable
				run = false
				@info "$(now()) -- "*ntw*": Solution is not stable => EXIT"
			end
		end
	end
	
	writedlm("data/ranks_1b1_"*ntw*"_$(P0)_"*ranking_type*"_"*ranking_measure*".csv",ranks,',')
	writedlm("data/rmvd_1b1_"*ntw*"_$(P0)_"*ranking_type*"_"*ranking_measure*".csv",rmvd,',')
	writedlm("data/cuts_1b1_"*ntw*"_$(P0)_"*ranking_type*"_"*ranking_measure*".csv",cuts,',')
	
	return ranks,rmvd,cuts
end


function rmv_1b1_2(ntw::String,ranking_type::String,ranking_measure::String,P0::Float64,M::Array{Float64,1},D::Array{Float64,1},max_iter::Int64=50,eps::Float64=1e-6,h::Float64=.1)
	
	Lsp = readdlm(ntw*"_data/"*ntw*"_lap_mat.csv",',')
	Bsp = readdlm(ntw*"_data/"*ntw*"_inc_mat.csv",',')
	
	L = sparse(Array{Int,1}(Lsp[:,1]),Array{Int,1}(Lsp[:,2]),vec(Lsp[:,3]))
	B = sparse(Array{Int,1}(Bsp[:,1]),Array{Int,1}(Bsp[:,2]),vec(Bsp[:,3]))
	w = Bsp[2*(1:Int(size(Bsp)[1]/2)),4]

	n,m = size(B)
		
	P = P0*vec(readdlm(ntw*"_data/P_"*ntw*".csv",','))
	P .-= mean(P)
	
	list_sync = readdir("sync_states")
	if length(intersect(list_sync,[ntw*"_sync_$P0.csv",])) == 0
		x1 = sync(ntw,P0,M,D)
	else
		x1 = vec(readdlm("sync_states/"*ntw*"_sync_$P0.csv",','))
	end
	
	if x1[1] == "nope"
		global rmvd,run
		@info "$(now()) -- "*ntw*": No sync possible!"
		run = false
	else
		ll = Array{Tuple{Int64,Int64},1}()
		mes = Array{Float64,1}()
		
		if ranking_measure == "Omega"
			Om = res_dist(L)
			for i in 1:m
				l = (Int(Bsp[2*i-1,1]),Int(Bsp[2*i,1]))
				push!(ll,l)
				push!(mes,-L[l[1],l[2]]*Om[l[1],l[2]])
			end
		elseif ranking_measure == "Womega"
			cdx = cos.(x1[1:n]*ones(1,n)-ones(n)*transpose(x1[1:n])).*(1 .- diagm(0 => ones(n)))
			Lw = L.*cdx
			Lw = diagm(0 => vec(sum(Lw,dims=2)))
			Om = res_dist(Lw)
			for i in 1:m
				l = (Int(Bsp[2*i-1,1]),Int(Bsp[2*i,1]))
				push!(ll,l)
				push!(mes,-Lw[l[1],l[2]]*Om[l[1],l[2]])
			end
		elseif ranking_measure == "dKf1"
			Om = res_dist(L)
			Om2 = res_dist(L,2)
			for i in 1:m
				l = (Int(Bsp[2*i-1,1]),Int(Bsp[2*i,1]))
				push!(ll,l)
				b = -L[l[1],l[2]]
				push!(mes,(b*Om2[l[1],l[2]])/(1-b*Om[l[1],l[2]]))
			end
		elseif ranking_measure == "dKf2"
			Om1 = res_dist(L)
			Om2 = res_dist(L,2)
			Om3 = res_dist(L,3)
			for i in 1:m
				l = (Int(Bsp[2*i-1,1]),Int(Bsp[2*i,1]))
				push!(ll,l)
				b = -L[l[1],l[2]]
				O1 = Om1[l[1],l[2]]
				O2 = Om2[l[1],l[2]]
				O3 = Om3[l[1],l[2]]
				push!(mes,(2*b*O3)/(1-b*O1)+(b*O2)^2/(1-b*O1)^2)
			end
		elseif ranking_measure == "load"
			dx1 = x1[1:n]*ones(1,n) - ones(n)*transpose(x1[1:n])
			load = L.*dx1
			load_rate = abs.(dx1./(pi/2))
			for i in 1:m
				l = (Int(Bsp[2*i-1,1]),Int(Bsp[2*i,1]))
				push!(ll,l)
				push!(mes,load_rate[l[1],l[2]])
			end
		elseif ranking_measure == "Omega+load"
			Om = res_dist(L)
			dx1 = x1[1:n]*ones(1,n) - ones(n)*transpose(x1[1:n])
			load = L.*dx1
			for i in 1:m
				l = (Int(Bsp[2*i-1,1]),Int(Bsp[2*i,1]))
				push!(ll,l)
				push!(mes,Om[l[1],l[2]]*load[l[1],l[2]]^2)
			end
		else
			for i in 1:m
				l = (Int(Bsp[2*i-1,1]),Int(Bsp[2*i,1]))
				push!(ll,l)
				@info "$(now()) -- "*ntw*": No clear ranking strategy..."
			end
		end
		
		rank0 = Array{Int64,1}(sortslices([mes 1:m],dims=1,rev=true)[:,2])
		rank = copy(rank0)
		run = true
		rmvd = Array{Int64,1}()
		cut = Array{Int64,1}()
		for i in 1:m
			Lt = copy(L)
			l = ll[i]
			Lt = rm_line(L,ll[i])
			if !isconnected(Lt)
				push!(cut,i)
			end
		end
		Lr = copy(L)
		Br = copy(B)
		wr = copy(w)
	end
		
	count = 0
		
	while run && count < m-n+1
		count += 1
		@info "$(now()) -- "*ntw*": round $count"
		
		used_rank = Array{Int64,1}()
		if ranking_type == "init"
			used_rank = rank0
		elseif ranking_type == "updated"
			used_rank = rank
		elseif ranking_type == "random"
			used_rank = Array{Int64,1}(sortslices([rand(length(rank)) 1:m],dims=1)[:,2])
		elseif ranking_type == "tini"
			used_rank = rank0[end:-1:1]
		elseif ranking_type == "detadpu"
			used_rank = rank[end:-1:1]
		end
		k = 1
		kmax = length(setdiff(used_rank,cut,rmvd))
		i = setdiff(used_rank,cut,rmvd)[k]
		Lt = rm_line(Lr,ll[i])
		cut0 = copy(cut)
		while (!isconnected(Lt)) && (k < kmax)
			push!(cut,i)
			k += 1
			i = setdiff(used_rank,cut0,rmvd)[k]
			Lt = rm_line(Lr,ll[i])
		end
		
		if k >= kmax
			@info "$(now()) -- "*ntw*": Network splits => EXIT"
			run = false
		else
			Lr = copy(Lt)
			Br,wr = rm_line_inc(Br,wr,ll[i])
			push!(rmvd,i)
			l = ll[i]
			
			xs,n_iter = NR_kuramoto(Br,wr,P)
#			xs,dxs,n_iter = kuramoto2(Lr,M,D,P,x1[1:n],x1[(n+1):(2*n)])
			stable = isstable(Lr,xs)
			
			
			mes = Array{Float64,1}()
			if ranking_measure == "Omega"
				Om = res_dist(Lr)
				for i in setdiff(1:m,rmvd)
					l = ll[i]
					push!(mes,-Lr[l[1],l[2]]*Om[l[1],l[2]])
				end
			elseif ranking_measure == "Womega"
				cdx = cos.(x1[1:n]*ones(1,n)-ones(n)*transpose(x1[1:n])).*(1 .- diagm(0 => ones(n)))
				Lw = Lr.*cdx
				Lw = diagm(0 => vec(sum(Lw,dims=2)))
				Om = res_dist(Lw)
				for i in setdiff(1:m,rmvd)
					l = ll[i]
					push!(mes,-Lw[l[1],l[2]]*Om[l[1],l[2]])
				end
			elseif ranking_measure == "dKf1"
				Om = res_dist(Lr)
				Om2 = res_dist(Lr,2)
				for i in setdiff(1:m,rmvd)
					l = ll[i]
					b = -Lr[l[1],l[2]]
					push!(mes,(b*Om2[l[1],l[2]])/(1-b*Om[l[1],l[2]]))
				end
			elseif ranking_measure == "dKf2"
				Om1 = res_dist(Lr)
				Om2 = res_dist(Lr,2)
				Om3 = res_dist(Lr,3)
				for i in setdiff(1:m,rmvd)
					l = ll[i]
					b = -Lr[l[1],l[2]]
					O1 = Om1[l[1],l[2]]
					O2 = Om2[l[1],l[2]]
					O3 = Om3[l[1],l[2]]
					push!(mes,(2*b*O3)/(1-b*O1)+(b*O2)^2/(1-b*O1)^2)
				end
			elseif ranking_measure == "load"
				ddx = xs[1:n]*ones(1,n) - ones(n)*transpose(xs[1:n])	
				load = Lr.*ddx
				load_rate = abs.(ddx./(pi/2))
				for i in setdiff(1:m,rmvd)
					l = ll[i]
					push!(mes,load_rate[l[1],l[2]])
				end
			elseif ranking_measure == "Omega+load"
				Om = res_dist(Lr)
				ddx = xs[1:n]*ones(1,n) - ones(n)*transpose(xs[1:n])	
				load = Lr.*ddx
				for i in setdiff(1:m,rmvd)
					l = ll[i]
					push!(mes,Om[l[1],l[2]]*load[l[1],l[2]]^2)
				end
			else 
				for i in setdiff(1:m,rmvd)
					l = ll[i]
					@info "$(now()) -- "*ntw*": No clear ranking strategy..."
				end
			end
			rank = Array{Int64,1}(sortslices([mes setdiff(1:m,rmvd)],dims=1,rev=true)[:,2])
			
			if n_iter >= max_iter
				run = false
				@info "$(now()) -- "*ntw*": No sync anymore => EXIT"
			elseif !stable
				run = false
				@info "$(now()) -- "*ntw*": Solution is not stable => EXIT"
			end
		end
	end
	
	writedlm("data/rmvd_1b1_"*ntw*"_$(P0)_"*ranking_type*"_"*ranking_measure*".csv",rmvd,',')
	
	return rmvd
end




