using DelimitedFiles,Statistics,LinearAlgebra,SparseArrays

include("rm_line.jl")
include("isconnected.jl")
include("kuramoto.jl")
include("res_dist.jl")

function rmv_1b1(ntw::String,ranking_type::String,ranking_measure::String,P0::Float64,M::Array{Float64,1},D::Array{Float64,1},eps::Float64=1e-6,max_iter::Int64=100000,h::Float64=.1)
	
	Asp = readdlm(ntw*"_adj_mat.csv",',')
	
	n = Int(maximum(Asp))
	m = Int(size(Asp)[1]/2)
	
	A = sparse(vec(Asp[:,1]),vec(Asp[:,2]),vec(Asp[:,3]))
	L = spdiagm(0 => vec(sum(A,dims=2))) - A
	
	Om = res_dist(L)
	
	P = P0*vec(readdlm("P_"*ntw*".csv"))
	P .-= mean(P)
	
	x1 = vec(readdlm("sync_states/"*ntw*"_sync_$P0.csv",','))
	dx1 = x[1:n]*ones(1,n) - ones(n)*transpose(x1[1:n])
	load = L.*dx1
	load_rate = abs.(dx1./(pi/2))
	
	ll = Array{Tuple{Int64,Int64},1}()
	mes = Array{Float64,1}()
	for i in 1:m
		l = (Int64(Asp[2*i-1,1]),Int(Asp[2*i-1,2]))
		push!(ll,l)
		if ranking_measure == "Omega"
			push!(mes,Om[l[1],l[2]])
		elseif ranking_measure == "load"
			push!(mes,load_rate[l[1],l[2]])
		elseif ranking_measure == "Omega+load"
			push!(mes,Om[l[1],l[2]]*load[l[1],l[2]]^2)
		else
			@info "$(now()) -- "*ntw*": No clear ranking strategy..."
		end
	end
	
	rank = Array{Int64,1}(sortslices([mes 1:m],dims=1,rev=true)[:,2])
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
	
	count = 0
	
	if x1[1] == "nope"
		global cuts,ranks,rmvd,run
		@info "$(now()) -- "*ntw*": No sync possible!"
		cuts = Array{Array{Int64,1},1}()
		ranks = Array{Array{Int64,1},1}()
		run = false
	end
	
	while run && count < m-n+1
		global count,Lr,run,cuts,rmvd,rank,ranks
		count += 1
		@info "$(now()) -- "*ntw*": round $count"
		
		used_rank = Array{Int64,1}()
		if
#		TODO: define used ranking...
		end
		k = 1
		i = setdiff(rank,cuts[end],rmvd)[k]
		Lt = rm_line(Lr,ll[i])
		cut = cuts[end]
		while !isconnected(Lt)
			push!(cut,i)
			k += 1
			i = setdiff(rank,cuts[end],rmvd)[k]
			Lt = rm_line(Lr,ll[i])
		end
		Lr = copy(Lt)
		push!(rmvd,i)
		push!(cuts,cut)
		l = ll[i]
		Om = res_dist(Lr)
		
		xs,dxs,n_iter = kuramoto2(Lr,M,D,P,x1[1:n],x1[(n+1):(2*n)])
		
		ddx = xs[1:n]*ones(1,n) - ones(n)*transpose(xs[1:n])
		
		load = Lr.*ddx
		load_rate = abs.(ddx./(pi/2))
		
		mes = Array{Float64,1}()
		for i in setdiff(1:m,rmvd)
			l = ll[i]
			if ranking_measure == "Omega"
				push!(mes,Om[l[1],l[2]])
			elseif ranking_measure == "load"
				push!(mes,load_rate[l[1],l[2]])
			elseif ranking_measure == "Omega+load"
				push!(mes,Om[l[1],l[2]]*load[l[1],l[2]]^2)
			else 
				@info "$(now()) -- "*ntw*": No clear ranking strategy..."
			end
		end
		rank = Array{Int64,1}(sortslices([mes setdiff(1:m,rmvd)],dims=1,rev=true)[:,2])
		push!(ranks,rank)
		
		if n_iter >= max_iter
			run = false
			@info "$(now()) -- "*ntw*": No sync anymore."
		end
	end
	
	writedlm("data/ranks_1b1_$(P0)_"*ranking_type*"_"*ranking_measure*".csv",ranks,',')
	writedlm("data/rmvd_1b1_$(P0)_"*ranking_type*"_"*ranking_measure*".csv",rmvd,',')
	writedlm("data/cuts_1b1_$(P0)_"*ranking_type*"_"*ranking_measure*".csv",cuts,',')
end



