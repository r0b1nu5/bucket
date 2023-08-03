###### NOT WORKING, BECAUSE STLSQ() DOES NOT PARALLELIZE #######################
@info "Warning: does not parallelize..."

using Distributed 

# Define the number of threads
n_thr = 3

if nworkers() < n_thr
	addprocs(n_thr-nworkers())
end

# If file "parallel" does not exist, then create it
if !isdir("parallel/")
	mkdir("parallel/")
end

@everywhere using DataDrivenDiffEq, ModelingToolkit, LinearAlgebra, DataDrivenSparse, LinearAlgebra, PyPlot, Combinatorics, Statistics

# Infers an hypergraph with knowledge of the states 'X' and of the derivatives 'Y'.
#
# INPUTS:
# 'X': time series of the system's state. Each row is the time series of the state of one agent. 
# 'Y': time series of the system's velocity. Each row is the time series of the velocity of on agent. 
# 'ooi': orders of interest. Vector of integers listing the orders of interactions that we analyze. 
# 'dmax': maximal degree to be considered in the Taylor expansion. Typically, dmax=maximum(ooi).
# 'thr_global': Threshold value to decide whether an edge exists or not.
# 'd': Internal dimension of the agents. It is assmued that the components of an agents are consecutive in the state vector.
#
# OUTPUTS:
# 'Ainf': dictionary associating the inferred coefficient of the Taylor expansion to a pair (node,hyperedge). Namely, 'Ainf[(i,h)]' is the coefficient corresponding to the hyperedge 'h' in the Taylor expansion of the dynamics of node 'i'. 
# 'coeff': matrix of coefficents obtained by SINDy. The row indices are the agents' indices and the columns indices are the indices of the monomials of the Taylor series.
# 'idx_o': dictionary of each monomial index in the matrix coeff. To each order 'o' is associated the list of the indices of monomial of order 'o' involving distinct agents.
# 'agents_o': dictionary of the agents involved in each monomial of 'idx_o'. Each element of 'idx_o[o]' is the index of a monomial, and each element of 'agents_o[o]' is the list of agents (their indices) involved in this monomial.

@everywhere function hyper_inf_par(X::Matrix{Float64}, Y::Matrix{Float64}, ooi::Vector{Int64}, dmax::Int64, thr_glob::Float64=.1)
	n,T = size(X)

	if size(X) != size(Y)
		@info "Dimensions of states and derivatives do not match."
		return nothing
	end

	args = Vector{Tuple{Matrix{Float64},Int64,Matrix{Float64},Int64}}()
	for i in 1:n
		push!(args,(X[[i,],:],i,Y,dmax))
	end

	pmap(infer_par,args)

	coeff = readdlm("parallel/coeff_1.csv",',')
	for i in 2:n
		coeff = [coeff readdlm("parallel/coeff_$i.csv")]
	end
	coeff = Matrix(transpose(coeff))

	# Retrieving the results of SINDy and doing the inference by comparing the identified coefficients with the threshold.
	idx_o = Dict{Int64,Vector{Int64}}()
	agents_o = Dict{Int64,Vector{Vector{Int64}}}()
	inf_o = Dict{Int64,Vector{Vector{Int64}}}()
	Ainf = Dict{Int64,Any}()
	Uinf = Dict{Int64,Any}(maximum(ooi) => Vector{Vector{Int64}}())
	for o in sort(ooi,rev=true)
		Ainf[o] = Dict{Tuple{Int64,Vector{Int64}},Float64}() # Inferred hyperedges of order o
		Uinf[o-1] = Vector{Vector{Int64}}() # Uninferrable hyperedges of order o-1
		idx_o[o],agents_o[o] = get_idx_o(o-1,x,prebasis)
		inf_o[o] = Vector{Int64}()
		for k in 1:length(idx_o[o])
			id = idx_o[o][k]
			agents = agents_o[o][k]
			cagents = setdiff(1:n,agents)
			y = coeff[id,cagents]
			ynz = y[Int64.(setdiff((1:length(y)).*(abs.(y) .> thr_glob),[0.,]))]
			yds = Int64.(setdiff(cagents.*(abs.(y) .> thr_glob),[0,]))
			for j in 1:length(yds)
				Ainf[o][(yds[j],sort([yds[j];agents]))] = ynz[j]
				for a in agents
					push!(Uinf[o-1],sort([yds[j];setdiff(agents,[a,])]))
				end
			end
		end
	end

	return Ainf, coeff, idx_o, agents_o
end

@everywhere function infer_par(args::Tuple{Matrix{Float64},Int64,Matrix{Float64},Int64})
	Xi,i,Y,dmax = args
	n,T = size(Y)

	# Setting up the problem
	@variables x[1:n]
	problem = DirectDataDrivenProblem(Xi,Y,name = :HyperInference)

	# Defining the basis of functions to use, i.e., the monomials up to order 'dmax'.
	prebasis = polynomial_basis([x[i] for i in 1:n],dmax)
	basis = Basis(prebasis,[x[i] for i in 1:n])

	res = solve(problem,basis,STLSQ())
#	coeff = res.out.Ξ[1,:,:]
	coeff = Vector(res.out[1].coefficients)

	writedlm("parallel/coeff_$i.csv")
end


# In case we want to infer only one order of hyperedge.
@everywhere function hyper_inf_par(X::Matrix{Float64}, Y::Matrix{Float64}, ooi::Int64, dmax::Int64, thr_glob::Float64=1.)
    return hyper_inf_par(X,Y,[ooi,],dmax,thr_glob)
end

# When the agents have a internal dimension d higher than 1.
@everywhere function hyper_inf_par(X::Matrix{Float64}, Y::Matrix{Float64}, ooi::Union{Int64,Vector{Int64}}, dmax::Int64, d::Int64=1, thr_glob::Float64=.1)
	Ainf, coeff, idx_o, agents_o = hyper_inf_par(X,Y,ooi,dmax,thr_glob)

	Ainf, AAinf = one2dim(Ainf,d)
end


# Transforms Ainf to agents with d-dimensional internal dynamics.
@everywhere function one2dim(Ainf::Dict{Int64,Any}, d::Int64=1)
	ooi = keys(Ainf)
	A = Dict{Int64,Any}()
	AA = Dict{Int64,Any}()
	for o in ooi
		A[o] = Dict{Tuple{Int64,Vector{Int64}},Float64}()
		AA[o] = Dict{Tuple{Int64,Vector{Int64},Vector{Int64}},Float64}() # Description of the keys: 1. index of the source node of the interaction, 2. list of the interacting agents, 3. list of the components interacting.
		for k in keys(Ainf[o])
			i,v = k
			j = ceil(Int64,i/d)
			u = ceil.(Int64,v./d)
			w = (v .- 1).%d .+ 1
			
			if length(union(u)) == length(u)
				A[o][(j,u)] = 1.
				AA[o][(j,u,w)] = Ainf[o][k]
			end
		end
	end

	return A,AA
end

# Retrieves the indices (in 'prebasis') of the monomials of order 'o' in the variables 'x', involving distincts agents.
@everywhere function get_idx_o(o::Int64, x::Symbolics.Arr{Num,1}, prebasis::Vector{Num})
	comb = combinations(1:length(x),o)
	idx = Int64[]
	agents = Vector{Int64}[]
	for c in comb
		mon = get_monomial(x,c)
		m,i = findmax(isequal.(prebasis,mon))
		if m > .1
			push!(idx,i)
			push!(agents,c)
		end
	end
	return idx,agents
end


# Constructing the monomial of variables in 'x' with indices of 'c'.
@everywhere function get_monomial(x::Symbolics.Arr{Num,1}, c::Vector{Int64})
	mon = 1
	for i in c
		mon *= x[i]
	end
	return mon
end

# Constructing all monomials of order 'o' in the variables 'x'. 
@everywhere function get_monomial(x::Symbolics.Arr{Num,1}, o::Int64)
    n = length(x)
    mon = Num[]
    comb = collect(combinations(1:n,o))

    for c in combi
	    push!(mon,get_monomial(x,c))
    end

    return mon
end

# Returns the inferred 2nd-order adjacency matrix
# Components of 'A2t_bool' are 1. iff the inferred coefficient is larger than 'thr'.
# Components of 'A2t_float' are are equal to the coefficient, irrespective of their magnitude.
@everywhere function inferred_adj_2nd(Ainf2::Dict{Tuple{Int64,Vector{Int64}},Float64}, n::Int64, thr::Float64=0.)
	A2t_bool = zeros(n,n)
	A2t_float = zeros(n,n)

	for k in keys(Ainf2)
		i = k[1]
		ic = setdiff(k[2],[i,])
		for p in permutations(ic)
			A2t_bool[i,p[1]] = (abs.(Ainf2[k]) > thr)
			A2t_float[i,p[1]] = Ainf2[k]
		end
	end

	return A2t_bool, A2t_float
end

# Returns the inferred 3rd-order adjacency tensor
# Components of 'A3t_bool' are 1. iff the inferred coefficient is larger than 'thr'.
# Components of 'A3t_float' are are equal to the coefficient, irrespective of their magnitude.
@everywhere function inferred_adj_3rd(Ainf3::Dict{Tuple{Int64,Vector{Int64}},Float64}, n::Int64, thr::Float64=0.)
	A3t_bool = zeros(n,n,n)
	A3t_float = zeros(n,n,n)

	for k in keys(Ainf3)
		i = k[1]
		ic = setdiff(k[2],[i,])
		for p in permutations(ic)
			A3t_bool[i,p[1],p[2]] = (abs.(Ainf3[k]) > thr)
			A3t_float[i,p[1],p[2]] = Ainf3[k]
		end
	end

	return A3t_bool, A3t_float
end

# Returns the inferred 4th-order adjacency tensor
# Components of 'A4t_bool' are 1. iff the inferred coefficient is larger than 'thr'.
# Components of 'A4t_float' are are equal to the coefficient, irrespective of their magnitude.
@everywhere function inferred_adj_4th(Ainf4::Dict{Tuple{Int64,Vector{Int64}},Float64}, n::Int64, thr::Float64=0.)
	A4t_bool = zeros(n,n,n,n)
	A4t_float = zeros(n,n,n,n)

	for k in keys(Ainf4)
		i = k[1]
		ic = setdiff(k[2],[i,])
		for p in permutations(ic)
			A4t_bool[i,p[1],p[2],p[3]] = (abs.(Ainf4[k]) > thr)
			A4t_float[i,p[1],p[2],p[3]] = Ainf4[k]
		end
	end

	return A4t_bool, A4t_float
end

# Checks the sensitivity and specificity of our inference.
# Considers unweighted interactions, i.e., summarizes the inference in a boolean array.
@everywhere function check_inference_bool(A2::Matrix{Float64}, A3::Array{Float64,3}, A4::Array{Float64,4}, Ainf::Dict{Int64,Any}, thr::Float64=0.)
    n = size(A2)[2]
    A2t = zeros(n,n)
    A3t = zeros(n,n,n)
    A4t = zeros(n,n,n,n)

    # Remove edges belonging to higher-order edges
    A2s34 = A2.*(1 .- ((sum(A3,dims=3)[:,:,1] + sum(A4,dims=[3,4])[:,:,1,1]) .> .1))
    A3s4 = A3.*(1 .- ((sum(A4,dims=4)[:,:,:,1]) .> .1))

    # Compute sensitivity and specificity for the inference of the 2nd order edges (if any).
    if 2 in keys(Ainf)
	    @info "================================================================"
	    @info "2nd-order satistics"
	    A2t = inferred_adj_2nd(Ainf[2],n,thr)[1]
	
	    p2 = 0 # positives
	    fp2 = 0 # false positives
	    n2 = 0 # negatives
	    fn2 = 0 # false negatives
	    for i in 1:n
	        for j in 1:n
	            if A2t[i,j] > .1
	                p2 += 1
#	                if A2[i,j] < .1
			if A2s34[i,j] < .1
	                    fp2 += 1
	                end
	            else
	                n2 += 1
#	                if A2[i,j] > .1
	                if A2s34[i,j] > .1
	                    fn2 += 1
	                end
	            end
	        end
	    end
	    tp2 = p2 - fp2 # true positives
	    tn2 = n2 - fn2 # true negatives
	
	    sen2 = tp2/(tp2 + fn2) # sensitivity
	    spe2 = tn2/(tn2 + fp2) # specificity
	
	    @info "Detection of 2-edges: sensitivity = $(round(sen2,digits=2)), specificity = $(round(spe2,digits=2))."
    else
	    (sen2,spe2,tp2,fp2,tn2,fn2) = (NaN,NaN,NaN,NaN,NaN,NaN)
    end

    # Compute sensitivity and specificity for the inference of the 3rd order edges (if any).
    if 3 in keys(Ainf)
	    @info "================================================================"
	    @info "3rd-order statistics"
	    A3t = inferred_adj_3rd(Ainf[3],n,thr)[1]
	    
	    p3 = 0 # positives
	    fp3 = 0 # false positives
	    n3 = 0 # negatives
	    fn3 = 0 # false negatives
	    for i in 1:n
	        for j in 1:n
	            for k in 1:n
	                if A3t[i,j,k] > .1
	                    p3 += 1
#	                    if A3[i,j,k] < .1
	                    if A3s4[i,j,k] < .1
	                        fp3 += 1
	                    end
	                else
	                    n3 += 1
#	                    if A3[i,j,k] > .1
	                    if A3s4[i,j,k] > .1
	                        fn3 += 1
	                    end
	                end
	            end
	        end
	    end
	    tp3 = p3 - fp3 # true positives
	    tn3 = n3 - fn3 # true negatives
	
	    sen3 = tp3/(tp3 + fn3) # sensitivity
	    spe3 = tn3/(tn3 + fp3) # specificity
	
	    @info "Detection of 3-edges: sensitivity = $(round(sen3,digits=2)), specificity = $(round(spe3,digits=2))."
    else
	    (sen3,spe3,tp3,fp3,tn3,fn3) = (NaN,NaN,NaN,NaN,NaN,NaN)
    end

    # Compute sensitivity and specificity for the inference of the 4th order edges (if any).
    if 4 in keys(Ainf)
	    @info "================================================================"
	    @info "4th-order statistics"
	    A4t = inferred_adj_4th(Ainf[4],n,thr)[1]
	    
	    p4 = 0 # positives
	    fp4 = 0 # false positives
	    n4 = 0 # negatives
	    fn4 = 0 # false negatives
	    for i in 1:n
		    for j in 1:n
			    for k in 1:n
				    for l in 1:n
					    if A4t[i,j,k,l] > .1
						    p4 += 1
						    if A4[i,j,k,l] < .1
							    fp4 += 1
						    end
					    else
						    n4 += 1
						    if A4[i,j,k,l] > .1
							    fn4 += 1
						    end
					    end
				    end
			    end
		    end
	    end

	    tp4 = p4 - fp4 # true positives
	    tn4 = n4 - fn4 # true negatives
	
	    sen4 = tp4/(tp4 + fn4) # sensitivity
	    spe4 = tn4/(tn4 + fp4) # specificity
	
	    @info "Detection of 4-edges: sensitivity = $(round(sen4,digits=2)), specificity = $(round(spe4,digits=2))."
    else
	    (sen4,spe4,tp4,fp4,tn4,fn4) = (NaN,NaN,NaN,NaN,NaN,NaN)
    end
    return (sen2,spe2,tp2,fp2,tn2,fn2), (sen3,spe3,tp3,fp3,tn3,fn3), (sen4,spe4,tp4,fp4,tn4,fn4)
end


# Checks the correlation between the non-zero elements of the inferred and actual adjacency tensors. 
# Each component is considered in the correlation if it is nonzero in at least one of the two tensors (actual or inferred). 
@everywhere function check_inference_float(A2::Matrix{Float64}, A3::Array{Float64,3}, Ainf::Dict{Int64,Any},thr::Float64=1e-6)
	n = size(A2)[2]

	if 2 in keys(iAinf)
		A2t = inferred_adj_2nd(Ainf,n,thr)[2]
		a2 = Float64[]
		a2t = Float64[]
		for i in 1:n
			for j in 1:n
				if min(abs(A2t[i,j]),abs(A2[i,j])) > thr
					push!(a2,A2[i,j])
					push!(a2t,A2t[i,j])
				end
			end
		end
		
		if length(a2) == 0
			@info "No 2-edge inferred (default r2 = 0.0)."
			r2 = 0.
		else
			r2 = cor(a2,a2t)
		end
	end

	if 3 in keys(inf_o)
		A3t = inferred_adj_3rd(Ainf,n,thr)[2]
		a3 = Float64[]
		a3t = Float64[]
		for i in 1:n
			for j in 1:n
				for k in 1:n
					if min(abs(A3t[i,j,k]),abs(A3[i,j,k])) > thr
						push!(a3,A3[i,j,k])
						push!(a3t,A3t[i,j,k])
					end
				end
			end
		end

		if length(a3) == 0
			@info "No 3-edge inferred (default r3 = 0.0)."
			r3 = 0.
		else
			r3 = cor(a3,a3t)
		end
	end

	return (r2,a2,a2t), (r3,a3,a3t)
end





