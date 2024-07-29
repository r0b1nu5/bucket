using DataDrivenDiffEq, ModelingToolkit, DataDrivenSparse, LinearAlgebra, PyPlot, Combinatorics, Statistics, Random

include("this.jl")

# Infers an hypergraph with knowledge of the states 'X' and of the derivatives 'Y'.
#
# INPUTS:
# 'X': time series of the system's state. Each row is the time series of the state of one agent. 
# 'Y': time series of the system's velocity. Each row is the time series of the velocity of on agent. 
# 'ooi': orders of interest. Vector of integers listing the orders of interactions that we analyze. 
# 'dmax': maximal degree to be considered in the Taylor expansion. Typically, dmax=maximum(ooi).
# 'sparse_thr': Sparsity threshold in SINDy. Threshold value to decide whether an edge exists or not.
# 'd': Internal dimension of the agents. It is assmued that the components of an agents are consecutive in the state vector.
#
# OUTPUTS:
# 'Ainf': dictionary associating the inferred coefficient of the Taylor expansion to a pair (node,hyperedge). Namely, 'Ainf[(i,h)]' is the coefficient corresponding to the hyperedge 'h' in the Taylor expansion of the dynamics of node 'i'. 
# 'coeff': matrix of coefficents obtained by SINDy. The row indices are the agents' indices and the columns indices are the indices of the monomials of the Taylor series.
# 'idx_o': dictionary of each monomial index in the matrix coeff. To each order 'o' is associated the list of the indices of monomial of order 'o' involving distinct agents.
# 'agents_o': dictionary of the agents involved in each monomial of 'idx_o'. Each element of 'idx_o[o]' is the index of a monomial, and each element of 'agents_o[o]' is the list of agents (their indices) involved in this monomial.

function hyper_inf(X::Matrix{Float64}, Y::Matrix{Float64}, ooi::Vector{Int64}, dmax::Int64, λ::Float64=.1, ρ::Float64=1., niter::Int64=10)
	n,T = size(X)
	zer0 = 1e-10

	if size(X) != size(Y)
		@info "Dimensions of states and derivatives do not match."
		return nothing
	end

 #= ####### THIS IS DONE IN THIS FOR NOW...
	# Defining the basis of functions to use, i.e., the monomials up to order 'dmax'.
	@variables x[1:n]
	@info "Prebasis..."
	prebasis = polynomial_basis([x[i] for i in 1:n],dmax)
	@info "Basis..."
	basis = Basis(prebasis,[x[i] for i in 1:n])
	l = length(basis.eqs)
# =#

#= # USES PRE-IMPLEMENTED SINDY
	# Setting up the problem
	problem = DirectDataDrivenProblem(X,Y,name = :HyperInference)

	@info "Problem is set."

	# Solving the problem using SINDy.
#=
	coeff = try 
		res = solve(problem,basis,STLSQ())
#		res.out.Ξ[1,:,:]
		Matrix(res.out[1].coefficients')
	catch e
		if isa(e,DimensionMismatch)
			@error "No interaction was inferred for some of the variables. Either some of them are completely disconnected from the rest of the system (in which case they need to be removed from the data), or the time series were too far from eachother and no Taylor expansion was valid (in which case, the spread of initial conditions should be reduced)."
			zeros(n,l)
		end
	end
=#

	res = solve(problem,basis,STLSQ(sparse_thr))
#	res = solve(problem,basis,SR3()); @info "Finished SR3."
#	res = solve(problem,basis,ADMM()); @info "Finished ADMM."

	@info "Problem is solved."

# Apparently, depending on some package version, either of the following lines can work. Choose your own and comment the other.
#	coeff = Matrix(res.out.Ξ[1,:,:]')
	coeff = Matrix(res.out[1].coefficients)
	err = res.residuals

#	@info "coeff = $coeff"
# =#

 #= RESTRICT THE MONOMIALS TO LARGE CORRELATION
	forbid = keep_correlated(X,.5)
	coeff,idx_mon,err,relerr = this(X,Y,ooi,dmax,forbid,λ,ρ,niter)

# =#

# #= USES MYSINDY
	coeff,idx_mon,err,relerr = this(X,Y,ooi,dmax,λ,ρ,niter)
# =#

# #=
	Ainf = Dict{Int64,Matrix{Float64}}(o => zeros(0,o+1) for o in 1:dmax+1)
	for id in keys(idx_mon)
		js = idx_mon[id]
		for i in setdiff(1:n,js)
			Ainf[length(js)+1] = [Ainf[length(js)+1];[i js' coeff[i,id]]]
		end
	end
# =#

 #=
	Ainf = Dict{Int64,Matrix{Float64}}(o => zeros(0,o+1) for o in 1:dmax+1)
	edges = Dict{Vector{Float64},Int64}()
	for i in 1:n
		ids = (1:length(coeff[i,:]))[abs.(coeff[i,:]) .> zer0]
		for id in ids
			js = idx_mon[id]
			edge = [i;sort(setdiff(js,[i,]))]
			o = length(edge)
			if edge in keys(edges)
				Ainf[o][edges[edge],o+1] += coeff[i,id]
			else
				Ainf[o] = [Ainf[o];[edge' coeff[i,id]]]
				edges[edge] = size(Ainf[o])[1]
			end
		end
	end
# =#

 #= ################## A BIT OUTDATED NOW...
	# Retrieving the results of SINDy and doing the inference by comparing the identified coefficients with the threshold.
	idx_o = Dict{Int64,Vector{Int64}}()
	agents_o = Dict{Int64,Vector{Vector{Int64}}}()
	Ainf_ = Dict{Int64,Any}()
	Uinf = Dict{Int64,Any}(maximum(ooi) => Vector{Vector{Int64}}())
	for o in sort(ooi,rev=true)
		Ainf_[o] = Dict{Tuple{Int64,Vector{Int64}},Float64}() # Inferred hyperedges of order o
		Uinf[o-1] = Vector{Vector{Int64}}() # Uninferrable hyperedges of order o-1
		idx_o[o],agents_o[o] = get_idx_o(o-1,x,prebasis)
		for k in 1:length(idx_o[o])
			id = idx_o[o][k]
			agents = agents_o[o][k]
			cagents = setdiff(1:n,agents)
			y = coeff[cagents,id]
			ynz = y[Int64.(setdiff((1:length(y)),[0.,]))]
			yds = Int64.(setdiff(cagents,[0,]))
			for j in 1:length(yds)
				Ainf_[o][(yds[j],sort([yds[j];agents]))] = ynz[j]
				for a in agents
					push!(Uinf[o-1],sort([yds[j];setdiff(agents,[a,])]))
				end
			end
		end
	end
# =#

##	return Ainf, coeff, idx_o, agents_o
	return Ainf, coeff, err, relerr
#	return Ainf, Ainf_, coeff, err, relerr
end

# Same with the tuning of the sparsity parameter λ in SINDy
function hyper_inf_sparsity(X::Matrix{Float64}, Y::Matrix{Float64}, ooi::Vector{Int64}, dmax::Int64, thr_glob::Float64=.1, λ::Float64=1e-1)
	n,T = size(X)

	if size(X) != size(Y)
		@info "Dimensions of states and derivatives do not match."
		return nothing
	end

	# Setting up the problem
	@variables x[1:n]
	problem = DirectDataDrivenProblem(X,Y,name = :HyperInference)

	# Defining the basis of functions to use, i.e., the monomials up to order 'dmax'.
	prebasis = polynomial_basis([x[i] for i in 1:n],dmax)
	basis = Basis(prebasis,[x[i] for i in 1:n])
	l = length(basis.eqs)

	@info "Problem is set."

	res = solve(problem,basis,STLSQ(λ)); @info "Finished SINDy."

	@info "Problem is solved."

# Apparently, depending on some package version, either of the following lines can work. Choose your own and comment the other.
#	coeff = Matrix(res.out.Ξ[1,:,:]')
	coeff = Matrix(res.out[1].coefficients)

	# Retrieving the results of SINDy and doing the inference by comparing the identified coefficients with the threshold.
	idx_o = Dict{Int64,Vector{Int64}}()
	agents_o = Dict{Int64,Vector{Vector{Int64}}}()
	Ainf = Dict{Int64,Any}()
	Uinf = Dict{Int64,Any}(maximum(ooi) => Vector{Vector{Int64}}())
	for o in sort(ooi,rev=true)
		Ainf[o] = Dict{Tuple{Int64,Vector{Int64}},Float64}() # Inferred hyperedges of order o
		Uinf[o-1] = Vector{Vector{Int64}}() # Uninferrable hyperedges of order o-1
		idx_o[o],agents_o[o] = get_idx_o(o-1,x,prebasis)
		for k in 1:length(idx_o[o])
			id = idx_o[o][k]
			agents = agents_o[o][k]
			cagents = setdiff(1:n,agents)
			y = coeff[cagents,id]
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

#	return Ainf, coeff, idx_o, agents_o
	return Ainf, coeff, res.residuals
end

# In case we want to infer only one order of hyperedge.
function hyper_inf(X::Matrix{Float64}, Y::Matrix{Float64}, ooi::Int64, dmax::Int64, thr_glob::Float64=1.)
    return hyper_inf(X,Y,[ooi,],dmax,thr_glob)
end

# When the agents have a internal dimension d higher than 1.
function hyper_inf(X::Matrix{Float64}, Y::Matrix{Float64}, ooi::Union{Int64,Vector{Int64}}, dmax::Int64, d::Int64=1, thr_glob::Float64=.1)
	Ainf, coeff, idx_o, agents_o = hyper_inf(X,Y,ooi,dmax,thr_glob)

	Ainf, AAinf = one2dim(Ainf,d)


end

# Transforms Ainf to agents with d-dimensional internal dynamics.
function one2dim(Ainf::Dict{Int64,Any}, d::Int64=1)
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

# Generates the adjacency tensors of systems with d internal dimension of the agents, from the Ainf blindly inferred from 'hyper_inf' without knowledge of the internal dimension of the agents. (Implemented for up to 3-order edges.) Note that the interactions are assumed symetric.
function adj_tensors(Ainf::Dict{Int64,Any}, n::Int64, d::Int64=1)
	A2 = zeros(n,n)
	A3 = zeros(n,n,n)

	if 2 in keys(Ainf)
		for k in keys(Ainf[2])
			i,v = k
			j = ceil(Int64,i/d)
			u = ceil.(Int64,v./d)
			w = (v .- 1).%d .+ 1

			if length(union(u)) == 2
				p,q = u
				A2[p,q] = max(A2[p,q],abs(Ainf[2][k]))
				A2[q,p] = max(A2[q,p],abs(Ainf[2][k]))
			end
		end
	end

	O3 = zeros(3,3,3)
	O3[1,2,3] = 1.
	O3[1,3,2] = 1.
	O3[2,1,3] = 1.
	O3[2,3,1] = 1.
	O3[3,1,2] = 1.
	O3[3,2,1] = 1.
	if 3 in keys(Ainf)
		for k in keys(Ainf[3])
			i,v = k
			j = ceil(Int64,i/d)
			u = ceil.(Int64,v./d)
			w = (v .- 1).%d .+ 1

			if length(union(u)) == 3
				p,q,r = u
				A3[u,u,u] = max(A3[p,q,r],abs(Ainf[3][k]))*O3
			end
		end
	end

	return A2,A3
end

# Retrieves the indices (in 'prebasis') of the monomials of order 'o' in the variables 'x', involving distincts agents.
function get_idx_o(o::Int64, x::Symbolics.Arr{Num,1}, prebasis::Vector{Num})
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
function get_monomial(x::Symbolics.Arr{Num,1}, c::Vector{Int64})
	mon = 1
	for i in c
		mon *= x[i]
	end
	return mon
end

# Constructing all monomials of order 'o' in the variables 'x'. 
function get_monomial(x::Symbolics.Arr{Num,1}, o::Int64)
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
function inferred_adj_2nd(Ainf2::Dict{Tuple{Int64,Vector{Int64}},Float64}, n::Int64, thr::Float64=0.)
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
function inferred_adj_3rd(Ainf3::Dict{Tuple{Int64,Vector{Int64}},Float64}, n::Int64, thr::Float64=0.)
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
function inferred_adj_4th(Ainf4::Dict{Tuple{Int64,Vector{Int64}},Float64}, n::Int64, thr::Float64=0.)
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

# Concatenates A2 and the slices of A3 into a n x (n^2+1) matrix.
function cat_As(A2::Matrix{Float64}, A3::Array{Float64,3})
	n = size(A2)[1]
	
	adj = A2
	for i in 1:n
		adj = [adj A3[:,:,i]]
	end

	return adj
end

# Checks the sensitivity and specificity of our inference.
# Considers unweighted interactions, i.e., summarizes the inference in a boolean array.
function check_inference_bool(A2::Matrix{Float64}, A3::Array{Float64,3}, A4::Array{Float64,4}, Ainf::Dict{Int64,Any}, thr::Float64=0.)
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
function check_inference_float(A2::Matrix{Float64}, A3::Array{Float64,3}, Ainf::Dict{Int64,Any},thr::Float64=1e-6)
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

# ROC curve for adjacenty lists
# A0 is the inferred adjacency, A is the ground truth (assumed boolean adjacency list).
# hyperedges are distinct by their first index, the ordering of the other indices does not matter.
function my_ROC(A0::Matrix{Float64}, A::Matrix{Float64}, n::Int64)
	m,o = size(A0)
	mm,o = size(A)
	o -= 1

	max_edges = n*binomial(n-1,o-1)
	max_edges = n*sum(binomial(n-i,o-i) for i in 1:o)

	A1 = sortslices([A0[:,o+1] A0[:,1:o]],dims=1,rev=true)
	I1 = [A1[i,2:o+1] for i in 1:m]
	V1 = [A1[i,1] for i in 1:m]
	I = [A[i,1:o] for i in 1:mm]

	tp = [0,]
	fp = [0,]
	
	for i in I1
		if i in I
			push!(tp,tp[end]+1)
			push!(fp,fp[end])
		else
			push!(tp,tp[end])
			push!(fp,fp[end]+1)
		end
	end

	# Completes the inference randomly
	mtp = length(I) - tp[end]
	mfp = max_edges-length(I)-fp[end]
	@info "$mtp, $mfp"
	t = shuffle([ones(mtp);zeros(mfp)])
	tt = 1 .- t
	s = [sum(t[1:i]) for i in 1:length(t)]
	ss = [sum(tt[1:i]) for i in 1:length(tt)]

	tp = [tp;(tp[end] .+ s)]
	fp = [fp;(fp[end] .+ ss)]

	#push!(tp,length(I))
	#push!(fp,max_edges-length(I))

	tpr = tp/length(I)
	fpr = fp/(max_edges-length(I))

	return tpr,fpr
end

function my_ROC(A01::Matrix{Float64}, A1::Matrix{Float64}, A02::Matrix{Float64}, A2::Matrix{Float64}, n::Int64)
	m1,o1 = size(A01)
	mm1,o1 = size(A1)
	o1 -= 1
	max_edges1 = n*sum(binomial(n-i,o1-i) for i in 1:o1)

	m2,o2 = size(A02)
	mm2,o2 = size(A2)
	o2 -= 1
	max_edges2 = n*sum(binomial(n-i,o2-i) for i in 1:o2)

	max_edges = max_edges1 + max_edges2

	a1 = sortslices([A01[:,o1+1] A01[:,1:o1]],dims=1,rev=true)
	I01 = [a1[i,2:o1+1] for i in 1:m1]
	V1 = [a1[i,1] for i in 1:m1]
	I1 = [A1[i,1:o1] for i in 1:mm1]

	a2 = sortslices([A02[:,o2+1] A02[:,1:o2]],dims=1,rev=true)
	I02 = [a2[i,2:o2+1] for i in 1:m2]
	V2 = [a2[i,1] for i in 1:m2]
	I2 = [A2[i,1:o2] for i in 1:mm2]

	I0 = [I01;I02]
	V = [V1;V2]
	I = [I1;I2]

	tp = [0,]
	fp = [0,]
	
	for i in I0
		if i in I
			push!(tp,tp[end]+1)
			push!(fp,fp[end])
		else
			push!(tp,tp[end])
			push!(fp,fp[end]+1)
		end
	end

	# Completes the inference randomly
	mtp = length(I) - tp[end]
	mfp = max_edges-length(I)-fp[end]
	@info "$mtp, $mfp"
	t = shuffle([ones(mtp);zeros(mfp)])
	tt = 1 .- t
	s = [sum(t[1:i]) for i in 1:length(t)]
	ss = [sum(tt[1:i]) for i in 1:length(tt)]

	tp = [tp;(tp[end] .+ s)]
	fp = [fp;(fp[end] .+ ss)]

	#push!(tp,length(I))
	#push!(fp,max_edges-length(I))

	tpr = tp/length(I)
	fpr = fp/(max_edges-length(I))

	return tpr,fpr
end

# returns the list of pairs of agents whose correlation is smaller than α.
function keep_correlated(X::Matrix{Float64}, α::Float64=.5)
	C = cor(X')
	n = size(C)[1]
	
	c = [abs(C[i,j]) for i in 1:n-1 for j in i+1:n]
	ids = [[i,j] for i in 1:n-1 for j in i+1:n]
	forbid = ids[c .< α]

	@info "$(n*(n-1)/2 - length(forbid)) pairs kept, $(length(forbid)) pairs discarded."

	return forbid
end


