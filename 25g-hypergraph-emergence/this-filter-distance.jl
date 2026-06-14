include("this-tools.jl")
include("this.jl")

# ================================================================================
function this_filter_distance(X::Matrix{Float64}, Y::Matrix{Float64}, ooi::Vector{Int64}, dmax::Int64, dist_keep::Union{Int64,Float64}, coord::Union{Matrix{Int64},Matrix{Float64}}, λ::Float64=.1, ρ::Float64=1., niter::Int64=10)
	if size(X) != size(Y)
		@info "Dimensions of states and derivatives do not match."
		return nothing
	end
	n,T = size(X)

	# Listing pairs of agents whose trajectories have a correlation coefficient above 'α'.
	keep = keep_distance(coord,dist_keep)

	# Running THIS ##############################################################
	d = get_d(n,dmax) # Lists the indices of the agents involved in each monomial.
	idx_mon = Dict{Int64,Vector{Int64}}()
	for i in 1:size(d)[1]
		mon = d[i,:][d[i,:] .!= 0]
#		if length(mon) == length(union(mon))
			idx_mon[i] = sort(mon)
#		end
	end
	d2i = Dict{Vector{Int64},Int64}(sort(d[i,:]) => i for i in 1:size(d)[1]) # Returns the index of each element of 'd'.
#	d2i = Dict{Vector{Int64},Int64}(idx_mon[k] => k for k in keys(idx_mon)) # Returns the index of each element of 'd'.

	i2keep = Dict{Int64,Vector{Int64}}(i => Int64[] for i in 1:n) 
	for k in 1:length(keep)
		i,j = keep[k]
		push!(i2keep[i],k)
	end

	coeff = Dict{Int64,Matrix{Float64}}()
	ids = Dict{Int64,Vector{Int64}}(i => Int64[] for i in 1:n) # For each agent, contains the indices of the relevant monomials (according to 'keep').
	err = 0.
	energy = 0.
	for i in 1:n
		θ = ones(1,T)
		id = [1,]
		for k in i2keep[i]
			j = setdiff(keep[k],[i,])[1]
			θ = vcat(θ,reduce(vcat,[prod(X[v,:],dims=1) for v in [setdiff([j,l],[0,]) for l in 0:n]]))
			append!(id,[d2i[v] for v in [sort([j,l]) for l in 0:n]])
		end

		iii = unique(a -> id[a], eachindex(id))
		id = id[iii]
		θ = θ[iii,:]

		@info "Running mySINDy for agent $i"
		if isempty(id)
			coef = zeros(0,0)
			er = norm(Y[i,:],1)
		else
			coef,er = mySINDy(θ,Y[[i,],:],λ,ρ,niter)
			ids[i] = id
		end

		coeff[i] = coef
		err += er
		energy += norm(Y[i,:],1)
	end

	relerr = err/energy

	@info "THIS completed."
	#############################################################################
	
	# Reconstructing the adjacency tensors ######################################
	Ainf = Dict{Int64,Matrix{Float64}}(o => zeros(0,o+1) for o in 1:dmax+1) # For each order o, associates a dictionary associating the pair (agents, hyperedge) to the inferred weight, as seen from the agent.
	for i in 1:n
		for j in 1:length(ids[i])
			id = ids[i][j]
			jj = idx_mon[id]
			o = length(jj)+1
			Ainf[o] = vcat(Ainf[o],[i jj' coeff[i][j]])
		end
	end

	@info "Dictionary of inferred adjacency tensors built."
	#############################################################################

	return Ainf, coeff, ids, relerr
end

# ================================================================================

function keep_distance(coord::Union{Matrix{Int64},Matrix{Float64}}, dist::Union{Int64,Float64})
	n,dim = size(coord)

	d = [norm(coord[i,:]-coord[j,:]) for i in 1:n-1 for j in i+1:n]
	ids = [[i,j] for i in 1:n-1 for j in i+1:n]
	keep = ids[d .< dist]

	@info "Filtering: $(length(keep))/$(Int(n*(n-1)/2)) pairs kept."

	return keep
end






