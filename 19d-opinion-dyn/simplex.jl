

function vertices(p::Int64)
	vertex = Dict{String,Vector{Float64}}()

	if p == 1
		vertex["1"] = [1.,]
	else
		v0 = vertices(p-1)
		ks = keys(v0)
		
		s = [zeros(p-1);1.]
		vertex["$p"] = s

		for k in ks
			vertex[k] = [v0[k];0.]
			l = length(k)
			vertex[k*"$p"] = (l*vertex[k]+s)/(l+1)
		end
	end

	return vertex
end

# x is a point in the rectangle simplex in R^{p-1}, i.e., such that its components sum to a value smaller than or equal to one. 
# x is sent to its counterpart in the simplex defined by the columns of S, which are in R^p. 

function rect2simplex(x::Vector{Float64}, S::Matrix{Float64})
	p = length(x) + 1
	V = (S - repeat(S[:,1],1,p))[:,2:p]

	return S[:,1] + V*x
end


# Draws a random point in the p-simplex uniformly, following the "sorting" method of Donald B. Rubin, The Bayesian bootstrap Ann. Statist. 9, 1981, 130-134.

function unif_simplex_unit(p::Int64)
	x = [0.;sort(rand(p-1));1.]

	return x[2:p+1]-x[1:p]
end


# Draws a random point uniformly from an arbitrary simplex, defined by the columns of S.
# Follows the idea presented in doi.org/10.13140/RG.2.1.3807.6968.

function unif_simplex_arb(S::Matrix{Float64})
	n,p = size(S)

	zs = [1.;rand(p-1);0.]
	λs = [1.,]
	for j in 2:p
		push!(λs,zs[j]^(1/(p-j+1)))
	end
	push!(λs,0.)

	pλs = cumprod(λs)

	V = (1-λs[2])*pλs[1]*S[:,1]
	for i in 2:p
		V = [V (1-λs[i+1])*pλs[i]*S[:,i]]
	end

	x = sum(V,dims=2)

	return x
end








