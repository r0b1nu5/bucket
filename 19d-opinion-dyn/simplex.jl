

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

	x = vec(sum(V,dims=2))

	return x
end

# Slides the point x (typically in the simplex) towards the centroid of the unitary p-simplex. 
# For α = 1, the summit does not move, and for α = 0, the summit reaches the centroid.

function slide_summit(x::Vector{Float64}, α::Float64)
	if α < 0.
		@info "Warning: α < 0."
	elseif α > 1.
		@info "Warning: α > 1."
	end

	p = length(x)
	c = 1/p*ones(p)

	return (1-α)*c + α*x
end

# Slides each column of xs towards the centroid of the unitary p-simplex by a factor α. 

function slide_summit(xs::Matrix{Float64}, α::Float64)
	p,n = size(xs)

	S = zeros(p,0)

	for i in 1:n
		x = xs[:,i]
		S = [S slide_summit(x,α)]
	end

	return S
end

# Slides the i-th column of xs towards the centroid of the unitary p-simplex by a factor αs[i]. 

function slide_summit(xs::Matrix{Float64}, αs::Vector{Float64})
	p,n = size(xs)
	S = zeros(p,0)

	for i in 1:n
		S = [S slide_summit(xs[:,i],αs[i])]
	end

	return S
end

# Slides the q-th summit of the unitary p-simplex towards its centroid by a factor α. 

function slide_summit(p::Int64, q::Int64, α::Float64)

	s = [zeros(q-1);1.;zeros(p-q)]
	
	return slide_summit(s,α)
end

# Slides each of the qs summits towards the centroid by a factor α.

function slide_summit(p::Int64, qs::Vector{Int64}, α::Float64)
	S = zeros(p,0)

	for q in qs
		S = [S slide_summit(p,q,α)]
	end

	return S
end

# Slides each of the qs summits towards the centroid by a factor αs[i] respectively.

function slide_summit(p::Int64, qs::Vector{Int64}, αs::Vector{Float64})
	S = zeros(p,0)

	for i in 1:length(qs)
		S = [S slide_summit(p,qs[i],αs[i])]
	end

	return S
end

# Selects a simplex in Ss at random with probablilities ps.
# ps is assumed to sum to one.

function rand_simplex(Ss::Vector{Matrix{Float64}}, ps::Vector{Float64})
	P = [0;cumsum(ps)]

	m,i = findmax((P.<rand()).*(1:length(P)))

	return Ss[i]
end

# Selects a simplex in Ss uniformly.

function rand_simplex(Ss::Vector{Matrix{Float64}})
	n = length(Ss)

	return rand_simplex(Ss,ones(n)./n)
end


# Defines the summits of the admissible opinion space.
# p: number of parties

function admissible_summits(p::Int64)
	# Generate the dictionnary of the vertices of the unitary simplex.
	V = vertices(p)

	c = "" # Name of the centroid
	for k in 1:p
		c *= "$k"
	end

	# Slide each of the summits to the appropriate position to equalize the volume of each party.
	for k in 1:p
		# Party 'k' has exactly 'binomial(p-1,k-1)' admissible orthoschemes. Therefore, extremists have only one. The ratio between the admissible volumes is the given by 'r', which is the reduction factor for each volume.
		r = 1/binomial(p-1,k-1) 
		V["$k"] = slide_summit(V["$k"],r)
	end

	return V
end

# Generates a dictionnary of the admissible simplices names, labelled by the winning party. 
#=
function admissible_simplices(p::Int64)
	s = Dict{Int64,Vector{String}}()
	for k in 1:p
		s[k] = "$k"
		# TODO


end
=#
	
# Draws one opinion uniformly in the admissible space.
# The list of summits is already given in 'V'. One can get 'V' by using the function 'admissible_summits'.

function rand_opinion(p::Int64, V::Dict{String,Vector{Float64}})
	k = rand(1:p) # Selects one of the parties
	perm = gen_rand_perm([zeros(k-1);ones(p-k)]) # Take one of the permutations of k-1 zeros and p-k ones.
	
	idx = zeros(Int64,p)
	left = Vector{Int64}(setdiff((1 .- perm).*(2:p),[0.,])) # Indices of the parties at the left of party k in the order of the opinion.
	right = Vector{Int64}(setdiff(perm.*(2:p),[0.,])) # Indices of the parties at the right of party k in the order of the opinions.
	idx[1] = k
	idx[left] = k-1:-1:1
	idx[right] = k+1:p

	# Get the summits of the orthoscheme as summits of the simplex where the opinion is to be drawn.
	S = zeros(p,0)
	str = ""
	for j in 1:p
		str *= "$(idx[j])"
		str = prod(sort([str[i] for i in 1:length(str)]))
		S = [S V[str]]
	end

	# Generate the opinion.
	x = unif_simplex_arb(S)

	return x
end

# Generates a random permutation of 1:n elements

function gen_rand_perm(n::Int64)
	p = Vector(1:n)
	for i in 1:n-1
		j = i+rand(1:n-i)
		p = [p[1:i-1];p[j];p[i+1:j-1];p[i];p[j+1:n]]
	end
	return p
end

# Generates a random permutation of x.

function gen_rand_perm(x::Vector{Float64})
	n = length(x)
	p = copy(x)
	for i in 1:n-1
		j = i+rand(0:n-i)
		if j != i
			p = [p[1:i-1];p[j];p[i+1:j-1];p[i];p[j+1:n]]
		end
	end
	return p
end








