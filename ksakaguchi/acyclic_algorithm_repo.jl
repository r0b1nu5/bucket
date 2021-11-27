include("tools_repo.jl")

"""
	run_acyclic_algorithm(B::Union{Matrix{Float64},SparseMatrixCSC{Float64,Int64}}, ω::Vector{Float64}, h::Function, γ::Tuple{Float64,Float64})
	run_acyclic_algorithm(B::Union{Matrix{Float64},SparseMatrixCSC{Float64,Int64}}, ω::Vector{Float64}, h::Vector{Function}, γ::Tuple{Float64,Float64})

Runs the algorithm to decide if a solution exist for the Dissipative Flow Problem [Delabays, Jafarpour, and Bullo (2021)] on an acyclic graph. The nodes and edges of the graph are first reordered in order to match requirements described in the paper. Then the original numbering is retrieved.

_INPUT_:\\
`B`: Incidence matrix of the graph considered, with n vertices and m edges. The graph is considered undirected and unweighted (weights can be incorporated in the coupling functions).\\
`ω`: Vector of natural frequencies of the oscillators. \\
`h`: Coupling function(s) and its (their) inverse(s) over the edges of the interaction graph. If a single function is given, the couplings are assumed to identical. Otherwise, a (2m)-vector of functions needs to be given. Edge indexing is such that for e = 1:m, the orientation of e corresponds to the one given by the e-th column of B, and for e = m+1:2*m, the orientation of e is the opposite of edge e-m.\\
`γ`: Tuple(s) of the lower and upper bounds on the argument of the coupling function(s) `h`, such that they are strictly increasing. Dimension is the same as `h`.

_OUTPUT_:\\
`exists`: Returns true only if a solution exists. If the solution exists, then it is unique [Delabays, Jafarpour, and Bullo (2021)].\\
`θ`: Solution (if it exists) of the Dissipative Flow Network problem, defined up to a constant angle shift. 
`ff`: Vector (of dimension 2m) of flows on the edges. For e = 1:m (resp. e = m+1:2*m), ff[e] is the flow over the edge defined by B0[:,e] (resp. -B0[:,e]).\\
`φ`: Synchronous frequency corresponding to the solution. 
"""
function run_acyclic_algorithm(B::Union{Matrix{Float64},SparseMatrixCSC{Float64,Int64}}, ω::Vector{Float64}, h::Tuple{Function,Function}, γ::Tuple{Float64,Float64})
	n,m = size(B)
	
	# For homogeneous couplings, construct a vector of identical transfer functions and give it the vectorial version.
	h1 = Vector{Tuple{Function,Function}}()
	γ1 = Vector{Tuple{Float64,Float64}}()
	for i in 1:2*m
		push!(h1,h)
		push!(γ1,γ)
	end
	
	return run_acyclic_algorithm(B,ω,h1,γ1)
end

function run_acyclic_algorithm(B::Union{Matrix{Float64},SparseMatrixCSC{Float64,Int64}}, ω::Vector{Float64}, h::Vector{Tuple{Function,Function}}, γ::Vector{Tuple{Float64,Float64}})
	n,m = size(B)
	
	# Dimension checks
	if m != n-1
		@info "The graph is not acyclic or not connected: aborting."
		return nothing
	end
	if length(ω) != n
		@info "Laplacian and natural frequencies do not have same dimension ($n,$(length(ω))): aborting."
		return nothing
	end
	if length(h) != 2*m
		@info "There not the same number of coupling functions and edges ($(length(H)),$(2*m)): aborting."
		return nothing
	end

	# Sort nodes and edges in order to satisfy the assumptions of Theorem XXX [Delabays, Jafarpour, and Bullo (2021)].
	B0,id,ed = reindex([B -B])
	B1 = B0[:,1:m]
	ω1 = ω[id]
	h1 = h[ed]
	γ1 = γ[ed]

	# Construct the transfer functions H[ij](f) = h[ij](h[ji]^{-1}(f)) and their domain.
	H1 = Vector{Function}()
	hγ1 = Vector{Tuple{Float64,Float64}}()
	for i in 1:m
		push!(H1,(f -> (h1[i][1](-h1[i+m][2](f)))))
		push!(hγ1,(h1[i][1](γ1[i][1]),h1[i][1](γ1[i][2])))
	end
	for i in m+1:2*m
		push!(H1,(f -> (h1[i][1](-h1[i-m][2](f)))))
		push!(hγ1,(h1[i][1](γ1[i][1]),h1[i][1](γ1[i][2])))
	end


	# Run the algorithm
	exist,f,φ = acyclic_algorithm(B1,H1,ω1,hγ1)


	# De-sort nodes
	if exist
		di,de = deindex(id,ed)
		Bf = B0[di,de]

		ff = f[de]

		Δ = [h[i][2](ff[i]) for i in 1:2*m]

		Bd = pinv(Matrix(B))'
		θ = Bd*Δ[1:m]
	else
		ff = f
		θ = zeros(n)
	end

	return exist, θ, ff, φ
end


"""
	acyclic_algorithm(B::Union{Matrix{Float64},SparseMatrixCSC{Float64,Int64}}, H::Vector{Function}, ω::Vector{Float64}, hγ::Vector{Tuple{Float64,Float64}}, ϵ::Float64=1e-10)
	
Recursive algorithm determining the existence of a unique solution for the Dissipative Flow Problem [Delabays, Jafarpour, and Bullo (2021)] on an acyclic graph. The algorithm requires the nodes and edges to be indexed as described in the original paper.

_INPUT_:\\
`B`: Incidence matrix of the graph considered. The nodes and edges are assumed to ordered according to the requirements of the original reference [Delabays, Jafarpour, and Bullo (2021)]. \\
`H`: Vector (of dimension 2m) of the transfer functions, relating the flows over the two orientations of the edges. For e = 1:m, H[e](f) = h[e](h[e+m]^{-1}(f)) and for e = m+1:2*m, H[e](f) = h[e](h[e-m]^{-1}(f)). \\
`ω`: Vector of natural frequencies of the oscillators.\\
`hγ`: Vector of tuples of the lower and upper bounds of the domain of the transfer functions.\\
`ϵ`: Correction parameter to avoid evaluating the transfer functions exactly at the boundary of their domain.

_OUTPUT_:\\
`exists`: Returns true if a solution exists. If it does, then it is unique [Delabays, Jafarpour, and Bullo (2021)].\\
`ff`: Vector (of dimension 2m) of flows on the edges, solving the Dissipative Flow Problem. For e = 1:m (resp. e = m+1:2*m), ff[e] is the
  flow over the edge defined by B0:,e (resp. -B0[:,e]).\\
`φ`: Synchronous frequency corresponding the solution.
"""
function acyclic_algorithm(B::Union{Matrix{Float64},SparseMatrixCSC{Float64,Int64}}, H::Vector{Function}, ω::Vector{Float64}, hγ::Vector{Tuple{Float64,Float64}}, ϵ::Float64=1e-10)
	# Initializing...
	n,m = size(B)

	keepon = true

	te,et = targets(B)

	l = Vector{Float64}()
	u = Vector{Float64}()
	push!(l,ω[1] - hγ[1][2] + ϵ)
	push!(u,ω[1] - hγ[1][1] - ϵ)

	F = Vector{Function}()
	for i in 1:2*m
		push!(F,(φ -> 0.))
	end
	
	F[1] = retro_function(1,ω,H,et)

	i = 1
	
	# At each step,...
	while keepon && i < n-1
		i += 1

		#... define the flow over edge i as a function of the synchronous frequency `φ`,...
		F[i] = retro_function(i,ω,H,et)
		
		#...determine if a solution exists in the tolerated interval of values for the synchronous frequency `φ`. If it does, adapt the the bounds on `φ` so the flow functions are all well-defined.
		Fl = F[i](l[i-1])
		Fu = F[i](u[i-1])

		if Fl < hγ[i][1] || Fu > hγ[i][2]
			@info "NO SOLUTION"
			keepon = false
		else
			if hγ[i][1] <= Fl <= hγ[i][2]
				le = l[end]
			else
				le = dichot((φ -> F[i](φ) - hγ[i][2]),l[end],u[end])
			end

			if hγ[i][1] <= Fu <= hγ[i][2]
				ue = u[end]
			else
				ue = dichot((φ -> F[i](φ) - hγ[i][1]),l[end],u[end])
			end
			push!(l,le + ϵ)
			push!(u,ue - ϵ)
		end
	end

	# Define the reversed flow functions for each edge.
	for i in 1:m
		F[i+m] = (φ -> H[i+m](F[i](φ)))
	end

	# If the solution exists, compute its sychronous frequency `φ`.
	if keepon
		jd = setdiff((B[n,1:m] .< -1e-8).*(1:m),[0,])
		if length(jd) == 0
			Fn = (φ -> ω[n] - φ)
		else
			Fn = (φ -> ω[n] - sum(F[j+m](φ) for j in jd) - φ)
		end
	
		Fl = Fn(l[end])
		Fu = Fn(u[end])
	
		if min(Fl,Fu) <= 0. <= max(Fl,Fu)
			φs = dichot(Fn,l[end],u[end])
		else
			@info "NO SOLUTION."
			keepon = false
		end
	end

	if keepon
		return true, [F[i](φs) for i in 1:2*m], φs
	else
		return false, [], Inf
	end
end


