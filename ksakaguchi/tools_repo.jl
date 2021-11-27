using LinearAlgebra, SparseArrays, Statistics

"""
	L2B(L::Union{Matrix{Float64},SparseMatrixCSC{Float64,Int64}})

Compute an incidence matrix associated to a Laplacian matrix.

_INPUT_:\\
`L`: Laplacian matrix (dense or sparse).

_OUTPUT_:\\
`B`: Incidence matrix of the graph (same sparsity as the input).\\
`w`: Vector of edge weights.\\
`Bt`: Transpose of B.
"""
function L2B(L::Matrix{Float64})
	n = size(L)[1]
	
	B = Array{Float64,2}(undef,n,0)
	w = Array{Float64,1}(undef,0)
	for i in 1:n-1
		for j in i+1:n
			if L[i,j] != 0.0
				ed = zeros(n)
				ed[i] = 1.0
				ed[j] = -1.0
				B = [B ed]
				push!(w,-L[i,j])
			end
		end
	end
	
	return B,w,Array(B')
end

function L2B(L::SparseMatrixCSC{Float64,Int})
	n = size(L)[1]
	
	B = Array{Float64,2}(undef,n,0)
	w = Array{Float64,1}(undef,0)
	Bt = Array{Float64,2}(undef,0,n)
	for i in 1:n-1
		for j in i+1:n
			if L[i,j] != 0.0
				ed = zeros(n)
				edt = zeros(1,n)
				ed[i] = 1.0
				edt[1,i] = 1.0
				ed[j] = -1.0
				edt[1,j] = -1.0
				B = [B ed]
				Bt = [Bt;edt]
				push!(w,-L[i,j])
			end
		end
	end
	
	return sparse(B),w,sparse(Bt)
end


"""
	reindex(B::Union{Matrix{Float64},SparseMatrixCSC{Float64,Int64}}, id::Vector{Int64}=Int64[], ed::Vector{Int64}=Int64[])

Re-orders the nodes and edges of the graph to satisfy the conditions for Theorem XXX [Delabays, Jafarpour, and Bullo (2021)] to apply. 

_INPUT_:\\
`B`: Incidence matrix of the _bidirected_ graph. \\
`id`: Indices of the nodes. If none is given, the indexing is assumed to be sequential, i.e., id = 1:n.\\
`ed`: Indices of the edges. If none is given, the indexing is assumed to be sequential, i.e., ed = 1:n-1.

_OUTPUT_:\\
`B2`: Re-ordered incidence matrix of the _bidirected_ graph. The graph structure is preserved, only the numbering is modified. \\
`id2`: Vector of the original node indexing, with their new position in matrix B2.\\
`ed2`: Vector of the original edge indexing, with their new position in matrix B2.
"""
function reindex(B::Union{Matrix{Float64},SparseMatrixCSC{Float64,Int64}}, id::Vector{Int64}=Int64[], ed::Vector{Int64}=Int64[])
	n,m2 = size(B)
	m = Int64(m2/2)

	if length(id) == 0
		id = Array(1:n)
		ed = Array(1:2*m)
	end

	if n == 2
		# If there are two nodes, the ordering is straighforward.
		B2 = [1. -1.;-1. 1.]
		if B[1,1] > 0.
			id2 = [id[1],id[2]]
			ed2 = [ed[1],ed[2]]
		else
			id2 = [id[2],id[1]]
			ed2 = [ed[1],ed[2]]
		end
	else
		# If there are more than two nodes, then...
		Babs = abs.(B[:,1:m])
	
		# ... identify a node with a unique neighbor, ...
		d = vec(sum(Babs,dims=2))
		mi,id0 = findmin(d)
		ma,ed0 = findmax(Babs[id0,:])
	
		# ... re-index all nodes so it is the first one, ...
		if B[id0,ed0] < 0.
			ED0 = [ed0+m;(1:ed0-1);(ed0+1:m);ed0;(m+1:m+ed0-1);(m+ed0+1:2*m)]
		else
			ED0 = [ed0;(1:ed0-1);(ed0+1:m);m+ed0;(m+1:m+ed0-1);(m+ed0+1:2*m)]
		end
		B1 = B[[id0;(1:id0-1);(id0+1:n)],ED0]
		id1 = [id[id0];id[1:id0-1];id[id0+1:n]]
		j1 = id1[findmax(abs.(B1[2:n,1]))[2] + 1]
		ed1 = ed[ED0]
		
		# ... and apply the same process to the submatrix with first row and column removed.
		b1,i1,e1 = reindex(B1[2:n,[2:m;m+2:2*m]],id1[2:n],ed1[[2:m;m+2:2*m]])

		j2 = findmin(abs.(i1 .- j1))[2]
		e2 = zeros(n-1)
		e2[j2] = -1.

		B2 = [[1.;e2] [zeros(1,m-1);b1[:,1:m-1]] [-1.;-e2] [zeros(1,m-1);b1[:,m:2*m-2]]]
		id2 = [id[id0];i1]
		ed2 = [ed1[1];e1[1:m-1];ed1[m+1];e1[m:end]]
	end

	return B2, id2, ed2
end

"""
	deindex(id::Vector{Int64}, ed::Vector{Int64})

Reverses the re-numbering operated by `reindex`.

_INPUT_:\\
`id`: Vector of the original node indexing, with their new position in the output matrix of `reindex`. This is the second output of `reindex`.\\
`ed`: Vector of the original edge indexing, with their new position in the output matrix of `reindex`. This is the third output of `reindex`.

_OUTPUT_:\\
`di`: Inverse indexing of nodes to recover the original ordering of nodes, before applying `reindex`.\\
`de`: Inverse indexing of edges to recover the original ordering of edges, berore applying `reindex`.
"""
function deindex(id::Vector{Int64}, ed::Vector{Int64})
	n = length(id)
	m = length(ed)

	di = Vector{Int64}()
	for i in 1:n
		push!(di,findmin(abs.(id .- i))[2])
	end
	de = Vector{Int64}()
	for i in 1:m
		push!(de,findmin(abs.(ed .- i))[2])
	end

	return di,de
end


"""
	targets(B::Union{Matrix{Float64},SparseMatrixCSC{Float64,Int64}})

Given an incidence matrix for an acyclic graph, with indexing satisfying the assumptions of Theorem XXX in [Delabays, Jafarpour, and Bullo (2021)], returns, for each node i:\\
- The unique node j>i to which it is connected;\\
- The list of nodes j<i that are connected to it.

_INPUT_:\\
`B`: Incidence matrix of the graph. The graph is assumed connected and acyclic, and its numbering is assumed to satisfy the assumptions of Theorem XXX in [Delabays, Jafarpour, and Bullo (2021)].

_OUTPUT_:\\
`te`: Dictionary associating to each node i = 1:n-1, the unique node with larger index to which it is connected.\\
`et`: Dictionary associating to each node i = 1:n, the list of nodes with smaller index to which it is connected.
"""
function targets(B::Union{Matrix{Float64},SparseMatrixCSC{Float64,Int64}})
	n,m = size(B)

	if n == 2
		te = Dict{Int64,Int64}(1 => 2)
		et = Dict{Int64,Vector{Int64}}(1 => Int64[], 2 => [1,])
	else
		k0,e0 = targets(B[2:n,2:n-1])
		te = Dict{Int64,Int64}(i => k0[i-1]+1 for i in 2:n-1)
		k1 = findmax(abs.([0.;B[2:n,1]]))[2]
		te[1] = k1
		et = Dict{Int64,Vector{Int64}}(i => e0[i-1] .+ 1 for i in 2:n)
		push!(et[k1],1)
		et[1] = Int64[]
	end

	return te,et
end


"""
 	dichot(F::Function, l::Float64, u::Float64, tol::Float64=1e-6)

Computes the zero of the function `F` in the interval [`l`,`u`] with tolerance `tol`, if it exits. The function `F` is assumed strictly monotone and well-defined on [`l`,`u`]. The algorithm proceeds by dichotomy.

_INPUT_:\\
`F`: Function whose zero needs to be found. It is assumed strictly monotone.\\
`l`: Lower bound of the interval where the zero needs to be found.\\
`u`: Upper bound of the interval where the zero needs to be found. \\
`tol`: Tolerance on the error to the actual solution.

_OUTPUT_:\\
`z`: Zero of the function `F` over [`l`,`u`], if it exists. 
"""
function dichot(F::Function, l::Float64, u::Float64, tol::Float64=1e-6)
	ll = l
	uu = u
	Fl = F(ll)
	Fu = F(uu)

	if Fl*Fu > 0.
		@info "NO SOLUTION."
		return nothing
	end

	s = sign(Fu-Fl)
	δ = (uu - ll)/10

	while δ > tol
		φ = LinRange(ll,uu,11)
		Fφ = [F(φ[i]) for i in 1:11]
		
		G = Fφ[1:10].*Fφ[2:11]
		mi,i = findmin(G)
		
		ll = φ[max(i-1,1)]
		uu = φ[min(i+1,11)]
		δ = uu - ll
	end

	return (ll + uu)/2
end


"""
	retro_function(i::Int64, ω::Vector{Float64}, H::Vector{Function}, et::Dict{Int64,Vector{Int64}}) 

Recursively defines the flow over an edge as a function of the synchronous frequency `φ`.

_INPUT_"\\
`i`: Index of the edge whose flow function has to be defined. Note that, according to indexing, `i` is also the index of the node at the source of edge `i`.\\
`ω`: Vector of natural frequencies of the oscillators. \\
`H`: Vector of the transfer functions over all (directed) edges. \\
`et`: Dictionary associating to each node i, the list of nodes with lower index, connected to it. This is the second output of `targets`.

_OUTPUT_:\\
`F`: Flow function over edge `i`, with respect to the synchronous frequency `φ`.
"""
function retro_function(i::Int64, ω::Vector{Float64}, H::Vector{Function}, et::Dict{Int64,Vector{Int64}})
	if length(et[i]) == 0
		# If node i has a unique neighbor, the relation between flow and synchronous frequency is straightforward.
		return (φ -> ω[i] - φ)
		# If node i has multiple neighbors, the flow depends on the sum of the neighboring flows, hence the sum over neighbors (except one).
	else
		return (φ -> ω[i] - sum(H[j+m](retro_function(j,ω,H,et)(φ)) for j in et[i]) - φ)
	end
end


