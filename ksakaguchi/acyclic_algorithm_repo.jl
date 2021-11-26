include("tools_repo.jl")

function run_acyclic_algorithm(L::Union{Matrix{Float64},SparseMatrixCSC{Float64,Int64}}, H0::Function, ω0::Vector{Float64}, hγ0::Tuple{Float64,Float64})
	B0,w = L2B(L)
	n,m = size(B0)
	
	H = Vector{Function}()
	hγ = Vector{Tuple{Float64,Float64}}()
	for i in 1:2*m
		push!(H,H0)
		push!(hγ,hγ0)
	end
	
	return run_acyclic_algorithm(L,H,ω0,hγ)
end

function run_acyclic_algorithm(L::Union{Matrix{Float64},SparseMatrixCSC{Float64,Int64}}, H0::Vector{Function}, ω0::Vector{Float64}, hγ0::Vector{Tuple{Float64,Float64}})
	B0,w = L2B(L)
	n,m = size(B0)
	
# Dimension checks
	if m != n-1
		@info "The graph is not acyclic of not connected: aborting."
		return nothing
	end
	if length(ω0) != n
		@info "Laplacian and natural frequencies do not have same dimension ($n,$(length(ω))): aborting."
		return nothing
	end
	if length(H0) != 2*m
		@info "There not the same number of coupling functions and edges ($(length(H0)),$(2*m)): aborting."
		return nothing
	end

# Sort nodes
	B,id,ed = order_nodes(B0)
	ω = ω0[id]
	H = H0[[ed;ed.+m]]
	hγ = hγ0[[ed;ed.+m]]


# Run the algorithm

	exist,f,φ,l,u = acyclic_algorithm(B,H,ω,hγ)


# De-sort nodes
	if exist
		di,de = desort_nodes(id,ed)
		Bf = B[di,de]
		t = [(minimum(B0[:,i].*Bf[:,]) < -.1) for i in 1:n-1]

		ff = [f[de + m*t];f[de +m*(1 .- t)]]
		lf = l[de]
		uf = u[de]
	else
		ff = f
		lf = l
		uf = u
	end

	return exist, ff, φ
end


function acyclic_algorithm(B::Union{Matrix{Float64},SparseMatrixCSC{Float64,Int64}}, H::Vector{Function}, ω::Vector{Float64}, hγ::Vector{Tuple{Float64,Float64}})
	n,m = size(B)

	keepon = true

	k = targets(B)
	ki = Dict{Int64,Array{Int64,1}}(i => [] for i in 1:n)
	for j in keys(k)
		push!(ki[k[j]],j)
	end

	l = Array{Float64,1}()
	u = Array{Float64,1}()
	push!(l,ω[1] - hγ[1][2])
	push!(u,ω[1] - hγ[1][1])

	F = Array{Function,1}()
	for i in 1:2*m
		push!(F,(φ -> 0.))
	end
	
	F[1] = retro_function(1,ω,H[1],k,ki)

	i = 1
	
	while keepon && i < n-1
		i += 1
		@info "$i"

		F[i] = retro_function(i,ω,H[i],k,ki)
		
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
			push!(l,le)
			push!(u,ue)
		end
	end

	for i in 1:m
		F[i+m] = (φ -> H[i+m](F[i](φ),α))
	end

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
		return true, [F[i](φs) for i in 1:2*m], φs, l, u
	else
		return false, [], Inf, Inf, -Inf
	end
end


function L2B(L::Array{Float64,2})
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
	
	return B,w
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

# Re-indexes the nodes and edges such that it satisfies the assumptions of the algorithm in the KS notes...

function reorder_nodes(B0::Union{Matrix{Float64},SparseMatrixCSC{Float64,Int64}})
	n,m = size(B0)
	id0 = Array(1:n)
	ed0 = Array(1:m)

	if n == 2
		if issparse(B0)
			B2 = sparse([1,2],[1,1],[1.,-1.])
		else
			B2 = [1.,-1.]
		end
		id2 = [id0[1],id0[2]]
		ed2 = [ed0[1],]
	else
		Babs = abs.(B0)
	
		d = vec(sum(Babs,dims=2))
		mi,id = findmin(d)
		ma,ed = findmax(Babs[id,:])
	
		B1 = B0[[id;(1:id-1);(id+1:n)],[ed;(1:ed-1);(ed+1:n-1)]]
		id1 = [id0[id];id0[1:id-1];id0[id+1:n]]
		j1 = id1[findmax(abs.(B1[2:n,1]))[2] + 1]
		ed1 = [ed0[ed];ed0[1:ed-1];ed0[ed+1:n-1]]

		b1,i1,e1 = reorder_nodes(B1[2:n,2:n-1],id1[2:n],ed1[2:n-1])

		j2 = findmin(abs.(i1 .- j1))[2]
		e2 = zeros(n-1)
		e2[j2] = -1.

		B2 = [abs.(B1[[1,],:]) ; e2 b1]
		id2 = [id0[id];i1]
		ed2 = [ed0[ed];e1]
	end

	return B2, id2, ed2
end

function reorder_nodes(B0::Union{Matrix{Float64},SparseMatrixCSC{Float64,Int64}}, id0::Vector{Int64}, ed0::Vector{Int64})
	n,m = size(B0)

	if n == 2
		B2 = [1.,-1.]
		id2 = [id0[1],id0[2]]
		ed2 = [ed0[1],]
	else
		Babs = abs.(B0)
	
		d = vec(sum(Babs,dims=2))
		mi,id = findmin(d)
		ma,ed = findmax(Babs[id,:])
	
		B1 = B0[[id;(1:id-1);(id+1:n)],[ed;(1:ed-1);(ed+1:n-1)]]
		id1 = [id0[id];id0[1:id-1];id0[id+1:n]]
		j1 = id1[findmax(abs.(B1[2:n,1]))[2] + 1]
		ed1 = [ed0[ed];ed0[1:ed-1];ed0[ed+1:n-1]]

		b1,i1,e1 = reorder_nodes(B1[2:n,2:n-1],id1[2:n],ed1[2:n-1])

		j2 = findmin(abs.(i1 .- j1))[2]
		e2 = zeros(n-1)
		e2[j2] = -1.

		B2 = [abs.(B1[[1,],:]) ; e2 b1]
		id2 = [id0[id];i1]
		ed2 = [ed0[ed];e1]
	end

	return B2, id2, ed2
end

function desort_nodes(id::Vector{Int64}, ed::Vector{Int64})
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

# Computes the target node of e_i according to the numbering in the KS notes...
# The incidence matrix B i assumed to be the output of "reorder_nodes"!!!

function targets(B::Union{Matrix{Float64},SparseMatrixCSC{Float64,Int64}})
	n,m = size(B)

	if n == 2
		k = Dict{Int64,Int64}(1 => 2)
	else
		k0 = targets(B[2:n,2:n-1])
		k = Dict{Int64,Int64}(i => k0[i-1]+1 for i in 2:n-1)
		k1 = findmax(abs.([0.;B[2:n,1]]))[2]
		k[1] = k1
	end

	return k
end


# Computes the zero of the STRICTLY MONOTONE function F in the interval [l,u], if it exists.

function dichot(F::Function, l0::Float64, u0::Float64, tol::Float64=1e-6)
	l = l0
	u = u0
	Fl = F(l)
	Fu = F(u)

	if Fl*Fu > 0.
		@info "NO SOLUTION."
		return nothing
	end

	s = sign(Fu-Fl)
	δ = (u - l)/10

	while δ > tol
		φ = LinRange(l,u,11)
		Fφ = [F(φ[i]) for i in 1:11]
		
		G = Fφ[1:10].*Fφ[2:11]
		mi,i = findmin(G)
		
		l = φ[max(i-1,1)]
		u = φ[min(i+1,11)]
		δ = u - l
	end

	return (l + u)/2
end

function retro_function(i::Int64, ω::Array{Float64,1}, H::Function, k::Dict{Int64,Int64}, ki::Dict{Int64,Array{Int64,1}})
	if length(ki[i]) == 0
		return (φ -> ω[i] - φ)
	else
		return (φ -> ω[i] - sum(H(retro_function(j,ω,H,k,ki)(φ)) for j in ki[i]) - φ)
	end
end


