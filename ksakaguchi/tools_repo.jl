using LinearAlgebra, SparseArrays, Statistics

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


function order_nodes(B0::Union{Matrix{Float64},SparseMatrixCSC{Float64,Int64}}, id0::Vector{Int64}=Int64[], ed0::Vector{Int64}=Int64[])
	n,m = size(B0)

	if length(id0) == 0
		id0 = Array(1:n)
		ed0 = Array(1:m)
	end

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
		
		b1,i1,e1 = order_nodes(B1[2:n,2:n-1],id1[2:n],ed1[2:n-1])

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


