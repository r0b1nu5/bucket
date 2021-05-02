include("tools.jl")

function acyc_algo(L::Array{Float64,2}, H::Function, ω0::Array{Float64,1}, hγ::Tuple{Float64,Float64}, α::Float64=.1)
	B0,w = L2B(L)
	n,m = size(B0)

	if m != n-1
		@info "The graph is not acyclic or not connected!"
	end

	B,id = reorder_nodes(B0,Array(1:n))
	ω = ω0[id]

	k = targets(B)
	ki = Dict{Int64,Array{Int64,1}}(i => [] for i in 1:n)
	for j in keys(k)
		push!(ki[k[j]],j)
	end

	l = Array{Float64,1}()
	u = Array{Float64,1}()
	push!(l,ω[1] - hγ[2])
	push!(u,ω[1] - hγ[1])

	F = Array{Function,1}()
	for i in 1:2*m
		push!(F,(φ -> 0.))
	end
	
	F[1] = retro_function(1,ω,H,α,k,ki)
#	F[1+m] = (φ -> H(retro_function(1,ω,H,α,k,ki)(φ),α))
#	F[1+m] = (φ -> H(F[1](φ),α))

	keepon = true
	i = 1
	
	while keepon && i < n-1
		i += 1

		F[i] = retro_function(i,ω,H,α,k,ki)
#		F[i+m] = (φ -> H(retro_function(i,ω,H,α,k,ki)(φ),α))
#		F[i+m] = (φ -> H(F[i](φ),α))
		
		Fl = F[i](l[i-1])
		Fu = F[i](u[i-1])

		if Fl < hγ[1] || Fu > hγ[2]
			@info "NO SOLUTION"
			keepon = false
		else
			if hγ[1] <= Fl <= hγ[2]
				le = l[end]
			else
				le = dichot((φ -> F[i](φ) - hγ[2]),l[end],u[end])
			end

			if hγ[1] <= Fu <= hγ[2]
				ue = u[end]
			else
				ue = dichot((φ -> F[i](φ) - hγ[1]),l[end],u[end])
			end
			push!(l,le)
			push!(u,ue)
		end
	end

	for i in 1:m
		F[i+m] = (φ -> H(F[i](φ),α))
	end

	if keepon
		jd = setdiff((B[n,1:m-1] .< -1e-8).*(1:m-1),[0,])
		if length(jd) == 0
			Fn = (φ -> φ)
		else
			Fn = (φ -> sum(F[j+m](φ) for j in jd) + φ)
		end
	
		Fl = Fn(l[end])
		Fu = Fn(u[end])
	
@info "l = $(l[end]), u = $(u[end])"
@info "Fl = $Fl, Fu = $Fu"

		if Fl <= ω[n] <= Fu
			φs = dichot(Fn,l[end],u[end])
		else
			@info "NO SOLUTION."
			keepon = false
		end
	end

	if keepon
		return [F[i](φs) for i in 1:2*m], φs
	else
		return nothing
	end
end

# Re-indexes the nodes and edges such that it satisfies the assumptions of the algorithm in the KS notes...

function reorder_nodes(B0::Array{Float64,2}, id0::Array{Int64,1})
	n,m = size(B0)

	if n == 2
		B2 = [1.,-1.]
		id2 = [id0[1],id0[2]]
	else
		Babs = abs.(B0)
	
		d = vec(sum(Babs,dims=2))
		mi,id = findmin(d)
		ma,ed = findmax(Babs[id,:])
	
		B1 = B0[[id;(1:id-1);(id+1:n)],[ed;(1:ed-1);(ed+1:n-1)]]
		id1 = [id0[id];id0[1:id-1];id0[id+1:n]]
		j1 = id1[findmax(abs.(B1[2:n,1]))[2] + 1]

		b1,i1 = reorder_nodes(B1[2:n,2:n-1],id1[2:n])

		j2 = findmin(abs.(i1 .- j1))[2]
		e2 = zeros(n-1)
		e2[j2] = -1.

		B2 = [abs.(B1[[1,],:]) ; e2 b1]
		id2 = [id0[id];i1]
	end

	return B2, id2
end

# Computes the target node of e_i according to the numbering in the KS notes...
# The incidence matrix B i assumed to be the output of "reorder_nodes"!!!

function targets(B::Array{Float64,2})
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

function retro_function(i::Int64, ω::Array{Float64,1}, H::Function, α::Float64, k::Dict{Int64,Int64}, ki::Dict{Int64,Array{Int64,1}})
	if length(ki[i]) == 0
		return (φ -> ω[i] - φ)
	else
		return (φ -> ω[i] - sum(H(retro_function(j,ω,H,α,k,ki)(φ),α) for j in ki[i]) - φ)
	end
end


