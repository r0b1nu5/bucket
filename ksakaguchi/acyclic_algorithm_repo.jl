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
		t = [(minimum(B0[:,i].*Bf[:,i]) < -.1) for i in 1:n-1]

		ff = [f[de + m*t];f[de + m*(1 .- t)]]
		lf = l[de]
		uf = u[de]
	else
		ff = f
		lf = l
		uf = u
	end

	return exist, ff, φ
end


function acyclic_algorithm(B::Union{Matrix{Float64},SparseMatrixCSC{Float64,Int64}}, H::Vector{Function}, ω::Vector{Float64}, hγ::Vector{Tuple{Float64,Float64}}, ϵ::Float64=1e-10)
	n,m = size(B)

	keepon = true

	te,et = targets(B)
#=	ki = Dict{Int64,Vector{Int64}}(i => [] for i in 1:n)
	for j in keys(k)
		push!(ki[k[j]],j)
	end
=#

	l = Vector{Float64}()
	u = Vector{Float64}()
	push!(l,ω[1] - hγ[1][2] + ϵ)
	push!(u,ω[1] - hγ[1][1] - ϵ)

	F = Vector{Function}()
	for i in 1:2*m
		push!(F,(φ -> 0.))
	end
	
	F[1] = retro_function(1,ω,H,te,et)

	i = 1
	
	while keepon && i < n-1
		i += 1
		@info "$i"

		F[i] = retro_function(i,ω,H,te,et)
		
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

	for i in 1:m
		F[i+m] = (φ -> H[i+m](F[i](φ)))
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


