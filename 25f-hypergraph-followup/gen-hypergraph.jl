using LinearAlgebra, Random

# Returns the adjacency tensor
# Only 3-edges
# n nodes
# p is the probability that an edge exists
# The hypergraph is connected

function rand_3_graph(n::Int64, p::Float64)
	m_max = binomial(n,3)
	m = sum(rand(m_max) .< p)
	
	A = zeros(0,4)
	E = Vector{Int64}[]
	B = zeros(n,0)
	connected = false
	c0 = 0
	while !connected && c0 < 1000
		c0 += 1

		A = zeros(0,4)
		a = zeros(n,n)
		E = Vector{Int64}[]
		B = zeros(n,0)
		for l in 1:m
			le = 0
			c1 = 0
			i,j,k = 0,0,0
			ijk = 0
			te = true
			while te && c1 < 1000
				c1 += 1
				ijk = sort(randperm(n)[1:3])
				if !(ijk in E)
					te = false
				end
			end
			if c1 == 1000
				@info "Not enough edges."
			else
				push!(E,ijk)
				i,j,k = ijk
				A = [A;[i j k 1.]]
				a[i,j]=a[i,k]=a[j,i]=a[j,k]=a[k,i]=a[k,j] = 1.
				B = [B can_bas(i,n) can_bas(j,n) can_bas(k,n)]
			end
		end
		l = adj2lap(a)
		λ2 = eigvals(l)[2]
		if abs(λ2) > 1e-4
			connected = true
		end
	end

	if c0 == 1000
		@info "WARNING: Hypergraph is not connected!"
	end

	return A,B,E
end

function can_bas(i::Int64, n::Int64)
	v = zeros(n)
	v[i] = 1
	return v
end

function adj2lap(A::Matrix{Float64})
	return diagm(0 => sum(A,dims=2)[:,1]) - A
end


