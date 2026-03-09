using LinearAlgebra

include("kmeans.jl")

function coarse_grain(A::Union{Matrix{Float64},SparseMatrixCSC{Float64,Int64}}, k::Int64, X::Matrix{Float64}, Y::Matrix{Float64})
	d = vec(sum(A,dims=1))
	
	N,T = size(X)
	n = floor(Int64,N/k)
	Atemp = A - diagm(0 => ones(N))

	r = collect(1:N)

	g = Vector{Vector{Int64}}()

	while length(r) > 1
		i = findmin(d)[2]
		j = findmax(Atemp[r[i],r])[2]
		push!(g,[r[i],r[j]])
		popat!(d,max(i,j))
		popat!(d,min(i,j))
		popat!(r,max(i,j))
		popat!(r,min(i,j))
	end
	if length(r) == 1
		push!(g,[r[1],])
	end

	AA = zeros(length(g),length(g))
	for i in 1:length(g)-1
		for j in i+1:length(g)
			if maximum(A[g[i],g[j]]) > 1e-2
				AA[i,j] = 1.
				AA[j,i] = 1.
			end
		end
	end

	X2 = zeros(0,T)
	Y2 = zeros(0,T)
	for k in g
		X2 = [X2;mean(X[k,:],dims=1)]
		Y2 = [Y2;mean(Y[k,:],dims=1)]
	end

	return g,X2,Y2,AA
end




function coarse_grain_1(A::Union{Matrix{Float64},SparseMatrixCSC{Float64,Int64}}, k::Int64, X::Matrix{Float64}, Y::Matrix{Float64})
        d = vec(sum(A,dims=1))
        L = diagm(0 => d) - A

        n,T = size(X)

        λs = eigvals(L)
        nz = sum(abs.(λs) .< 1e-10)

        us = eigvecs(L)
        u2 = us[:,nz+1]
        u3 = us[:,nz+2]

        c0 = 2*rand(k,2) .- 1
        g,c = my_kmeans([u2 u3],c0)

        X2 = zeros(0,T)
        Y2 = zeros(0,T)
        for i in 1:k
            X2 = [X2; mean(X[g .== i,:],dims=1)]
            Y2 = [Y2; mean(Y[g .== i,:],dims=1)]
        end
        
        return g,X2,Y2
end
