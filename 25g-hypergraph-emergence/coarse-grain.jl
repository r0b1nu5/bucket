using LinearAlgebra

include("kmeans.jl")

function coarse_grain(A::Union{Matrix{Float64},SparseMatrixCSC{Float64,Int64}}, k::Int64, X::Matrix{Float64}, Y::Matrix{Float64})
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
