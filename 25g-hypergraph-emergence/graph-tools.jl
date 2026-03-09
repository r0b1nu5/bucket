using LinearAlgebra, Random, SparseArrays, PyPlot

function gen_full_incidence(n::Int64)
	if n > 3
		B1 = sparse([ones(n-1);2:n],[1:n-1;1:n-1],[ones(n-1);-ones(n-1)])
		B = [B1 [spzeros(1,Int64((n-1)*(n-2)/2));gen_full_incidence(n-1)]]
	elseif n == 3
		B = sparse([1,2,1,3,2,3],[1,1,2,2,3,3],[1.,-1.,1.,-1.,1.,-1.])
	end

	return B
end

	
function gen_rand_graph(n::Int64, p::Float64)
	Bfull = gen_full_incidence(n)

	w = rand(Int64(n*(n-1)/2))
	B = Bfull[:,w .< p]

	L = B*B'
	A = spdiagm(diag(L)) - L

	return A,B
end


function plot_graph(A::Union{Matrix{Float64}, Matrix{Int64}, SparseMatrixCSC{Float64,Int64}})
    n = size(A)[1]

    d = vec(sum(A,dims=1))
    L = diagm(0 => d) - A

    λs = eigvals(L)
    nz = sum(abs.(λs) .< 1e-10)

    us = eigvecs(L)
    u2 = us[:,nz+1] + .1*rand(n)
    u3 = us[:,nz+2] + .1*rand(n)

    for i in 1:n-1
        for j in i+1:n
            if abs(A[i,j]) > 1e-10
                PyPlot.plot(u2[[i,j]],u3[[i,j]],"-k",lw=2)
            end
        end
    end
    PyPlot.plot(u2,u3,"ok")
end


function plot_hypergraph(A::Union{Matrix{Float64}, Matrix{Int64}, SparseMatrixCSC{Float64,Int64}}, Ainf::Dict{Int64,Matrix{Float64}})
    n = size(A)[1]

    d = vec(sum(A,dims=1))
    L = diagm(0 => d) - A

    λs = eigvals(L)
    nz = sum(abs.(λs) .< 1e-10)

    us = eigvecs(L)
    u2 = us[:,nz+1] + .1*rand(n)
    u3 = us[:,nz+2] + .1*rand(n)

    for k in setdiff(sort(collect(keys(Ainf))),[1,])
        AA = Ainf[k]
        for i in 1:size(AA)[1]
            js = Int64.(AA[i,1:end-1])
            PyPlot.fill(u2[js],u3[js],"-",color="C0",alpha=.5)
        end
    end

    PyPlot.plot(u2,u3,"ok")
end

