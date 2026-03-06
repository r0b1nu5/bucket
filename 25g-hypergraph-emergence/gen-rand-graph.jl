using LinearAlgebra, Random, SparseArrays

function gen_rand_graph(n::Int64, p::Float64)
	Bfull = gen_full_incidence(n)

	w = rand(Int64(n*(n-1)/2))
	B = Bfull[:,w .< p]

	L = B*B'
	A = spdiagm(diag(L)) - L

	return A,B
end

function gen_full_incidence(n::Int64)
	if n > 3
		B1 = sparse([ones(n-1);2:n],[1:n-1;1:n-1],[ones(n-1);-ones(n-1)])
		B = [B1 [spzeros(1,Int64((n-1)*(n-2)/2));gen_full_incidence(n-1)]]
	elseif n == 3
		B = sparse([1,2,1,3,2,3],[1,1,2,2,3,3],[1.,-1.,1.,-1.,1.,-1.])
	end

	return B
end

	



