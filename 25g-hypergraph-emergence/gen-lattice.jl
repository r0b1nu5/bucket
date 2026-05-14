using LinearAlgebra, Random, SparseArrays, PyPlot

function gen_square_lattice(n::Int64, m::Int64)
	B1 = sparse([1:n-1;2:n],[1:n-1;1:n-1],[ones(n-1);-ones(n-1)])
	I = spdiagm(0 => ones(n))

	B = spzeros(n*m,2*(n-1)*(m-1)+(n-1)+(m-1))
	for i in 1:m
		B[(1:n) .+ (i-1)*n,(1:n-1) .+ (i-1)*(n-1)] = B1
	end
	for i in 1:m-1
		B[(1:n) .+ (i-1)*n,(1:n) .+ (n-1)*m .+ n*(i-1)] = I
		B[(1:n) .+ i*n,(1:n) .+ (n-1)*m .+ n*(i-1)] = -I
	end

	L = B*B'
	A = spdiagm(diag(L)) - L
	
	return A,B
end

function gen_square_lattice(n::Int64)
	return gen_square_lattice(n,n)
end




