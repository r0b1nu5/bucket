using LinearAlgebra, Graphs

function gen_wER(n::Int64, p::Float64)
	λ2 = -1000.
	A = zeros(n,n)

	while λ2 < 1e-10
		A = Matrix(adjacency_matrix(erdos_renyi(n,p)))
		D = diagm(0 => vec(sum(A,dims=1)))
		L = D - A
		λ2 = eigvals(L)[2]
	end

	w = rand(n,n)
	W = w - w'

	M = W.*A

	return M
end



