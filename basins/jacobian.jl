using LinearAlgebra, SparseArrays

function jacobian(L::SparseMatrixCSC{Float64,Int64}, θ::Vector{Float64})
	A = spdiagm(0 => diag(L)) - L
	I,J,V = findnz(A)

	C = Vector{Float64}()

	for k in 1:length(I)
		i = I[k]
		j = J[k]
		push!(C,cos(θ[i]-θ[j]))
	end

	preJ = sparse(I,J,C)
	Jsp = preJ - spdiagm(0 => vec(sum(preJ,dims=1)))
	J = Symmetric(Matrix(Jsp))

	return J
end




