using LinearAlgebra

function Leps(x0::Array{Float64,1}, eps::Float64)
	n = length(x0)

	A = Float64.((0 .< abs.(repeat(x0,1,n) - repeat(x0',n,1)) .< eps))
	D = diagm(0 => vec(sum(A,dims=1)))
	L = D - A

	return L, A, D
end


