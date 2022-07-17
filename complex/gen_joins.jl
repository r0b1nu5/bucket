using LinearAlgebra

function gen_joins(A::Matrix{Float64}, As::Vector{Matrix{Float64}})
	ns = [size(As[i])[1] for i in 1:length(As)]
	n = length(ns)
	N = sum(ns)

	AA = zeros(N,N)

	for i in 1:n
		m = sum(ns[1:i-1])
		AA[m+1:m+ns[i],m+1:m+ns[i]] = As[i]
		for j in i+1:n
			mi = sum(ns[1:i-1])
			mj = sum(ns[1:j-1])
			AA[mi+1:mi+ns[i],mj+1:mj+ns[j]] = A[i,j]*ones(ns[i],ns[j])
			AA[mj+1:mj+ns[j],mi+1:mi+ns[i]] = A[j,i]*ones(ns[j],ns[i])
		end
	end

	return AA, ns
end


function broadcast(θ::Vector{Float64}, ns::Vector{Int64})
	θθ = Float64[]
	for i in 1:length(ns)
		θθ = [θθ;θ[i]*ones(ns[i])]
	end

	return θθ
end

function A2L(A::Matrix{Float64})
	return diagm(0 => vec(sum(A,dims=1))) - A
end

function er_ntw(n::Int64, p::Float64)
	A = zeros(n,n)

	for i in 1:n-1
		for j in i+1:n
			if rand() < p
				A[i,j] = 1.
				A[j,i] = 1.
			end
		end
	end

	return A
end




