using LinearAlgebra, Graphs, Statistics 

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


function jac_lv_bunin(N::Vector{Float64}, A::Matrix{Float64}, κ::Vector{Float64}, μ::Float64=5., σ::Float64=2.7)
	S = length(N)
	J = diagm(0 => N)*(-σ/sqrt(S)*A .- μ/S)
	J += diagm(0 => (κ - 2*N - σ/sqrt(S)*A*N) .- μ/S*sum(N))

	return J
end

function jac_lv_bunin(N::Vector{Float64}, A::Matrix{Float64}, κ::Float64=1., μ::Float64=5., σ::Float64=2.7)
	return jac_lv_bunin(N,A,κ*ones(length(N)),μ,σ)
end


function get_survivors(N::Vector{Float64}, zer0::Float64=1e-15)
	return vec(1:length(N))[N .> zer0]
end


function analyze_limit_cycle(Ns::Matrix{Float64}, A::Matrix{Float64})
	Nm = vec(mean(Ns,dims=2))
	J = jac_lv_bunin(Nm,A)
	λs = eigvals(J)
	us = eigvecs(J)

	return J,λs,us
end




