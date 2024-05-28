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


function jac_lv_bunin(N::Vector{Float64}, A::Matrix{Float64}, κ::Vector{Float64}, μsS::Float64=5., σsS::Float64=2.7)
	J = diagm(0 => N)*(-σsS*A .- μsS)
	J += diagm(0 => (κ - 2*N - σsS*A*N) .- μsS*sum(N))

	return J
end

function jac_lv_bunin(N::Vector{Float64}, A::Matrix{Float64}, κ::Float64=1., μsS::Float64=5., σsS::Float64=2.7)
	return jac_lv_bunin(N,A,κ*ones(length(N)),μsS,σsS)
end


function get_survivors(N::Vector{Float64}, zer0::Float64=1e-15)
	return vec(1:length(N))[N .> zer0]
end


function analyze_jac(N::Vector{Float64}, A::Matrix{Float64}, κ::Vector{Float64}, μsS::Float64=5., σsS::Float64=2.7)
	J = jac_lv_bunin(N,A,κ,μsS,σsS)
	eig = eigen(J)
	
	return J,eig.values,eig.vectors
end

function analyze_jac(N::Vector{Float64}, A::Matrix{Float64}, κ::Float64=1., μsS::Float64=5., σsS::Float64=2.7)
	return analyze_jac(N,A,κ*ones(length(N)),μsS,σsS)
end


function analyze_limit_cycle(Ns::Matrix{Float64}, A::Matrix{Float64}, κ::Vector{Float64}, μsS::Float64=5., σsS::Float64=2.7)
	Nm = vec(mean(Ns,dims=2))
	return analyze_jac(Nm,A,κ,μsS,σsS)
end

function analyze_limit_cycle(Ns::Matrix{Float64}, A::Matrix{Float64}, κ::Float64=1., μsS::Float64=5., σsS::Float64=2.7)
	return analyze_limit_cycle(Ns,A,κ*ones(size(A)[1]),μsS,σsS)
end





