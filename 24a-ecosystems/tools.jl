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

function gen_A_correl(n::Int64, γ::Float64)
	sγ = sign(γ)
	aγ = abs(γ)

	A = zeros(n,n)
	if aγ < 1
		β = sqrt(aγ/(1-aγ))*sγ
		for i in 1:n
			for j in i+1:n
				x = 2*sqrt(3)*(rand() - .5)
				y = 2*sqrt(3)*(rand() - .5)
				z = 2*sqrt(3)*(rand() - .5)
				A[i,j] = (x + β*z)/sqrt(1+β^2)
				A[j,i] = (x + abs(β)*z)/sqrt(1+β^2)
			end
			A[i,i] = 2*sqrt(3)*(rand() - .5)
		end
	end

	return A
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


# Matches each eigenvalue of λ2 to the closest eigenvalue of λ1
function associate_eigvals(λ1::Vector{Complex{Float64}}, λ2::Union{Vector{Float64},Vector{Complex{Float64}}})
	λ3 = Complex{Float64}[]
	λt = copy(λ2)
	for λ in λ1
		l,i = findmin(abs.(λ2 .- λ))
		push!(λ3,λ2[i])
		popat!(λ2,i)
	end
	return λ3
end


# Plots a curve with increasing color density
function plot_density(x::Vector{Float64}, y::Vector{Float64}, col::Any, pα::Int64=1)
	n = length(x)
	αs = LinRange(0,1,n)

	for i in 4:5:n-3
		PyPlot.plot(x[i-3:i+3],y[i-3:i+3],color=col,alpha=αs[i]^pα)
	end

	return nothing
end
