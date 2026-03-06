using LinearAlgebra, SparseArrays

function kuramoto(θ0::Vector{Float64}, ω::Vector{Float64}, B::Union{Matrix{Float64},SparseMatrixCSC{Float64,Int64}}, K::Vector{Float64}, ϕ::Union{Float64,Vector{Float64}}, h::Float64=.01, maxiter::Int64=1000, tol::Float64=1e-6)
	n = length(θ0)

	θ = copy(θ0)
	θs = copy(θ0)
	dθ = 1000.
	dθs = zeros(n,0)
	
	iter = 0

	while iter < maxiter && maximum(abs.(dθ)) > tol
		iter += 1

		k1 = f_kuramoto(θ,ω,B,K,ϕ)
		k2 = f_kuramoto(θ + h/2*k1,ω,B,K,ϕ)
		k3 = f_kuramoto(θ + h/2*k2,ω,B,K,ϕ)
		k4 = f_kuramoto(θ + h*k3,ω,B,K,ϕ)

		dθ = (k1+2*k2+2*k3+k4)./6
		θ += h*dθ

		θs = [θs θ]
		dθs = [dθs dθ]
	end

	return θs,dθs
end

function f_kuramoto(θ::Vector{Float64}, ω::Vector{Float64}, B::Matrix{Float64}, K::Vector{Float64}, ϕ::Vector{Float64})
	return ω - B*diagm(0 => K)*(sin.(B'*θ - ϕ) - sin.(ϕ))
end

function f_kuramoto(θ::Vector{Float64}, ω::Vector{Float64}, B::SparseMatrixCSC{Float64,Int64}, K::Vector{Float64}, ϕ::Vector{Float64})
	return ω - B*spdiagm(K)*(sin.(B'*θ - ϕ) - sin.(ϕ))
end

function f_kuramoto(θ::Vector{Float64}, ω::Vector{Float64}, B::Union{Matrix{Float64},SparseMatrixCSC{Float64,Int64}}, K::Float64, ϕ::Float64)
	n = length(θ)
	
	return f_kuramoto(θ,ω,B,K*ones(n),ϕ*ones(n))
end

function f_kuramoto(θ::Matrix{Float64}, ω::Vector{Float64}, B::Union{Matrix{Float64},SparseMatrixCSC{Float64,Int64}}, K::Union{Float64,Vector{Float64}}, ϕ::Union{Float64,Vector{Float64}})
	n,T = size(θ)
	X = zeros(n,0)
	
	for t in 1:T
		X = [X f_kuramoto(θ[:,t],ω,B,K,ϕ)]
	end

	return X
end



