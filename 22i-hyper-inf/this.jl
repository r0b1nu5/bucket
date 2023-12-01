using DataDrivenDiffEq, ModelingToolkit, LinearAlgebra, DataDrivenSparse, LinearAlgebra, PyPlot, Combinatorics, Statistics

function this(X::Matrix{Float64}, Y::Matrix{Float64}, ooi::Vector{Int64}, dmax::Int64, λ::Float64=.1, ρ::Float64=1., niter::Int64=10)
	if size(X) != size(Y)
		@info "Time series' sizes do not match."
	end

	θ = get_θ(X,dmax)

	return mySINDy(θ,Y,λ,ρ,niter)
end
	


# Adapted from [https://github.com/eurika-kaiser/SINDY-MPC/blob/master/utils/sparsifyDynamics.m]
function mySINDy(θ::Matrix{Float64}, Y::Matrix{Float64}, λ::Float64=.1, ρ::Float64=1., niter::Int64=10)
	n,T = size(Y)
	m,T = size(θ)
	energy = sum(abs.(Y))

	Ξ = Y*θ'*pinv(θ*θ' + ρ*Id(m)) # Least square with Tikhonov regularization

	nz = 1e6
	k = 1

	while k < niter && sum(abs.(Ξ) .> 1e-6) != nz
		k += 1
		nz = sum(abs.(Ξ) .> 1e-6)
#		err = sum((Y - Ξ*θ).^2) + ρ*sum(Ξ.^2)
		err = sum((Y - Ξ*θ).^2)
		@info "iter $k: $(nz) nonzero coefficients, error = $(round(err)), rel-err = $(round(err/energy*100,digits=2))%"

		smallinds = (abs.(Ξ) .< λ)
		Ξ[smallinds] .= 0.
		for i in 1:n
			biginds = .~smallinds[i,:]
			Ξ[i,biginds] = Y[[i,],:]*θ[biginds,:]'*pinv(θ[biginds,:]*θ[biginds,:]' + ρ*Id(sum(biginds)))
		end
	end

	err = sum((Y - Ξ*θ).^2) + ρ*sum(Ξ.^2)

	return Ξ, err, err/energy
end

function Id(n::Int64)
	return diagm(0 => ones(n))
end

function get_θ(X::Matrix{Float64}, dmax::Int64)
	n,T = size(X)

	@variables x[1:n]
	prebasis = polynomial_basis([x[i] for i in 1:n],dmax)
	basis = Basis(prebasis,[x[i] for i in 1:n])

	θ = zeros(length(basis),0)
	for t in 1:T
		θ = [θ basis(X[:,t])]
	end

	return θ
end
