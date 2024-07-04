using DataDrivenDiffEq, ModelingToolkit, LinearAlgebra, DataDrivenSparse, LinearAlgebra, PyPlot, Combinatorics, Statistics

function this(X::Matrix{Float64}, Y::Matrix{Float64}, ooi::Vector{Int64}, dmax::Int64, λ::Float64=.1, ρ::Float64=1., niter::Int64=10)
	if size(X) != size(Y)
		@info "Time series' sizes do not match."
	end

	@info "THIS: getting the θs."
	θ,idx_mon = get_θ(X,dmax)
	@info "THIS: got the θs."

	coeff,err,relerr = mySINDy(θ,Y,λ,ρ,niter)
	
	return coeff, idx_mon, err, relerr
end
	


# Adapted from [https://github.com/eurika-kaiser/SINDY-MPC/blob/master/utils/sparsifyDynamics.m]
function mySINDy(θ::Matrix{Float64}, Y::Matrix{Float64}, λ::Float64=.1, ρ::Float64=1., niter::Int64=10)
	n,T = size(Y)
	m,T = size(θ)
	energy = sum(Y.^2)

	@info "SINDy: (pseudo)inverting the matrix."
#	Ξ = Y*θ'*pinv(θ*θ' + ρ*Id(m)) # Least square with Tikhonov regularization
	Ξ = Y*θ'*inv(θ*θ' + ρ*Id(m)) # Least square with Tikhonov regularization
	@info "SINDy: matrix (pseudo)inverted."

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
#			Ξ[i,biginds] = Y[[i,],:]*θ[biginds,:]'*pinv(θ[biginds,:]*θ[biginds,:]' + ρ*Id(sum(biginds)))
			Ξ[i,biginds] = Y[[i,],:]*θ[biginds,:]'*inv(θ[biginds,:]*θ[biginds,:]' + ρ*Id(sum(biginds)))
		end
	end

#	cost = sum((Y - Ξ*θ).^2) + ρ*sum(Ξ.^2)
	err = sum((Y - Ξ*θ).^2)

	return Ξ, err, err/energy
end

function Id(n::Int64)
	return diagm(0 => ones(n))
end

function get_θ_old(X::Matrix{Float64}, dmax::Int64)
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

function get_θ(X::Matrix{Float64}, dmax::Int64)
	n,T = size(X)

	θ = zeros(0,T)
	idx_mon = Dict{Int64,Vector{Int64}}()

	if dmax < 0
		@info "dmax < 0..."
	elseif dmax == 0
		θ = ones(1,T)
		idx_mon = Dict{Int64,Vector{Int64}}(1 => Int64[])
	else
		θ = [ones(1,T);X]
		idx_mon = Dict{Int64,Vector{Int64}}(i+1 => [i,] for i in 1:n)
		idx_mon[1] = Int64[]

		c = n+1
		for d in 2:dmax
			for comb in combinations(1:n,d)
				θ = [θ;prod(X[comb,:],dims=1)]
				c += 1
				idx_mon[c] = comb
			end
		end
	end

	return θ, idx_mon
end


