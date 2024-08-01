using DataDrivenDiffEq, ModelingToolkit, DataDrivenSparse, LinearAlgebra, PyPlot, Combinatorics, Statistics

function this(X::Matrix{Float64}, Y::Matrix{Float64}, ooi::Vector{Int64}, dmax::Int64, λ::Float64=.1, ρ::Float64=1., niter::Int64=10)
	if size(X) != size(Y)
		@info "Time series' sizes do not match."
	end

#	@info "THIS: getting the θs."
	θ,d = get_θ(X,dmax)
	idx_mon = Dict{Int64,Vector{Int64}}()
	for i in 1:size(d)[1]
		mon = d[i,:][d[i,:] .!= 0]
		if length(mon) == length(union(mon))
			idx_mon[i] = sort(mon)
		end
	end
#	@info "THIS: got the θs."

	coeff,err,relerr = mySINDy(θ,Y,λ,ρ,niter)
	
	return coeff, idx_mon, err, relerr
end

# Version of THIS with restriction on the monomials used.
# Forbidden pairs of variables are in 'forbid'.
function this(X::Matrix{Float64}, Y::Matrix{Float64}, ooi::Vector{Int64}, dmax::Int64, forbid::Vector{Vector{Int64}}, λ::Float64=.1, ρ::Float64=1., niter::Int64=10)
	if size(X) != size(Y)
		@info "Time series' sizes do not match."
	end

#	@info "THIS: getting the θs."
	θ,d = get_θ(X,dmax)
	ids = Int64[]
	for i in 1:size(θ)[1]
		if length(intersect(collect(combinations(sort(setdiff(d[i,:],[0,])),2)),forbid)) == 0
			push!(ids,i)
		end
	end
	θ = θ[ids,:]
	d = d[ids,:]
	idx_mon = Dict{Int64,Vector{Int64}}()
	for i in 1:size(d)[1]
		mon = d[i,:][d[i,:] .!= 0]
		if length(mon) == length(union(mon))
			idx_mon[i] = sort(mon)
		end
	end
#	@info "THIS: got the θs."

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
#	@info "SINDy: matrix (pseudo)inverted."

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
#			@info "i/n = $i/$n"
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

function get_θ_oldold(X::Matrix{Float64}, dmax::Int64)
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

function get_θ_old(X::Matrix{Float64}, dmax::Int64)
	@assert dmax > 0

	n,T = size(X)
	XX = [X[i,:] for i in 1:n]

	θ = zeros(0,T)
	idx_mon = Dict{Int64,Vector{Int64}}()

# Inspired by the 'polynomial_basis' generator in "~/.julia/packages/DataDrivenDiffEq/iP5FS/src/utils/basis_generators.jl".	
	n_x = n
	n_c = binomial(n_x + dmax,dmax) # Black magic!
	θ = Matrix{Float64}(undef,n_c,T)
	# Creating the iterators for powers
	_check_degree(x) = sum(x) <= dmax ? true : false
	itr = Base.Iterators.product([0:dmax for i in 1:n_x]...)
	itr_ = Base.Iterators.Stateful(Base.Iterators.filter(_check_degree,itr))
	filled = false
	for i in 1:n_c
		θ[i,:] = ones(1,T)
		filled = true
		vars = Int64[]
		id = 0
		for (Xi,ci) in zip(XX,popfirst!(itr_))
			id += 1
			vars = [vars;id*ones(Int64,ci)]
			if !iszero(ci)
				filled ? θ[i,:] = Xi.^ci : θ[i,:] .*= Xi.^ci
				filled = false
			end
		end
		if length(vars) == length(union(vars))
			idx_mon[i] = vars
		end
	end


	return θ, idx_mon
end


function get_θ(X::Matrix{Float64}, dmax::Int64, i0::Int64=1)
	n,T = size(X)

	θ = ones(1,T)
	d = zeros(Int64,1,dmax)

	if dmax == 0
		return θ,d
	else
		for i in 1:n
			θ0,d0 = get_θ(X[i:n,:],dmax-1,i)
			#d0 += (i-1)*(d0 .> 0)
			θ = [θ; repeat(X[[i,],:],size(θ0)[1],1).*θ0]
			d = [d; d0 i*ones(Int64,size(d0)[1],1)]
			#d = [d; (d0.+i.-1) i*ones(Int64,size(d0)[1],1)]
		end

		d += (i0-1)*(d .> 0)
		
		return θ,d
	end
end


