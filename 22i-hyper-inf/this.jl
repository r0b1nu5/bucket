using DataDrivenDiffEq, ModelingToolkit, DataDrivenSparse, LinearAlgebra, PyPlot, Combinatorics, Statistics, Distributed

function this(X::Matrix{Float64}, Y::Matrix{Float64}, ooi::Vector{Int64}, dmax::Int64, λ::Float64=.1, ρ::Float64=1., niter::Int64=10)
	if size(X) != size(Y)
		@info "Time series' sizes do not match."
	end

#	@info "THIS: getting the θs."
	θ,d = get_θd(X,dmax)
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
function this_filter(X::Matrix{Float64}, Y::Matrix{Float64}, ooi::Vector{Int64}, dmax::Int64, keep::Vector{Vector{Int64}}, λ::Float64=.1, ρ::Float64=1., niter::Int64=10)
	if size(X) != size(Y)
		@info "Time series' sizes do not match."
	end

	n,T = size(X)

	
#	@info "THIS: getting the θs."
	d = get_d(n,dmax)
	d2i = Dict{Vector{Int64},Int64}(sort(d[i,:]) => i for i in 1:size(d)[1])

	i2keep = Dict{Int64,Vector{Int64}}(i => Int64[] for i in 1:n)
	for k in 1:length(keep)
		i,j = keep[k]
		push!(i2keep[i],k)
		push!(i2keep[j],k)
	end


#=
	ids = Dict{Int64,Vector{Int64}}(i => Int64[] for i in 1:n) # for each agent index returns the indices of relevant monomials (according to 'keep')
	θs = Dict{Int64,Matrix{Float64}}(i => zeros(0,T) for i in 1:n) # for each agent index returns the values of the relevant monomials (according to 'keep')

		append!(ids[i],[d2i[v] for v in [sort([j,k]) for k in 0:n]])
		io = open("data/thetas_$i.csv","a")
		writedlm(io,reduce(vcat,[prod(X[v,:],dims=1) for v in [setdiff([j,k],[0,]) for k in 0:n]]))
		close(io)
#		θs[i] = vcat(θs[i],reduce(vcat,[prod(X[v,:],dims=1) for v in [setdiff([j,k],[0,]) for k in 0:n]]))
		append!(ids[j],[d2i[v] for v in [sort([i,k]) for k in 0:n]])
		io = open("data/thetas_$j.csv","a")
		writedlm(io,reduce(vcat,[prod(X[v,:],dims=1) for v in [setdiff([i,k],[0,]) for k in 0:n]]))
		close(io)
#		θs[j] = vcat(θs[j],reduce(vcat,[prod(X[v,:],dims=1) for v in [setdiff([i,k],[0,]) for k in 0:n]]))
	end

=#	


	coeff = Dict{Int64,Matrix{Float64}}()
	ids = Dict{Int64,Vector{Int64}}(i => Int64[] for i in 1:n) # for each agent index returns the indices of relevant monomials (according to 'keep')
	for i in 1:n
		θ = ones(1,T)
		id = [1,]
		for k in i2keep[i]
			j = setdiff(keep[k],[i,])[1]
			θ = vcat(θ,reduce(vcat,[prod(X[v,:],dims=1) for v in [setdiff([j,l],[0,]) for l in 0:n]]))
			append!(id,[d2i[v] for v in [sort([j,l]) for l in 0:n]])
		end

		iii = unique(a -> id[a], eachindex(id))
		id = id[iii]
		θ = θ[iii,:]

		@info "============ mySINDy for agent $i =============="
		if isempty(id)
			coef = zeros(0,0)
		else
			coef,err,relerr = mySINDy(θ,Y[[i,],:],λ,ρ,niter)
			ids[i] = id
		end

		coeff[i] = coef
	end

	return coeff, ids, 1., 1.
end	

# Parallel version
function this_par_filter(X::Matrix{Float64}, Y::Matrix{Float64}, ooi::Vector{Int64}, dmax::Int64, keep::Vector{Vector{Int64}}, λ::Float64=.1, ρ::Float64=1., niter::Int64=10)
	if size(X) != size(Y)
		@info "Time series' sizes do not match."
	end

	n,T = size(X)

	
	d = get_d(n,dmax)
	d2i = Dict{Vector{Int64},Int64}(sort(d[i,:]) => i for i in 1:size(d)[1])

	i2keep = Dict{Int64,Vector{Int64}}(i => Int64[] for i in 1:n)
	for k in 1:length(keep)
		i,j = keep[k]
		push!(i2keep[i],k)
		push!(i2keep[j],k)
	end

	coeff = Dict{Int64,Matrix{Float64}}()
	ids = Dict{Int64,Vector{Int64}}(i => Int64[] for i in 1:n) # for each agent index returns the indices of relevant monomials (according to 'keep')
	
	args = [(i,i2keep[i],keep,X,Y[[i,],:],d2i,λ,ρ,niter) for i in 1:n]
	pmap(mySINDy_par,args)

	for i in 1:n
		coeff[i] = readdlm("data/temp-sindy-Xi-$i.csv",',')
		rm("data/temp-sindy-Xi-$i.csv")
		ids[i] = Int64.(vec(readdlm("data/temp-sindy-id-$i.csv",',')))
		rm("data/temp-sindy-id-$i.csv")
	end

	return coeff,ids
end
#=
	for i in 1:n
		θ = ones(1,T)
		id = [1,]
		for k in i2keep[i]
			j = setdiff(keep[k],[i,])[1]
			θ = vcat(θ,reduce(vcat,[prod(X[v,:],dims=1) for v in [setdiff([j,l],[0,]) for l in 0:n]]))
			append!(id,[d2i[v] for v in [sort([j,l]) for l in 0:n]])
		end

		iii = unique(a -> id[a], eachindex(id))
		id = id[iii]
		θ = θ[iii,:]

		@info "============ mySINDy for agent $i =============="
		if isempty(id)
			coef = zeros(0,0)
		else
			coef,err,relerr = mySINDy_par(θ,Y,λ,ρ,niter)
			ids[i] = id
		end

		coeff[i] = coef
	end

	return coeff, ids, 1., 1.
end	
=#

function mySINDy_par(arg::Tuple{Int64,Vector{Int64},Vector{Vector{Int64}},Matrix{Float64},Matrix{Float64},Dict{Vector{Int64},Int64},Float64,Float64,Int64})
	i,i2keepi,keep,X,Y,d2i,λ,ρ,niter = arg
	
	n,T = size(X)
	θ = ones(1,T)
	id = [1,]
	for k in i2keepi
		j = setdiff(keep[k],[i,])[1]
		θ = vcat(θ,reduce(vcat,[prod(X[v,:],dims=1) for v in [setdiff([j,l],[0,]) for l in 0:n]]))
		append!(id,[d2i[v] for v in [sort([j,l]) for l in 0:n]])
	end

	iii = unique(a -> id[a], eachindex(id))
	id = id[iii]
	θ = θ[iii,:]
	m,T = size(θ)

	@info "-------------- mySINDy for agent $i ---------------"
	if isempty(id)
		coeff = zeros(0,0)
	else
		k = 1
		@info "Y: $(size(Y)), θ: $(size(θ))"
		Ξ = Y*θ'*inv(θ*θ' + ρ*Id(m))
		@info "Ξ: $(size(Ξ))"
		while k < niter
			k += 1
			smallinds = (abs.(Ξ) .< λ)
			Ξ[smallinds] .= 0.
			biginds = (.~smallinds)[1,:]
			@info "$(size(smallinds)), $(size(biginds))"
			Ξ[1,biginds] = Y*θ[biginds,:]'*inv(θ[biginds,:]*θ[biginds,:]' + ρ*Id(sum(biginds)))
		end
	end

	writedlm("data/temp-sindy-Xi-$i.csv",Ξ,',')
	writedlm("data/temp-sindy-id-$i.csv",id,',')
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
@info "Itertions done"

#	cost = sum((Y - Ξ*θ).^2) + ρ*sum(Ξ.^2)
	err = sum((Y - Ξ*θ).^2)
@info "Error computed"

	return Ξ, err, err/energy
end

#=
# Parallel version
function mySINDy_par(θ::Matrix{Float64}, Y::Matrix{Float64}, λ::Float64=.1, ρ::Float64=1., niter::Int64=10)
	n,T = size(Y)
	m,T = size(θ)
	energy = sum(Y.^2)

	nz = 1e6
	k = 0
	
	Ξ = 2*λ*ones(n,m)

	@info "n = $n"

	while k < niter# && sum(abs.(Ξ) .> 1e-6) != nz
		k += 1
		nz = sum(abs.(Ξ) .> 1e-6)
		err = sum((Y - Ξ*θ).^2)
		@info "iter $k: $(nz) nonzero coefficients, error = $(round(err)), rel-err = $(round(err/energy*100,digits=2))%"

		smallinds = (abs.(Ξ) .< λ)
		biginds = .~smallinds
		Ξ[smallinds] .= 0.
		
		args = [(Y[[i,],:],θ[biginds[i,:],:],ρ,"data/temp-sindy-$i.csv") for i in 1:n]
		pmap(sindy_step,args)

		for i in 1:n
			Ξ[i,biginds[i,:]] = readdlm("data/temp-sindy-$i.csv",',')
			rm("data/temp-sindy-$i.csv")
		end
	end
@info "Itertions done"

#	cost = sum((Y - Ξ*θ).^2) + ρ*sum(Ξ.^2)
	err = sum((Y - Ξ*θ).^2)
@info "Error computed"
@info "$(typeof(Ξ))"

	return Ξ, err, err/energy
end

function sindy_step(args::Tuple{Matrix{Float64},Matrix{Float64},Float64,String})
	Y,θ,ρ,file = args
	n = size(θ)[1]
	writedlm(file,Y*θ'*inv(θ*θ' + ρ*Id(n)),',')
end
=#

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


function get_θd(X::Matrix{Float64}, dmax::Int64, i0::Int64=1)
	n,T = size(X)

	θ = ones(1,T)
	d = zeros(Int64,1,dmax)

	if dmax == 0
		return θ,d
	else
		for i in 1:n
			θ0,d0 = get_θd(X[i:n,:],dmax-1,i)
			#d0 += (i-1)*(d0 .> 0)
			θ = [θ; repeat(X[[i,],:],size(θ0)[1],1).*θ0]
			d = [d; d0 i*ones(Int64,size(d0)[1],1)]
			#d = [d; (d0.+i.-1) i*ones(Int64,size(d0)[1],1)]
		end

		d += (i0-1)*(d .> 0)
		
		return θ,d
	end
end

function get_θ(X::Matrix{Float64}, dmax::Int64, i0::Int64=1)
	n,T = size(X)

	θ = ones(1,T)

	if dmax == 0
		return θ
	else
		for i in 1:n
			θ0 = get_θ(X[i:n,:],dmax-1,i)
			θ = [θ; repeat(X[[i,],:],size(θ0)[1],1).*θ0]
		end

		return θ
	end
end

function get_d(n::Int64,dmax::Int64, i0::Int64=1)
	d = zeros(Int64,1,dmax)

	if dmax == 0
		return d
	else
		for i in 1:n
			d0 = get_d(n-i+1,dmax-1,i)
			d = [d;d0 i*ones(Int64,size(d0)[1],1)]
		end

		d += (i0-1)*(d .> 0)

		return d
	end
end




