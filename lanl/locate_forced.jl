using JuMP, AmplNLWriter, Ipopt, LinearAlgebra, DelimitedFiles, Dates, FFTW

include("res_dist.jl")

## Identification of the network Laplacian, using the covariance matrix method, proposed in Lokhov18

function system_identification_correl(X::Array{Float64,2}, dt::Float64, l::Float64=.01)
	@info "$(now()) -- Start..."
	
	nn,T = size(X)
	n = Int(nn/2)

	S0 = X[:,1:T-1]*X[:,1:T-1]' ./ (T-1)
	S1 = X[:,2:T]*X[:,1:T-1]' ./ (T-1)
	
	Ah = S1*inv(S0)
	Adh = (Ah - diagm(0 => ones(2*n)))./dt
	
	@info "$(now()) -- Stop."
	
	return Ah,Adh
end


function find_forcing_freq(X::Array{Float64,2}, dt::Float64)
	n,T = size(X)
	
	for i in 1:n
	end
end


## INPUT
# Xs: time series of measurements at the PMUs.
# pmu_idx: indices of nodes where measurements are available, note that the first one in the list will be used as reference.
# L: weighted Laplacian of the network.
# w: frequency of the forcing.
# n_fourier_modes: number of fourier modes kept to clean the data.

## OUTPUT
# err_id: Array of the errors (1st column) and index of the corresponding forced node.

function locate_forcing_full_gen(Xs::Array{Float64,2}, pmu_idx::Array{Int64,1}, L::Array{Float64,2}, w::Float64, T::Int64, dt::Float64, n_fourier_modes::Int64=5)
	n = size(L)[1]
	n_idx = length(pmu_idx)

 ##=	
	fXs = zeros(Complex{Float64},n_idx,T)
	fXst = zeros(Complex{Float64},n_idx,T)
	Xst = zeros(n_idx,T)
	for i in 1:n_idx
		fXs[i,:] = fft(Xs[i,:])
		ids = Int.(sortslices([norm.(fXs[i,2:end]) Array(2:T)],dims=1,rev=true)[1:2*n_fourier_modes,2])
		fXst[i,[1;ids]] = fXs[i,[1;ids]]
		Xst[i,:] = real.(ifft(fXst[i,:]))
	end
# =#
 #=
	Xst = Xs
# =#
	
	evv = eigen(L)
	us = evv.vectors
	ls = evv.values
	
	bs = (ls./(ls.^2 .+ w^2))*cos.(w*dt*(0:T-1)') - ((ls./(ls.^2 .+ w^2))*ones(1,T)).*(exp.(-ls*dt*(0:T-1)')) + (w./(ls.^2 .+ w^2))*sin.(w*dt*(0:T-1)')
	
	errs = Array{Float64,1}()
	
	print("|")
	for i in 1:n
		if i%10 == 0
			print("|")
		else
			print("-")
		end
				
		cs = (us[i,:]*ones(1,T)) .* bs
		
		err = 0.
		
		for j in 1:n_idx-1
			for k in j+1:n_idx
				jj = pmu_idx[j]
				kk = pmu_idx[k]
				Yj = us[jj,:]' * cs
				Yk = us[kk,:]' * cs
				for l in 1:n
					err += (Xst[j,l]*Yk[l] - Xst[k,l]*Yj[l])^2
				end
	#			err += norm(Xst[j,:].*Yk - Xst[k,:].*Yj)
			end
		end
		push!(errs,err)
	end
	print("\n")	
	
	err_id = sortslices([errs 1:n], dims=1)
	
	return err_id
end
		
	
		
function locate_forcing_slow(Xs::Array{Float64,2}, pmu_idx::Array{Int64,1}, L::Array{Float64,2})
	n = size(L)[1]
	n_idx,T = size(Xs)
	
	n_fourier_modes = 5
	
	fXs = zeros(Complex{Float64},n_idx,T)
	fXst = zeros(Complex{Float64},n_idx,T)
	Xst = zeros(n_idx,T)
	for i in 1:n_idx
		fXs[i,:] = fft(Xs[i,:])
		ids = Int.(sortslices([norm.(fXs[i,2:end]) Array(2:T)],dims=1,rev=true)[1:2*n_fourier_modes,2])
		fXst[i,[1;ids]] = fXs[i,[1;ids]]
		Xst[i,:] = real.(ifft(fXst[i,:]))
	end
	
	Xst = Xst - repeat(Xst[1,:]',n_idx,1)
	
	Ld = pinv(L)
	
	errs = Array{Float64,1}()
	
	print("|")
	for i in 1:n
		if i%10 == 0
			print("|")
		else
			print("-")
		end
		
		err = 0.
		for j in 2:n_idx-1
			for k in j+1:n_idx
				err += norm(Xst[j,:].*(Ld[i,pmu_idx[k]] - Ld[i,pmu_idx[1]]) - Xst[k,:].*(Ld[i,pmu_idx[j]] - Ld[i,pmu_idx[1]]))^2
			end
		end
		push!(errs,err)
	end
	print("\n")
	
	err_id = sortslices([errs 1:n], dims=1)
	
	return err_id
end

function locate_forcing_ipopt(X::Array{Float64,2}, A::Array{Float64,2}, dt::Float64, nf::Int64 = 1, l::Float64=.1)
	@info "$(now())"," -- Start..."
	
	nn,T = size(X)
	n = Int(nn/2)
	
	AX = [zeros(n,T); A[n+1:2*n,:]*X]
	
#	locate_f = Model(with_optimizer(Ipopt.Optimizer, mumps_mem_percent=5))	
	locate_f = Model(with_optimizer(AmplNLWriter.Optimizer, "ipopt"))

## Variables
	@variable(locate_f, 0 <= c[i = 1:n] <= 5)
	@variable(locate_f, f[i = 1:n])
	@variable(locate_f, 0 <= phi[i = 1:n] <= 2*pi)
	
## Constraints
#	@NLconstraint(locate_f, sum((c[i] > 0) for i = 1:n) == nf)
	
## Objective
	@NLobjective(locate_f, Min, sum((X[n+i,t+1] - AX[n+i,t] - c[i]*cos(f[i]*t + phi[i]))^2 for i = 1:n for t = 1:T-1) + l*sum(c[i] for i=1:n))
#	@NLobjective(locate_f, Min, sum((X[n+i,t+1] - AX[n+i,t] - c[i]*(1-(f[i]*t + phi[i])*(1-f[i]*t + phi[i])))*(X[n+i,t+1] - AX[n+i,t] - c[i]*(1-(f[i]*t + phi[i])*(1-f[i]*t + phi[i]))) for i = 1:n for t = 1:T-1) + l*sum(c[i] for i = 1:n))
		
	JuMP.optimize!(locate_f)
	
	@info "$(now())"," -- Stop."
	
	return (
		c = value.(c),
		f = value.(f),
		phi = value.(phi)
	)
end





## Identification of the network Laplacian by an optimization formulation (Ipopt).
## !!!!!!!!!!!!!!!!!!! FOR TOO LARGE NETWORK/DATASET, IT FILLS THE RAM AND CRASHES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!¨

function system_identification_ipopt(X::Array{Float64,2}, m::Array{Float64,1}, d::Array{Float64,1}, dt::Float64)
	@info "$(now()) -- Start..."
	
	nn,T = size(X)
	n = Int(nn/2)
	
	Mi = diagm(0 => 1 ./ m)
	D = diagm(0 => d)

#	system_id = Model(with_optimizer(Ipopt.Optimizer, mumps_mem_percent=5))
	system_id = Model(with_optimizer(Ipopt.Optimizer))
		
## Variables
	@variable(system_id, L[i = 1:n, j = 1:n])
	
## Constraints
	for i in 1:n-1
		for j in i+1:n
			@constraint(system_id, L[i,j] <= 0)
			@constraint(system_id, L[i,j] == L[j,i])
		end
	end
	
	for i in 1:n
#		@constraint(system_id, sum(L[i,:]) == 0)
		@constraint(system_id, L[i,i] == -sum(L[i,j] for j = 1:i-1) - sum(L[i,j] for j = i+1:n))
	end
	
## Objective
	@NLexpression(system_id, err[i = 1:n, t = 1:T-1], X[n+i,t+1] - X[n+i,t] - dt * (sum(-L[i,k]*X[k,t] for k = 1:n)/m[i] - d[i]*X[n+i,t]/m[i]))
	
	@NLobjective(system_id, Min, sum(err[i,t]^2 for i in 1:n for t in 1:T-1)) 

	optimize!(system_id)
	
	@info "$(now()) -- Stop."
	
	return value.(L)
end

## Identification of the network Laplacian by an optimization formulation (Ipopt).
## !!!!!!!!!!!!!!!!!!! FOR TOO LARGE NETWORK/DATASET, IT FILLS THE RAM AND CRASHES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!¨



function system_identification_ipopt(X::Array{Float64,2}, dt::Float64)
	@info "$(now()) -- Start..."
	
	nn,T = size(X)
	n = Int(nn/2)
	
	l = .01
		
	system_id = Model(with_optimizer(Ipopt.Optimizer, print_level=2, mumps_mem_percent=5))
	
## Variables
	@variables(system_id, begin
		L[i = 1:n, j = 1:n]
		m[i = 1:n] >= 0
		d[i = 1:n] >= 0
#		c[i = 1:n] 
#		f[i = 1:n] 
#		phi[i = 1:n]
	end)
	
## Constraints
	for i in 1:n-1
		for j in i+1:n
			@constraint(system_id, L[i,j] <= 0)
			@constraint(system_id, L[i,j] == L[j,i])
		end
	end
	
	for i in 1:n
		@constraint(system_id, sum(L[i,:]) == 0)
	end
	
## Objective
	@NLexpression(system_id, err[i = 1:n, t = 1:T-1], X[n+i,t+1] - X[n+i,t] - dt * (sum(-L[i,k]*X[k,t] for k = 1:n)/m[i] - d[i]*X[n+i,t]/m[i]))
	
	@NLobjective(system_id, Min, sum(err[i,t]^2 for i = 1:n for t = 1:T-1) - l * sum(L[i,j] for i = 1:n-1 for j = i+1:n)) 

	
	optimize!(system_id)
	
	@info "$(now()) -- Stop."
	
	return (
		L = value.(L),
		m = value.(m),
		d = value.(d),
#		c = value.(c),
#		f = value.(f),
#		phi = value.(phi)
	)
end
	




