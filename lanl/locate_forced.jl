using JuMP, Ipopt, LinearAlgebra, DelimitedFiles, Dates, FFTW

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

function locate_forcing_ipopt(X::Array{Float64,2}, A::Array{Float64,2}, dt::Float64, l::Float64 = .1)
	@info "$(now())"," -- Start..."
	
	nn,T = size(X)
	n = Int(nn/2)
	
	AX = A*X
	
	locate_f = Model(with_optimizer(Ipopt.Optimizer))	

## Variables
	@variable(locate_f, 0 <= a[i = 1:n] <= 10)
	@variable(locate_f, f[i = 1:n])
	@variable(locate_f, 0 <= phi[i = 1:n] <= 2*pi)
	
## Objective
@NLobjective(locate_f, Min, sum((X[n+i,t+1] - AX[n+i,t] - dt*a[i]*cos(f[i]*2*pi*t*dt + phi[i]))*(X[n+i,t+1] - AX[n+i,t] - dt*a[i]*cos(f[i]*2*pi*t*dt + phi[i])) for i in 1:n for t in 1:T-1) + l*sum(a[i] for i in 1:n))
	
	JuMP.optimize!(locate_f)
	
	@info "$(now())"," -- Stop."

 	
	return (
		a = value.(a),
		f = value.(f),
		phi = value.(phi)
	)

end


function location_freq_forcing_ipopt(X::Array{Float64,2}, A::Array{Float64,2}, dt::Float64, nf::Int64 = 1, l::Float64 = 1.)
	@info "$(now())"," -- Start..."
	
	nn,T = size(X)
	n = Int(nn/2)
	
	AX = A*X
	
	locate_f = Model(with_optimizer(Ipopt.Optimizer))

## Variables
	@variable(locate_f, 0 <= c[i = 1:n] <= 10)
	@variable(locate_f, f[i = 1:n])
#	@variable(locate_f, 0 <= phi[i = 1:n] <= 2*pi)
	
## Constraints
#	@NLconstraint(locate_f, sum((c[i] > 0) for i = 1:n) == nf)
#=	for i in 2:n
		@constraint(locate_f, c[i] == 0)
	end
=#
	@NLexpression(locate_f, co[i,t] = cos(f[i]*2*pi*t*dt) for i in 1:n for t in 1:T-1)
## Objective
	@NLobjective(locate_f, Min, sum((X[n+i,t+1] - AX[n+i,t] - dt*c[i]*co[i,t])*(X[n+i,t+1] - AX[n+i,t] - dt*c[i]*co[i,t]) for i in 1:n for t in 1:T-1) + l*sum(c[i] + f[i] for i in 1:n))
	
#	@NLobjective(locate_f, Min, sum((X[n+i,t+1] - AX[n+i,t] - c[i]*cos(f[i]*t + phi[i]))^2 for i = 1:n for t = 1:T-1) + l*sum(c[i] for i=1:n))
#	@NLobjective(locate_f, Min, sum((X[n+i,t+1] - AX[n+i,t] - c[i]*(1-(f[i]*t + phi[i])*(1-f[i]*t + phi[i])))*(X[n+i,t+1] - AX[n+i,t] - c[i]*(1-(f[i]*t + phi[i])*(1-f[i]*t + phi[i]))) for i = 1:n for t = 1:T-1) + l*sum(c[i] for i = 1:n))
		
	JuMP.optimize!(locate_f)
	
	@info "$(now())"," -- Stop."

 	
	return (
		c = value.(c),
		f = value.(f),
#		phi = value.(phi)
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



function system_identification_ipopt(X::Array{Float64,2}, dt::Float64, l::Float64=.01)
	@info "$(now()) -- Start..."
	
	nn,T = size(X)
	n = Int(nn/2)
		
	system_id = Model(with_optimizer(Ipopt.Optimizer))
	
## Variables
	@variables(system_id, begin
		Lm[i = 1:n, j = 1:n]
		dm[i = 1:n] >= 0
		a[i = 1:n] >= 0 
		f[i = 1:n] >= 0
		0 <= phi[i = 1:n] <= 2pi
	end)

 ##=
## Constraints
	for i in 1:n-1
		for j in i+1:n
			@constraint(system_id, Lm[i,j] <= 0)
			@constraint(system_id, Lm[j,i] <= 0)
		end
	end
	
	for i in 1:n
		@constraint(system_id, sum(Lm[i,:]) == 0)
	end
# =#

## Objective
@NLexpression(system_id, err[i = 1:n, t = 1:T-1], X[n+i,t+1] - X[n+i,t] - dt * (sum(-Lm[i,k]*X[k,t] for k = 1:n) - dm[i]*X[n+i,t] - a[i]*cos(2*pi*f[i]*t*dt + phi[i])))

@NLobjective(system_id, Min, sum(err[i,t]^2 for i = 1:n for t = 1:T-1) - l * sum(Lm[i,j] + Lm[j,i] for i = 1:n-1 for j = i+1:n)) 
	
	optimize!(system_id)
	
	@info "$(now()) -- Stop."
	
	return (
		Lm = value.(Lm),
		dm = value.(dm),
		a = value.(a),
		f = value.(f),
		phi = value.(phi)
	)
end
	

function system_identification_ipopt_fullA(X::Array{Float64,2}, dt::Float64, l::Float64=.01)
	@info "$(now()) -- Start..."
	
	nn,T = size(X)
	n = Int(nn/2)
		
	system_id = Model(with_optimizer(Ipopt.Optimizer))
	
## Variables
	@variables(system_id, begin
		A[i = n+1:2*n, j = 1:2*n]
		a[i = 1:n] >= 0 
		f[i = 1:n] >= 0
		0 <= phi[i = 1:n] <= 2pi
	end)

## Constraints
	for i in 1:n-1
		for j in i+1:n
			@constraint(system_id, A[i+n,j] >= 0)
			@constraint(system_id, A[j+n,i] >= 0)
			@constraint(system_id, A[i+n,j+n] == 0)
			@constraint(system_id, A[j+n,i+n] == 0)
		end
	end
	
	for i in 1:n
		@constraint(system_id, sum(A[i+n,k] for k = 1:n) == 0)
		@constraint(system_id, A[i+n,i+n] >= 0)
	end

## Objective
@NLexpression(system_id, err[i = 1:n, t = 1:T-1], X[n+i,t+1] - X[n+i,t] - sum(A[i+n,k]*X[k,t] for k = 1:2*n) - dt * (a[i]*cos(2*pi*f[i]*t*dt + phi[i])))

@NLobjective(system_id, Min, sum(err[i,t]^2 for i = 1:n for t = 1:T-1) - l * sum(A[i+n,j] + A[j+n,i] for i = 1:n-1 for j = i+1:n)) 
@info "Objective done."
	
	optimize!(system_id)
	
	@info "$(now()) -- Stop."
	
	return (
		A = value.(A),
		a = value.(a),
		f = value.(f),
		phi = value.(phi)
	)
end
	

function fake_l0(Xs::Array{Float64,2}, dt::Float64, l::Float64=.01)
	nn,T = size(Xs)
	n = Int(nn/2)

	sols = Array{Any,1}()
	objs = Array{Float64,1}()

	for i in 1:n
		@info "i = $i"

		sol,obj = single_for_ipopt(i,Xs,dt,l)

		push!(sols,sol)
		push!(objs,obj)
	end

	objective,id = findmin(objs)

	return id,sols[id],objective
end

function single_for_ipopt(id::Int64, Xs::Array{Float64,2}, dt::Float64, l::Float64=.01)
	nn,T = size(Xs)
	n = Int(nn/2)

	system_id = Model(with_optimizer(Ipopt.Optimizer))
	
## Variables
	@variables(system_id, begin
		Lm[i = 1:n, j = 1:n]
		dm[i = 1:n] >= 0
		a >= 0 
		f >= 0
		0 <= phi <= 2pi
	end)

 ##=
## Constraints
	for i in 1:n-1
		for j in i+1:n
			@constraint(system_id, Lm[i,j] <= 0)
			@constraint(system_id, Lm[j,i] <= 0)
		end
	end
	
	for i in 1:n
		@constraint(system_id, sum(Lm[i,:]) == 0)
	end
# =#

## Objective
	@NLexpression(system_id, err[i = [1:id-1;id+1:n], t = 1:T-1], Xs[n+i,t+1] - Xs[n+i,t] - dt * (sum(-Lm[i,k]*Xs[k,t] for k = 1:n) - dm[i]*Xs[n+i,t]))
	@NLexpression(system_id, erri[t = 1:T-1], Xs[n+id,t+1] - Xs[n+id,t] - dt * (sum(-Lm[id,k]*Xs[k,t] for k = 1:n) - dm[id]*Xs[n+id,t] - a*cos(2*pi*f*t*dt + phi)))

	@NLobjective(system_id, Min, sum(err[i,t]^2 for i = [1:id-1;id+1:n] for t = 1:T-1) + sum(erri[t]^2 for t = 1:T-1) - l * sum(Lm[i,j] + Lm[j,i] for i = 1:n-1 for j = i+1:n)) 
	
	optimize!(system_id)
	
	@info "$(now()) -- Stop."
	
	return (
		Lm = value.(Lm),
		dm = value.(dm),
		a = value(a),
		f = value(f),
		phi = value(phi)
		), objective_value(system_id)

end


function system_identification_iterate(Xs::Array{Float64,2}, dt::Float64, l::Float64=.01)
	nn,T = size(Xs)
	n = Int(nn/2)

	Ah,Adh = system_identification_correl(Xs,dt,l)
	
	frcs = Array{Tuple{Array{Float64,1},Array{Float64,1},Array{Float64,1}},1}()
	Ahs = Array{Array{Float64,2},1}()
	push!(Ahs,Ah)

	c = 0
	while c < 10
		c += 1
		@info "iter = $c"

		frc = get_frc(Xs,Ah,dt,l)
		Ah = get_A(Xs,frc,dt,l)

		push!(frcs,frc)
		push!(Ahs,Ah)
	end

	return Ahs,frcs
end

function get_frc(Xs::Array{Float64,2}, Ah::Array{Float64,2}, dt::Float64, l::Float64=.01)
	nn,T = size(Xs)
	n = Int(nn/2)
	
	AX = Ah*Xs

	frc_id = Model(with_optimizer(Ipopt.Optimizer))
	
## Variables
	@variables(frc_id, begin
		a[i = 1:n] >= 0 
		f[f = 1:n] >= 0
		0 <= phi[i = 1:n] <= 2pi
	end)

## Objective
	@NLexpression(frc_id, err[i = 1:n, t = 1:T-1], Xs[n+i,t+1] - AX[n+i,t] - dt*a[i]*cos(2*pi*f[i]*t*dt + phi[i]))

	@NLobjective(frc_id, Min, sum(err[i,t]^2 for i = [1:id-1;id+1:n] for t = 1:T-1) + l * sum(a[i] for i = 1:n-1 for j = i+1:n)) 
	
	optimize!(frc_id)
	
	return (value.(a), value.(f), value.(phi))
end

function get_A(Xs::Array{Float64,2}, frc::Tuple{Array{Float64,1},Array{Float64,1},Array{Float64,1}}, dt::Float64, l::Float64=.01)
	nn,T = size(Xs)
	n = Int(nn/2)
	
	a,f,phi = frc

	A_id = Model(with_optimizer(Ipopt.Optimizer))

	@variables(A_id, begin
			   Lm[i = 1:n, j = 1:n]
			   dm[i = 1:n]
		   end)

	for i in 1:n-1
		for j in i+1:n
			@constraint(A_id, Lm[i,j] <= 0)
			@constraint(A_id, Lm[j,i] <= 0)
		end
	end
	
	for i in 1:n
		@constraint(A_id, sum(Lm[i,:]) == 0)
	end
	
	@objective(A_id, Min, sum((Xs[n+i,t+1] - Xs[n+i,t] - dt * (sum(-Lm[i,k]*Xs[k,t] for k = 1:n) - dm[i]*Xs[n+i,t]) - a[i]*cos(2*pi*f[i]*dt*t + phi[i]))*(Xs[n+i,t+1] - Xs[n+i,t] - dt * (sum(-Lm[i,k]*Xs[k,t] for k = 1:n) - dm[i]*Xs[n+i,t] - a[i]*cos(2*pi*f[i]*dt*t + phi[i]))) for i = 1:n for t = 1:T-1) + l * sum(Lm[i,j] + Lm[j,i] for i in 1:n-1 for j in i+1:n))
	
	optimize!(A_id)

	Id = diagm(0 => ones(n))

	return [Id dt*Id; -dt*value.(Lm) -dt*diagm(0 => value.(dm))]
end







