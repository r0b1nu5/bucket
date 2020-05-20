using Distributed

n_thr = 4

if nworkers() < n_thr
	addprocs(n_thr - nworkers())
end

@everywhere include("final_new.jl")

@everywhere function locate_new_l0(id::String, Xs::Array{Float64,2}, tau::Float64, ks::Tuple{Int64,Int64,Int64}, ls::Tuple{Int64,Int64,Int64}, mu::Float64=1e-1, bp::Float64=1e-1)
	nn,NN = size(Xs)
	n = Int.(n/2)
	N = NN - 1

	kmin,kmax,dk = ks
	lmin,lmax,dl = ls

# Computing the needed inputs (time series, discrete derivative, and their Fourier transforms).
	x = Xs[:,1:end-1]
	Dx = (Xs[(n+1):(2*n),2:end] - Xs[(n+1):(2*n),1:end-1])/tau
	xt = Array{Complex{Float64},2}(undef,nn,N)
	for i in 1:nn
		xt[i,:] = ifft(x[i,:])*sqrt(N)
	end
	Dxt = Array{Complex{Float64},2}(undef,n,N)
	for i in 1:n
		Dxt[i,:] = ifft(Dx[i,:])*sqrt(N)
	end

	for k in kmin:dk:kmax
		for l in lmin:dl:lmax
			push!(args,(id,x,Dx,xt,Dxt,l,k,mu,bp))
		end
	end

	pmap(Lmin_l0_par,args)
end

@everywhere function locate_new_l2(id::String, Xs::Array{Float64,2}, tau::Float64, ks::Tuple{Int64,Int64,Int64}, mu::Float64=1e-1, bp::Float64=1e-1)
	nn,NN = size(Xs)
	n = Int.(n/2)
	N = NN - 1

	kmin,kmax,dk = ks

# Computing the needed inputs (time series, discrete derivative, and their Fourier transforms).
	x = Xs[:,1:end-1]
	Dx = (Xs[(n+1):(2*n),2:end] - Xs[(n+1):(2*n),1:end-1])/tau
	xt = Array{Complex{Float64},2}(undef,nn,N)
	for i in 1:nn
		xt[i,:] = ifft(x[i,:])*sqrt(N)
	end
	Dxt = Array{Complex{Float64},2}(undef,n,N)
	for i in 1:n
		Dxt[i,:] = ifft(Dx[i,:])*sqrt(N)
	end

	for k in kmin:dk:kmax
		push!(args,(id,x,Dx,xt,Dxt,k,mu,bp))
	end

	pmap(Lmin_l2_par,args)
end

"""
	Lmin_l0_par(id::String x::Array{Float64,2}, Dx::Array{Float64,2}, xt::Array{Complex{Float64},2}, Dxt::Array{Complex{Float64,2}, l::Int64, k::Int64, mu::Float64=1e-1, bp::Float64=1e-1)

Minimizes the quadratic error in the estimation of the forced trajectory, for a fixed frequency (k) and location (l) of the forcing. The optimization parameters are the dynamics matrix (A1), the damings (a2), and the forcing amplitude (gamma). 

_INPUT_:
`id`: Name of the time series.
`x`: Time series of the phase angles.
`Dx`: Time series of the phase frequencies. 
`xt`: (Inverse) Fourier transform of x.
`Dxt`: (Inverse) Fourier transform of Dx.
`l`: Fixed location of the forcing.
`k`: Fixed index of the forcing frequency (ν = k/T).
`mu`: Initial value of the barrier parameter (in IPOPT).
`bp`: Initial value of the bound_push parameter (in IPOPT).

_OUTPUT_:
`objective`: Value of the optimized objective function.
`A1`: Best estimate of the dynamcis matrix.
`a2`: Best estimate of the dampings.
`gamma`: Best estimate of the forcing amplitude.
"""
function Lmin_l0(id::String, x::Array{Float64,2}, Dx::Array{Float64,2}, xt::Array{Complex{Float64},2}, Dxt::Array{Complex{Float64},2}, l::Int64, k::Int64, mu::Float64=1e-1, bp::Float64=1e-1)
	nn,N = size(x)
	n = Int(nn/2)

# Definition of the needed parameters.
	Sigma0 = (x*x')/N
	Sigma1 = (x*Dx')/N

	xtk = xt[:,k]
	Dxtlk = Dxt[l,k]
	
	Fk = real.(xtk*xtk')
	flk = real.(Dxtlk*xtk')
	glk = norm(Dxtlk)^2

# Definition of the optimization problem.
	system_id = Model(with_optimizer(Ipopt.Optimizer, mu_init = mu, bound_push = bp))

	@variable(system_id, A1[i = 1:n, j = 1:n])
	for i in 1:n-1
		for j in i+1:n
			@constraint(system_id, A1[i,j] == A1[j,i])
		end
	end
	@variable(system_id, a2[i = 1:n])
	for i in 1:n
		@constraint(system_id, a2[i] >= 0.)
	end
	@variable(system_id, gamma >= 0.)

	@NLexpression(system_id, AtAS0[i = 1:n], sum(A1[j1,i]*A1[j1,j2]*Sigma0[j2,i] for j1 = 1:n for j2 = 1:n) + sum(A1[j,i]*a2[j]*Sigma0[n+j,i] for j = 1:n))
	@NLexpression(system_id, AtAS00[i = (n+1):(2*n)], sum(a2[i-n]*A1[i-n,j]*Sigma0[j,i] for j = 1:n) + a2[i-n]*a2[i-n]*Sigma0[i,i])
	@NLexpression(system_id, T1, sum(AtAS0[i] for i in 1:n) + sum(AtAS00[i] for i = (n+1):(2*n)))

	@NLexpression(system_id, AS1[i = 1:n], sum(A1[i,j]*Sigma1[j,i] for j = 1:n) + a2[i]*Sigma1[n+i,i])
	@NLexpression(system_id, T2, sum(AS1[i] for i in 1:n))

	@NLexpression(system_id, AAF[i = 1:n], A1[l,i]*(sum(A1[l,j]*Fk[j,i] for j = 1:n) + a2[l]*Fk[n+l,i]))
	@NLexpression(system_id, AAFF, a2[l]*(sum(A1[l,j]*Fk[j,n+l] for j = 1:n) + a2[l]*Fk[n+l,n+l]))
	@NLexpression(system_id, T3, sum(AAF[i] for i = 1:n) + AAFF)

	@NLexpression(system_id, T4, sum(flk[i]*A1[l,i] for i in 1:n) + flk[n+l]*a2[l])

	@NLexpression(system_id, g2, gamma^2)

	@NLobjective(system_id, Min, T1 + 2*T2 + .5*g2 - 2*gamma/sqrt(N)*sqrt(T3 + 2*T4 + glk))

	optimize!(system_id)

	writedlm("data/"*id*"_l0_$(l).$(k)_obj.csv",objective_value(system_id),',')
	writedlm("data/"*id*"_l0_$(l).$(k)_A1.csv",value.(A1),',')
	writedlm("data/"*id*"_l0_$(l).$(k)_a2.csv",value.(a2),',')
	writedlm("data/"*id*"_l0_$(l).$(k)_gamma.csv",value(gamma),',')
end

"""
	Lmin_l2_par(id::String, x::Array{Float64,2}, Dx::Array{Float64,2}, xt::Array{Complex{Float64},2}, Dxt::Array{Complex{Float64,2}, k::Int64, mu::Float64=1e-1, bp::Float64=1e-1)

Minimizes the quadratic error in the estimation of the forced trajectory, for a fixed frequency (k) of the forcing. The optimization parameters are the dynamics matrix (A1), the damings (a2), and the forcing amplitude (gamma). 

_INPUT_:
`id`: Name of the time series. 
`x`: Time series of the phase angles.
`Dx`: Time series of the phase frequencies. 
`xt`: (Inverse) Fourier transform of x.
`Dxt`: (Inverse) Fourier transform of Dx.
`k`: Fixed index of the forcing frequency (ν = k/T).
`mu`: Initial value of the barrier parameter (in IPOPT).
`bp`: Initial value of the bound_push parameter (in IPOPT).

_OUTPUT_:
`objective`: Value of the optimized objective function.
`A1`: Best estimate of the dynamcis matrix.
`a2`: Best estimate of the dampings.
`gamma`: Best estimate of the forcing amplitude at each possible location.
"""
function Lmin_l2_par(id::String, x::Array{Float64,2}, Dx::Array{Float64,2}, xt::Array{Complex{Float64},2}, Dxt::Array{Complex{Float64},2}, k::Int64, mu::Float64=1e-1, bp::Float64=1e-1)
	nn,N = size(x)
	n = Int(nn/2)

# Definition of the needed parameters.
	Sigma0 = (x*x')/N
	Sigma1 = (x*Dx')/N

	xtk = xt[:,k]
	
	Fk = real.(xtk*xtk')
	flk = [real.(Dxt[l,k]*xtk') for l = 1:n]
	glk = [norm(Dxt[l,k])^2 for l = 1:n]

# Definition of the optimization problem. 
	system_id = Model(with_optimizer(Ipopt.Optimizer, mu_init = mu, bound_push = bp))

	@variable(system_id, A1[i = 1:n, j = 1:n])
	for i in 1:n-1
		for j in i+1:n
			@constraint(system_id, A1[i,j] == A1[j,i])
		end
	end
	@variable(system_id, a2[i = 1:n])
	@variable(system_id, gamma[i = 1:n])
	for i in 1:n
		@constraint(system_id, a2[i] >= 0.)
		@constraint(system_id, gamma[i] >= 0.)
	end

	@NLexpression(system_id, AtAS0[i = 1:n], sum(A1[j1,i]*A1[j1,j2]*Sigma0[j2,i] for j1 = 1:n for j2 = 1:n) + sum(A1[j,i]*a2[j]*Sigma0[n+j,i] for j = 1:n))
	@NLexpression(system_id, AtAS00[i = (n+1):(2*n)], sum(a2[i-n]*A1[i-n,j]*Sigma0[j,i] for j = 1:n) + a2[i-n]*a2[i-n]*Sigma0[i,i])
	@NLexpression(system_id, T1, sum(AtAS0[i] for i in 1:n) + sum(AtAS00[i] for i = (n+1):(2*n)))

	@NLexpression(system_id, AS1[i = 1:n], sum(A1[i,j]*Sigma1[j,i] for j = 1:n) + a2[i]*Sigma1[n+i,i])
	@NLexpression(system_id, T2, sum(AS1[i] for i in 1:n))

	@NLexpression(system_id, g2, sum(gamma[i]^2 for i = 1:n))

	@NLexpression(system_id, AAF[l = 1:n, i = 1:n], A1[l,i]*(sum(A1[l,j]*Fk[j,i] for j = 1:n) + a2[l]*Fk[n+l,i]))
	@NLexpression(system_id, AAFF[l = 1:n], a2[l]*(sum(A1[l,j]*Fk[j,n+l] for j = 1:n) + a2[l]*Fk[n+l,n+l]))
	@NLexpression(system_id, T3[l = 1:n], sum(AAF[l,i] for i = 1:n) + AAFF[l])

	@NLexpression(system_id, T4[l = 1:n], sum(flk[l][i]*A1[l,i] for i in 1:n) + flk[l][n+l]*a2[l])

	@NLobjective(system_id, Min, T1 + 2*T2 + .5*g2 - 2/sqrt(N)*sum(gamma[l]*sqrt(T3[l] + 2*T4[l] + glk[l]) for l = 1:n))

	optimize!(system_id)

	writedlm("data/"*id*"_l2_$(l).$(k)_obj.csv",objective_value(system_id),',')
	writedlm("data/"*id*"_l2_$(l).$(k)_A1.csv",value.(A1),',')
	writedlm("data/"*id*"_l2_$(l).$(k)_a2.csv",value.(a2),',')
	writedlm("data/"*id*"_l2_$(l).$(k)_gamma.csv",value.(gamma),',')
end





