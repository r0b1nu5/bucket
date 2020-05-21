using DelimitedFiles, PyPlot, LinearAlgebra, JuMP, Ipopt, Distributed, FFTW

"""
	run_new_l0(Xs::Array{Float64,2}, tau::Float64, ks::Tuple{Int64,Int64,Int64}, plot::Bool=false, mu::Float64=1e-1, bp::Float64=1e-1)

Identifies dynamics and forcing characteristics baded only on measurements. Approximation of the forcing's frequency is discretized as 2*pi*k/T, with k integer.

_INPUT_:
`Xs`: Time series of the phase angles (rows 1:n) and of the phase frequencies (rows n+1:2*n).
`tau`: Time step size.
`ks = (kmin,kmax,dk)`: Values of k to be tried, kmin:dk:kmax.
`plot`: If true, generates the plots of the objective function vs. k.
`mu`: Initial value of the barrier parameter (in IPOPT).
`bp`: Initial value of the bound_push parameter (in IPOPT).

_OUTPUT_:
`Ls_l0`: values of the objective function for the various values of k and l, under l0.
`(L_l0, A_l0, d_l0, gamma_l0, k_l0, l_l0)`: results under l0.
	`L_l0`: Minimal value of the objective.
	`A_l0`: Estimate of the dynamics matrix.
	`d_l0`: Estimate of the damings.
	`gamma_l0`: Estimate of the forcing amplitude.
	`k_l0`: Estimate frequency index (see theory).
	`l_l0`: Estimate of the forcing location.
"""
function run_new_l0(Xs::Array{Float64,2}, tau::Float64, ks::Tuple{Int64,Int64,Int64}, plot::Bool=false, mu::Float64=1e-1, bp::Float64=1e-1)
	nn,NN = size(Xs)
	n = Int(nn/2)
	N = NN-1

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

# Compute warm start
	XXX,A1h,a2h = get_Ah_correl_new(Xs,tau)

# Run the optimizations
	Ls_l0 = zeros(n,length(kmin:dk:kmax))
	L_l0 = 1000.
	A_l0 = zeros(n,nn)
	d_l0 = zeros(n)
	gamma_l0 = 0.
	l_l0 = 0
	k_l0 = 0
	for l in 1:n
		c = 0
		for k in kmin:dk:kmax
			c += 1
			@info "l0: l = $l, k = $k"
			Lt = Lmin_l0(x,Dx,xt,Dxt,l,k,A1h,a2h,mu,bp)
			Ls_l0[l,c] = Lt[1]
			if Lt[1] < L_l0
				L_l0,A_l0,d_l0,gamma_l0 = Lt
				l_l0 = l
				k_l0 = k
			end
		end
	end

	@info "Best l0: l = $(l_l0), ω = $(2*pi*(k_l0-1)/(tau*N)), L = $(L_l0)."	

	if plot
		T = N*tau
		figure("l0, T = $(N*tau)")
		plot_new_l0(Ls_l0, L_l0, l_l0, k_l0, ks, T, n)
	end

	return Ls_l0, (L_l0,A_l0,d_l0,gamma_l0,l_l0,k_l0)
end


"""
	run_new_l2(Xs::Array{Float64,2}, tau::Float64, ks::Tuple{Int64,Int64,Int64}, plot::Bool=false, mu::Float64=1e-1, bp::Float64=1e-1)

Identifies dynamics and forcing characteristics baded only on measurements. Approximation of the forcing's frequency is discretized as 2*pi*k/T, with k integer.

_INPUT_:
`Xs`: Time series of the phase angles (rows 1:n) and of the phase frequencies (rows n+1:2*n).
`tau`: Time step size.
`ks = (kmin,kmax,dk)`: Values of k to try, (kmin:dk:kmax).
`plot`: If true, generates the plots of the objective function vs. k.
`mu`: Initial value of the barrier parameter (in IPOPT).
`bp`: Initial value of the bound_push parameter (in IPOPT).

_OUTPUT_:
`Ls_l2`: values of the objective function for the various values of k and l, under l2.
`(L_l2, A_l2, d_l2, gamma_l2, k_l2, l_l2)`: results under l2.
	`L_l2`: Minimal value of the objective.
	`A_l2`: Estimate of the dynamics matrix.
	`d_l2`: Estimate of the damings.
	`gamma_l2`: Estimate of the forcing amplitude.
	`k_l2`: Estimate frequency index (see theory).
	`l_l2`: Estimate of the forcing location.
"""
function run_new_l2(Xs::Array{Float64,2}, tau::Float64, ks::Tuple{Int64,Int64,Int64}, plot::Bool=false, mu::Float64=1e-1, bp::Float64=1e-1)
	nn,NN = size(Xs)
	n = Int(nn/2)
	N = NN-1

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

# Compute warm start
	XXX,A1h,a2h = get_Ah_correl_new(Xs,tau)

# Run the optimizations
	Ls_l2 = zeros(length(kmin:dk:kmax))
	L_l2 = 1000.
	A_l2 = zeros(n,nn)
	d_l2 = zeros(n)
	gamma_l2 = zeros(n)
	k_l2 = 0
	l_l2 = 0
	c = 1
	for k in kmin:dk:kmax
		c += 1
		@info "l2: k = $k"
		Lt = Lmin_l2(x,Dx,xt,Dxt,k,A1h,a2h,mu,bp)
		Ls_l2[c] = Lt[1]
		if Lt[1] < L_l2
			L_l2,A_l2,d_l2,gamma_l2 = Lt
			k_l2 = k
			l_l2 = findmax(gamma_l2)[2]
		end
	end

	@info "Best l2: ω = $(2*pi*(k_l2-1)/(tau*N)), L = $(L_l2)."

	if plot
		T = N*tau
		figure("l2, T = $(N*tau)")
		subplot(1,2,1)
		plot_new_l2_1(Ls_l2, L_l2, gamma_l2, l_l2, k_l2, ks, T)
		subplot(1,2,2)
		plot_new_l2_2(Ls_l2, L_l2, gamma_l2, l_l2)
	end

	return Ls_l2, (L_l2,A_l2,d_l2,gamma_l2,l_l2,k_l2)
end


"""
	run_new(Xs::Array{Float64,2}, tau::Float64, kmax::Int64, plot::Bool=false, mu::Float64=1e-1, bp::Float64=1e-1)

Identifies dynamics and forcing characteristics baded only on measurements. Approximation of the forcing's frequency is discretized as 2*pi*k/T,  integer.

_INPUT_:
`Xs`: Time series of the phase angles (rows 1:n) and of the phase frequencies (rows n+1:2*n).
`tau`: Time step size.
`kmax`: Maximal value of k to be tried.
`plot`: If true, generates the plots of the objective function vs. k.
`mu`: Initial value of the barrier parameter (in IPOPT).
`bp`: Initial value of the bound_push parameter (in IPOPT).

_OUTPUT_:
`Ls_l0`: values of the objective function for the various values of k and l, under l0.
`(L_l0, A_l0, d_l0, gamma_l0, k_l0, l_l0)`: results under l0.
	`L_l0`: Minimal value of the objective.
	`A_l0`: Estimate of the dynamics matrix.
	`d_l0`: Estimate of the damings.
	`gamma_l0`: Estimate of the forcing amplitude.
	`k_l0`: Estimate frequency index (see theory).
	`l_l0`: Estimate of the forcing location.
`Ls_l2`: values of the objective function for the various values of k and l, under l2.
`(L_l2, A_l2, d_l2, gamma_l2, k_l2, l_l2)`: results under l2.
	`L_l2`: Minimal value of the objective.
	`A_l2`: Estimate of the dynamics matrix.
	`d_l2`: Estimate of the damings.
	`gamma_l2`: Estimate of the forcing amplitude.
	`k_l2`: Estimate frequency index (see theory).
	`l_l2`: Estimate of the forcing location.
"""
function run_new(Xs::Array{Float64,2}, tau::Float64, kmax::Int64, plot::Bool=false, mu::Float64=1e-1, bp::Float64=1e-1)
	nn,NN = size(Xs)
	n = Int(nn/2)
	N = NN-1

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

# Compute warm start
	XXX,A1h,a2h = get_Ah_correl_new(Xs,tau)

# Run the l0 optimizations
	Ls_l0 = zeros(n,kmax)
	L_l0 = 1000.
	A_l0 = zeros(n,nn)
	d_l0 = zeros(n)
	gamma_l0 = 0.
	l_l0 = 0
	k_l0 = 0
	for l in 1:n
		for k in 1:kmax
			@info "l0: l = $l, k = $k"
			Lt = Lmin_l0(x,Dx,xt,Dxt,l,k,A1h,a2h,mu,bp)
			Ls_l0[l,k] = Lt[1]
			if Lt[1] < L_l0
				L_l0,A_l0,d_l0,gamma_l0 = Lt
				l_l0 = l
				k_l0 = k
			end
		end
	end

	@info "Best l0: l = $(l_l0), ω = $(2*pi*(k_l0-1)/(tau*N)), L = $(L_l0)."
	
# Run the l0 optimizations
	Ls_l2 = zeros(kmax)
	L_l2 = 1000.
	A_l2 = zeros(n,nn)
	d_l2 = zeros(n)
	gamma_l2 = zeros(n)
	k_l2 = 0
	l_l2 = 0
	for k in 1:kmax
		@info "l2: k = $k"
		Lt = Lmin_l2(x,Dx,xt,Dxt,k,A1h,a2h,mu,bp)
		Ls_l2[k] = Lt[1]
		if Lt[1] < L_l2
			L_l2,A_l2,d_l2,gamma_l2 = Lt
			k_l2 = k
			l_l2 = findmax(gamma_l2)[2]
		end
	end

	@info "Best l2: ω = $(2*pi*(k_l2-1)/(tau*N)), L = $(L_l2)."

	if plot
		ws = Array(2*pi*(0:kmax-1)/(tau*N))
		
		figure("l0, T = $(N*tau)")
		plot_new_l0(Ls_l0, L_l0, l_l0, k_l0, ws, n)
		figure("l2, T = $(N*tau)")
		subplot(1,2,1)
		plot_new_l2_1(Ls_l2, L_l2, gamma_l2, l_l2, k_l2, ws)
		subplot(1,2,2)
		plot_new_l2_2(Ls_l2, L_l2, gamma_l2, l_l2, k_l2, ws)
	end

	return Ls_l0, (L_l0,A_l0,d_l0,gamma_l0,l_l0,k_l0), Ls_l2, (L_l2,A_l2,d_l2,gamma_l2,l_l2,k_l2)
end


"""
	Lmin_l0(x::Array{Float64,2}, Dx::Array{Float64,2}, xt::Array{Complex{Float64},2}, Dxt::Array{Complex{Float64,2}, l::Int64, k::Int64, A1h::Array{Float64,2}, a1h::Array{Float64,1}, mu::Float64=1e-1, bp::Float64=1e-1)

Minimizes the quadratic error in the estimation of the forced trajectory, for a fixed frequency (k) and location (l) of the forcing. The optimization parameters are the dynamics matrix (A1), the damings (a2), and the forcing amplitude (gamma). 

_INPUT_:
`x`: Time series of the phase angles.
`Dx`: Time series of the phase frequencies. 
`xt`: (Inverse) Fourier transform of x.
`Dxt`: (Inverse) Fourier transform of Dx.
`l`: Fixed location of the forcing.
`k`: Fixed index of the forcing frequency (ν = k/T).
`A1h`: Warm start for A1.
`a2h`: Warm start for a2.
`mu`: Initial value of the barrier parameter (in IPOPT).
`bp`: Initial value of the bound_push parameter (in IPOPT).

_OUTPUT_:
`objective`: Value of the optimized objective function.
`A1`: Best estimate of the dynamcis matrix.
`a2`: Best estimate of the dampings.
`gamma`: Best estimate of the forcing amplitude.
"""
function Lmin_l0(x::Array{Float64,2}, Dx::Array{Float64,2}, xt::Array{Complex{Float64},2}, Dxt::Array{Complex{Float64},2}, l::Int64, k::Int64, A1h::Array{Float64,2}, a2h::Array{Float64,1}, mu::Float64=1e-1, bp::Float64=1e-1)
	t0 = time()
	
	nn,N = size(x)
	n = Int(nn/2)

# Definition of the needed parameters.
	S0 = (x*x')/N
	S1 = (x*Dx')/N

	xtk = xt[:,k]
	Dxtlk = Dxt[l,k]
	
	Fk = real.(xtk*xtk')
	flk = real.(Dxtlk*xtk')
	glk = norm(Dxtlk)^2

# Definition of the optimization problem.
	system_id = Model(with_optimizer(Ipopt.Optimizer, mu_init = mu, bound_push = bp))

	@variable(system_id, A1[i = 1:n, j = i:n])
	for i = 1:n
		for j = i:n
			set_start_value(A1[i,j],A1h[i,j])
		end
	end
	@variable(system_id, a2[i = 1:n])
	set_start_value.(a2,a2h)
	@variable(system_id, gamma >= 0.)
	set_start_value(gamma,1.)

# tr(L^2*\th*\th^T)
	@NLexpression(system_id, AtAS01, 
		      sum(A1[i,j]*A1[j,k]*S0[k,i] for i = 1:n for j = i:n for k = j:n) + 
		      sum(A1[i,j]*A1[k,j]*S0[k,i] for i = 1:n for j = i:n for k = 1:j-1) + 
		      sum(A1[j,i]*A1[j,k]*S0[k,i] for i = 1:n for j = 1:i-1 for k = j:n) + 
		      sum(A1[j,i]*A1[k,j]*S0[k,i] for i = 1:n for j = 1:i-1 for k = 1:j-1))

# tr(-L*D*\om*\th^T) = tr(-D*L*\th*\om^T)
	@NLexpression(system_id, AtAS02, 
		      sum(A1[i,j]*a2[j]*S0[n+j,i] for i = 1:n for j = i:n) + 
		      sum(A1[j,i]*a2[j]*S0[n+j,i] for i = 1:n for j = 1:i-1))

# tr(D^2*\om*\om^T)
	@NLexpression(system_id, AtAS03, sum(a2[i]^2*S0[n+i,n+i] for i = 1:n))
	
# tr(A^T*A*S0)
	@NLexpression(system_id, T1, AtAS01 + 2*AtAS02 + AtAS03)


# tr(-L*\th*(D\om)^T)
	@NLexpression(system_id, AS11, 
		      sum(A1[i,j]*S1[j,i] for i = 1:n for j = i:n) + 
		      sum(A1[j,i]*S1[j,i] for i =1:n for j = 1:i-1))

# tr(-D*\om*(D\om)^T)
	@NLexpression(system_id, AS12, sum(a2[i]*S1[n+i,i] for i = 1:n))

# tr(A*S1)
	@NLexpression(system_id, T2, AS11 + AS12)


# tr((L_{l,:})^T*L_{l,:}*F(k))
	@NLexpression(system_id, AAF1, 
		      sum(A1[l,i]*A1[l,j]*Fk[j,i] for i = l:n for j = l:n) + 
		      2*sum(A1[i,l]*A1[l,j]*Fk[j,i] for i = 1:l-1 for j = l:n) + 
		      sum(A1[i,l]*A1[j,l]*Fk[j,i] for i = 1:l-1 for j = 1:l-1))

# tr((L_{l,:})^T*D_{l,:}*F(k)) = tr((D_{l,:})^T*L_{l,:}*F(k))
	@NLexpression(system_id, AAF2, 
		      sum(A1[l,i]*a2[l]*Fk[n+l,i] for i = l:n) + 
		      sum(A1[i,l]*a2[l]*Fk[n+l,i] for i = 1:l-1))

# tr((D_{l,:})^T*D_{l,:}*F(k))
	@NLexpression(system_id, AAF3, a2[l]^2*Fk[n+l,n+l])

# tr((A_{l,:})^T*A_{l,:}*F(k))
	@NLexpression(system_id, T3, AAF1 + 2*AAF2 + AAF3)

	
# (f_l(k))^T*A_{l,:}
	@NLexpression(system_id, T4, sum(flk[i]*A1[l,i] for i = l:n) + sum(flk[i]*A1[i,l] for i = 1:l-1) + flk[n+l]*a2[l])


	@NLexpression(system_id, g2, gamma^2)

	
	@NLobjective(system_id, Min, T1 + 2*T2 + .5*g2 - 2*gamma/sqrt(N)*sqrt(T3 + 2*T4 + glk))

	optimize!(system_id)

	mL = zeros(n,n)
	for i in 1:n
		for j in i+1:n
			mL[i,j] = value(A1[i,j])
			mL[j,i] = value(A1[i,j])
		end
		mL[i,i] = value(A1[i,i])
	end

	@info "Full optimization took $(time() - t0)''."

	return objective_value(system_id), mL, value.(a2), value(gamma)
end

"""
Lmin_l2(x::Array{Float64,2}, Dx::Array{Float64,2}, xt::Array{Complex{Float64},2}, Dxt::Array{Complex{Float64,2}, k::Int64, A1h::Array{Float64,2}, a2h::Array{Float64,1}, mu::Float64=1e-1, bp::Float64=1e-1)

Minimizes the quadratic error in the estimation of the forced trajectory, for a fixed frequency (k) of the forcing. The optimization parameters are the dynamics matrix (A1), the damings (a2), and the forcing amplitude (gamma). 

_INPUT_:
`x`: Time series of the phase angles.
`Dx`: Time series of the phase frequencies. 
`xt`: (Inverse) Fourier transform of x.
`Dxt`: (Inverse) Fourier transform of Dx.
`k`: Fixed index of the forcing frequency (ν = k/T).
`A1h`: Warm start for A1.
`a2h`: Warm starrt for a2.
`mu`: Initial value of the barrier parameter (in IPOPT).
`bp`: Initial value of the bound_push parameter (in IPOPT).

_OUTPUT_:
`objective`: Value of the optimized objective function.
`A1`: Best estimate of the dynamcis matrix.
`a2`: Best estimate of the dampings.
`gamma`: Best estimate of the forcing amplitude at each possible location.
"""
function Lmin_l2(x::Array{Float64,2}, Dx::Array{Float64,2}, xt::Array{Complex{Float64},2}, Dxt::Array{Complex{Float64},2}, k::Int64, A1h::Array{Float64,2}, a2h::Array{Float64,1}, mu::Float64=1e-1, bp::Float64=1e-1)
	t0 = time()
	
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
#	for i in 1:n-1
#		for j in i+1:n
#			@constraint(system_id, A1[i,j] == A1[j,i])
#		end
#	end
	set_start_value.(A1,A1h)
	@variable(system_id, a2[i = 1:n])
	@variable(system_id, gamma[i = 1:n])
	for i in 1:n
#		@constraint(system_id, a2[i] >= 0.)
		@constraint(system_id, gamma[i] >= 0.)
	end
	set_start_value.(a2,a2h)
	set_start_value.(gamma,ones(n))

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

	@info "Full optimization took $(time() - t0)''."
	
	return objective_value(system_id), value.(A1), value.(a2), value.(gamma)
end

"""
	objective(Xs::Array{Float64,2}, tau::Float64, A1::Array{Float64,2}, d::Array{Float64,1}, gamma::Float64, l::Int64, k::Int64)

Computes the value of the objective function for the given paramters.

_INPUT_:
`Xs`: Time series of the phase angles (rows 1:n) and of the phase frequencies (rows n+1:2*n).
`tau`: Time step size.
`A1`: Dynamics matrix.
`d`: Dampings.
`gamma`: Forcing amplitude.
`l`: Forcing location.
`k`: Index of the forcing frequency (ν = k/T).

_OUTPUT_:
`obj`: Value of the objective.
"""
function objective(Xs::Array{Float64,2}, tau::Float64, A1::Array{Float64,2}, a2::Array{Float64,1}, gamma::Float64, l::Int64, k::Int64)
	x = Xs[:,1:end-1]
	nn,N = size(x)
	n = Int(nn/2)

	Dx = (Xs[(n+1):(2*n),2:end] - Xs[(n+1):(2*n),1:end-1])/tau
	xt = Array{Complex{Float64},2}(undef,nn,N)
	for i in 1:nn
		xt[i,:] = ifft(x[i,:])*sqrt(N)
	end
	Dxt = Array{Complex{Float64},2}(undef,n,N)
	for i in 1:n
		Dxt[i,:] = ifft(Dx[i,:])*sqrt(N)
	end
	
	Sigma0 = (x*x')/N
	Sigma1 = (x*Dx')/N

	xtk = xt[:,k]
	Dxtlk = Dxt[l,k]
	
	Fk = real.(xtk*xtk')
	flk = real.(Dxtlk*xtk')
	glk = norm(Dxtlk)^2

	A = [A1 diagm(0 => a2)]

	obj = tr(transpose(A)*A*Sigma0) + 2*tr(A*Sigma1) + .5*gamma^2 - 2*gamma/sqrt(N)*sqrt(tr(A[l,:]*transpose(A[l,:])*Fk) + (2*flk*A[l,:])[1] + glk)

	return obj
end

"""
	plot_new_l0(Ls_l0::Array{Float64,2}, L_l0::Float64, l_l0::Int64, k_l0::Int64, ks::Tuple{Int64,Int64,Int64}, T::Union{Float64,Int64}, n::Int64)

Plots the value of the objective vs. the estimated (discretized) forcing frequencies by the l0 approach.

_INPUT_:
`Ls_l0`: Values of the objective for the various values of l (rows) and k (columns). 
`L_l0`: Minimal objective value.
`l_l0`: Corresponding value of l.
`k_l0`: Corresponding value of k.
`ks = (kmin,kmax,dk)`: Values of k assessed (kmin:dk:kmax).
`T`: Observation time.
`n`: Number of agents.
"""
function plot_new_l0(Ls_l0::Array{Float64,2}, L_l0::Float64, l_l0::Int64, k_l0::Int64, ks::Tuple{Int64,Int64,Int64}, T::Union{Float64,Int64}, n::Int64)
	kmin,kmax,dk = ks
	ws = 2*pi*(kmin:dk:kmax)/T
	for l in 1:n
		PyPlot.plot(ws,Ls_l0[l,:],"o",label="l = $l")
	end
	xlabel("ω")
	ylabel("Log-likelihood")
	title("Best l0: l = $(l_l0), ω = $(round(2*pi*k_l0/T,sigdigits=5)), L = $(round(L_l0,sigdigits=5))")
end

"""
	plot_new_l2_1(Ls_l2::Array{Float64,1}, L_l2::Float64, gamma_l2::Array{Float64,1}, l_l2::Int64, k_l0::Int64, ks::Tuple{Int64,Int64,Int64}, T::Union{Float64,Int64})

Plots the value of the objective vs. the estimated (discretized) forcing frequencies by the l2 approach.

_INPUT_:
`Ls_l2`: Values of the objective for the various values of k. 
`L_l2`: Minimal objective value.
`gamma_l2`: Corresponding forcing amplitudes.
`l_l2`: Index with largest amplitude.
`k_l2`: Corresponding value of k.
`ks = (kmin,kmax,dk)`: Values of assessed (kmin:dk:kmax).
"""
function plot_new_l2_1(Ls_l2::Array{Float64,1}, L_l2::Float64, gamma_l2::Array{Float64,1}, l_l2::Int64, k_l2::Int64, ks::Tuple{Int64,Int64,Int64}, T::Union{Float64,Int64})
	kmin,kmax,dk = ks
	ws = 2*pi*(kmin:dk:kmax)/T
	PyPlot.plot(ws,Ls_l2,"o")
	xlabel("ω")
	ylabel("Log-likelihood")
	title("Best l2: ω = $(round(2*pi*k_l2/T,sigdigits=5)), L = $(round(L_l2,sigdigits=5))")
end

"""
	plot_new_l2_2(Ls_l2::Array{Float64,1}, L_l2::Float64, gamma_l2::Array{Float64,1}, l_l2::Int64)

Plots the value of the estimated forcing amplitude vs. the agent index, obtained by the l2 approach.

_INPUT_:
`Ls_l2`: Values of the objective for the various values of k. 
`L_l2`: Minimal objective value.
`gamma_l2`: Corresponding forcing amplitudes.
`l_l2`: Index with largest amplitude.
"""
function plot_new_l2_2(Ls_l2::Array{Float64,1}, L_l2::Float64, gamma_l2::Array{Float64,1}, l_l2::Int64)
	PyPlot.plot(1:length(gamma_l2),gamma_l2,"o")
	xlabel("node index")
	ylabel("Estimated amplitude")
	title("Source: l = $(l_l2), γ_l = $(round(gamma_l2[l_l2],sigdigits=4))")
end


"""
	TODO
"""
function plot_new_error(Lh::Array{Array{Float64,2},1}, L::Array{Float64,2}, dh::Array{Array{Float64,1},1}, d::Array{Float64,1}, gammah::Array{Array{Float64,1},1}, gamma::Array{Float64,1}, wh::Array{Float64,1}, w::Float64, Xss::Array{Array{Float64,2},1}, taus::Array{Float64,1}, ls::Array{Int64,1}, ks::Array{Int64,1})
	n = length(Lh)
	
	bs = Array(LinRange(-.4,.4,n+1))
	xs = (bs[2:end] + bs[1:end-1])/2
	wi = bs[2] - bs[1]
	a = Array(LinRange(1.,.3,n+1))
	
	figure(333)
	subplot(1,3,1)
	for i in 1:n
		PyPlot.bar(1 + xs[i], maximum(abs.(Lh[i] - L)./max.(abs.(L),1e-8).*(abs.(L) .> 1e-8)), wi, color="C0", alpha=a[i])
		PyPlot.bar(2 + xs[i], maximum(abs.(dh[i] - d)./max.(abs.(d),1e-8).*(abs.(d) .> 1e-8)), wi, color="C1", alpha=a[i])
		PyPlot.bar(3 + xs[i], maximum(abs.(gammah[i] - gamma)./max.(abs.(gamma),1e-8).*(abs.(gamma) .> 1e-8)), wi, color="C2", alpha=a[i])
		PyPlot.bar(4 + xs[i], abs(wh[i] - w)/abs(w), wi, color="C3", alpha=a[i])
	end

	xticks([1,2,3,4],["L","d","γ","ω"])
	ylabel("Max. relative error")

	subplot(1,3,2)
	for i in 1:n
		PyPlot.bar(1 + xs[i], norm(Lh[i] - L)/norm(L), wi, color="C0", alpha=a[i])
		PyPlot.bar(2 + xs[i], norm(dh[i] - d)/norm(d), wi, color="C1", alpha=a[i])
		PyPlot.bar(3 + xs[i], norm(gammah[i] - gamma)/norm(gamma), wi, color="C2", alpha=a[i])
		PyPlot.bar(4 + xs[i], abs(wh[i] - w)/abs(w), wi, color="C3", alpha=a[i])
	end

	xticks([1,2,3,4],["L","d","γ","ω"])
	ylabel("Relative error")

	bs = Array(LinRange(-.4,.4,n+2))
	xs = (bs[2:end] + bs[1:end-1])/2
	wi = bs[2] - bs[1]

	o0s = [objective(Xss[i],taus[i],L,d,maximum(gamma),ls[i],ks[i]) for i in 1:n]
	os = [objective(Xss[i],taus[i],Lh[i],dh[i],maximum(gammah[i]),ls[i],ks[i]) for i in 1:n]
	
	subplot(1,3,3)
	PyPlot.plot(1:n,os-o0s,"-o",color="C4")

	xticks(1:n)
	ylabel("Log-likelihood")
end

"""
    get_Ah_correl_new(Xs::Array{Float64,2}, dt::Float64)

Estimates the dynamics matrix bases on Lokhov18.

_INPUT_:
`Xs`: Time series of the phase angles (rows 1:n) and of the phase frequencies (n+1:2*n).
`dt`: Time step.

_OUTPUT_: 
`Ah`: Estimate of the full dynamics matrix.
`Lh`: Estimate of the Laplacian matrix normalized by the inertias (M^{-1}*L).
`dh`: Estimate of the damping over inertia ratios.
"""
function get_Ah_correl_new(Xs::Array{Float64,2}, dt::Float64)
	nn,T = size(Xs)
	n = Int(nn/2)

	S0 = Xs*Xs' ./ T
	S1 = Xs[:,2:T]*Xs[:,1:T-1]' ./ (T-1)

	Ah = S1*inv(S0)
	Lh = -Ah[n+1:2*n,1:n]/dt
	dh = (ones(n) - diag(Ah[n+1:2*n,n+1:2*n]))./dt

	return Ah,Lh,dh
end

