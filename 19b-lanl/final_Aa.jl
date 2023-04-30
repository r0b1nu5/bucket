using DelimitedFiles, PyPlot, LinearAlgebra, JuMP, Ipopt, FFTW

"""
	run_l0(Xs::Array{Float64,2}, τ::Float64, ls::Array{Int64,1}, ks::Array{Int64,1}, is_laplacian::Bool, plot::Bool=false, b::Float64=0., μ::Float64=1e-1, bp::Float64=1e-1)

Identifies dynamics and forcing characteristics baded only on measurements. Approximation of the forcing's frequency is discretized as 2*π*k/T, with k integer.

_INPUT_:
`Xs`: Time series of the phase angles (rows 1:n) and of the phase frequencies (rows n+1:2*n).
`τ`: Time step size.
`ls`: Values of l to be tried.
`ks`: Values of k to be tried.
`is_laplacian`: If true, assumes that the dynamics matrix (A1) is Laplacian with nonpositive off-diagonal terms. 
`plot`: If true, generates the plots of the objective function vs. k.
`b`: Regularization parameter to avoid overfitting.
`μ`: Initial value of the barrier parameter (in IPOPT).
`bp`: Initial value of the bound_push parameter (in IPOPT).

_OUTPUT_:
`Ls_l0`: values of the objective function for the various values of k and l, under l0.
`(L_l0, A_l0, d_l0, γ_l0, k_l0, l_l0)`: results under l0.
	`L_l0`: Minimal value of the objective.
	`A_l0`: Estimate of the dynamics matrix.
	`d_l0`: Estimate of the dampings.
	`γ_l0`: Estimate of the forcing amplitude.
	`k_l0`: Estimate frequency index (see theory).
	`l_l0`: Estimate of the forcing location.
"""
function run_l0_Aa(Xs::Array{Float64,2}, τ::Float64, A::Matrix{Float64}, a::Vector{Float64}, ls::Array{Int64,1}, ks::Array{Int64,1}, plot::Bool=false, b::Float64=0., μ::Float64=1e-1, bp::Float64=1e-1)
	writedlm("data/times.csv",[0.,0.],',')

	nn,NN = size(Xs)
	n = Int(nn/2)
	N = NN-1

# Computing the needed inputs (time series, discrete derivative, and their Fourier transforms).
	x = Xs[:,1:end-1]
	Dx = (Xs[(n+1):(2*n),2:end] - Xs[(n+1):(2*n),1:end-1])/τ
	xt = Array{Complex{Float64},2}(undef,nn,N)
	for i in 1:nn
		xt[i,:] = ifft(x[i,:])*sqrt(N)
	end
	Dxt = Array{Complex{Float64},2}(undef,n,N)
	for i in 1:n
		Dxt[i,:] = ifft(Dx[i,:])*sqrt(N)
	end

# Run the optimizations
	Ls_l0 = zeros(length(ls),length(ks))
	L_l0 = 1000.
	γ_l0 = 0.
	l_l0 = 0
	k_l0 = 0
	cl = 0
	for l in ls
		cl += 1
		ck = 0
		for k in ks
			ck += 1
			@info "l0: l = $l, k = $k"
			Lt = Lmin_l0_Aa(x,Dx,xt,Dxt,l,k,A,a,b,μ,bp)
			Ls_l0[cl,ck] = Lt[1]
			if Lt[1] < L_l0
				L_l0,γ_l0 = Lt
				l_l0 = l
				k_l0 = k
			end
		end
	end

	@info "Best l0: l = $(l_l0), ω = $(2*π*(k_l0)/(τ*N)), L = $(L_l0)."	

	if plot
		T = N*τ
		figure("l0-Aa, T = $(N*τ)")
		plot_l0(Ls_l0, L_l0, l_l0, k_l0, ks, T, n)
	end

	rm("data/times.csv")

	return Ls_l0, (L_l0,A,a,γ_l0,l_l0,k_l0)
end


"""
	run_l1(Xs::Array{Float64,2}, τ::Float64, ks::Array{Int64,1}, is_laplacian::Bool, plot::Bool=false, b::Float64=0., μ::Float64=1e-1, bp::Float64=1e-1)

Identifies dynamics and forcing characteristics baded only on measurements. Approximation of the forcing's frequency is discretized as 2*π*k/T, with k integer.

_INPUT_:
`Xs`: Time series of the phase angles (rows 1:n) and of the phase frequencies (rows n+1:2*n).
`τ`: Time step size.
`ks`: Values of k to try.
`is_laplacian`: If true, assumes that the dynamics matrix (A1) is Laplacian with nonpositive off-diagonal terms. 
`plot`: If true, generates the plots of the objective function vs. k.
`b`: Regularization parameters to avoid overfitting.
`μ`: Initial value of the barrier parameter (in IPOPT).
`bp`: Initial value of the bound_push parameter (in IPOPT).

_OUTPUT_:
`Ls_l1`: values of the objective function for the various values of k and l, under l1.
`(L_l1, A_l1, d_l1, γ_l1, k_l1, l_l1)`: results under l1.
	`L_l1`: Minimal value of the objective.
	`A_l1`: Estimate of the dynamics matrix.
	`d_l1`: Estimate of the damings.
	`γ_l1`: Estimate of the forcing amplitude.
	`k_l1`: Estimate frequency index (see theory).
	`l_l1`: Estimate of the forcing location.
"""
function run_l1_Aa(Xs::Array{Float64,2}, τ::Float64, A::Matrix{Float64}, a::Vector{Float64}, ks::Array{Int64,1}, plot::Bool=false, b::Tuple{Float64,Float64}=(0.,0.), μ::Float64=1e-1, bp::Float64=1e-1)
	writedlm("data/times.csv",[0.,0.],',')

	nn,NN = size(Xs)
	n = Int(nn/2)
	N = NN-1

# Computing the needed inputs (time series, discrete derivative, and their Fourier transforms).
	x = Xs[:,1:end-1]
	Dx = (Xs[(n+1):(2*n),2:end] - Xs[(n+1):(2*n),1:end-1])/τ
	xt = Array{Complex{Float64},2}(undef,nn,N)
	for i in 1:nn
		xt[i,:] = ifft(x[i,:])*sqrt(N)
	end
	Dxt = Array{Complex{Float64},2}(undef,n,N)
	for i in 1:n
		Dxt[i,:] = ifft(Dx[i,:])*sqrt(N)
	end

# Run the optimizations
	Ls_l1 = zeros(length(ks))
	L_l1 = 1000.
	γ_l1 = zeros(n)
	k_l1 = 0
	l_l1 = 0
	c = 0
	for k in ks
		c += 1
		@info "l1: k = $k"
		Lt = Lmin_l1_Aa(x,Dx,xt,Dxt,k,A,a,b,μ,bp)
		Ls_l1[c] = Lt[1]
		if Lt[1] < L_l1
			L_l1,γ_l1 = Lt
			k_l1 = k
			l_l1 = findmax(γ_l1)[2]
		end
	end

	@info "Best l1: ω = $(2*π*(k_l1)/(τ*N)), L = $(L_l1)."

	if plot
		T = N*τ
		figure("l1-Aa, T = $(N*τ)")
		subplot(1,2,1)
		plot_l1_1(Ls_l1, L_l1, γ_l1, l_l1, k_l1, ks, T)
		subplot(1,2,2)
		plot_l1_2(Ls_l1, L_l1, γ_l1, l_l1)
	end

	rm("data/times.csv")

	return Ls_l1, (L_l1,A,a,γ_l1,l_l1,k_l1)
end


"""
	Lmin_l0(x::Array{Float64,2}, Dx::Array{Float64,2}, xt::Array{Complex{Float64},2}, Dxt::Array{Complex{Float64,2}, l::Int64, k::Int64, A1h::Array{Float64,2}, a1h::Array{Float64,1}, b::Float64=0., μ::Float64=1e-1, bp::Float64=1e-1)

Minimizes the quadratic error in the estimation of the forced trajectory, for a fixed frequency (k) and location (l) of the forcing. The optimization parameters are the dynamics matrix (A1), the damings (a2), and the forcing amplitude (γ). 

_INPUT_:
`x`: Time series of the phase angles.
`Dx`: Time series of the phase frequencies. 
`xt`: (Inverse) Fourier transform of x.
`Dxt`: (Inverse) Fourier transform of Dx.
`l`: Fixed location of the forcing.
`k`: Fixed index of the forcing frequency (ν = k/T).
`A1h`: Warm start for A1.
`a2h`: Warm start for a2.
`b`: Regularization parameter to avoid overfitting.
`μ`: Initial value of the barrier parameter (in IPOPT).
`bp`: Initial value of the bound_push parameter (in IPOPT).

_OUTPUT_:
`objective`: Value of the optimized objective function.
`A1`: Best estimate of the dynamcis matrix.
`a2`: Best estimate of the dampings.
`γ`: Best estimate of the forcing amplitude.
"""
function Lmin_l0_Aa(x::Array{Float64,2}, Dx::Array{Float64,2}, xt::Array{Complex{Float64},2}, Dxt::Array{Complex{Float64},2}, l::Int64, k::Int64, A::Array{Float64,2}, a::Array{Float64,1}, b::Float64=0., μ::Float64=1e-1, bp::Float64=1e-1)
	t0 = time()
	
	nn,N = size(x)
	n = Int(nn/2)

# Definition of the needed parameters.
	S0 = (x*x')/N
	S1 = (x*Dx')/N

	xtk = xt[:,k+1]		# THE FIRST COLUMN IS f=0, WHICH WE DON'T WANT TO TREAT.
	
	Fk = real.(xtk*xtk')
	flk = real.(Dxt[l,k+1]*xtk') # THE FIRST COLUMN IS f=0, WHICH WE DON'T WANT TO TREAT.

	glk = norm(Dxt[l,k+1])^2 # THE FIRST COLUMN IS f=0, WHICH WE DON'T WANT TO TREAT.


# Definition of the optimization problem.
	system_id = Model(optimizer_with_attributes(Ipopt.Optimizer, "mu_init" => μ, "bound_push" => bp))
	@variable(system_id, γ >= 0.)
	set_start_value(γ,1.)

# tr(L^2*θ*θ^T)
	@NLexpression(system_id, AtAS01, 
		      sum(A[i,j]*A[j,m]*S0[m,i] for i = 1:n for j = i:n for m = j:n) + 
		      sum(A[i,j]*A[m,j]*S0[m,i] for i = 1:n for j = i:n for m = 1:j-1) + 
		      sum(A[j,i]*A[j,m]*S0[m,i] for i = 1:n for j = 1:i-1 for m = j:n) + 
		      sum(A[j,i]*A[m,j]*S0[m,i] for i = 1:n for j = 1:i-1 for m = 1:j-1))

# tr(-L*D*ω*θ^T) = tr(-D*L*θ*ω^T)
	@NLexpression(system_id, AtAS02, 
		      sum(A[i,j]*a[j]*S0[n+j,i] for i = 1:n for j = i:n) + 
		      sum(A[j,i]*a[j]*S0[n+j,i] for i = 1:n for j = 1:i-1))

# tr(D^2*ω*ω^T)
	@NLexpression(system_id, AtAS03, sum(a[i]^2*S0[n+i,n+i] for i = 1:n))
	
# tr(A^T*A*S0)
	@NLexpression(system_id, T1, AtAS01 + 2*AtAS02 + AtAS03)


# tr(-L*θ*(Dω)^T)
	@NLexpression(system_id, AS11, 
		      sum(A[i,j]*S1[j,i] for i = 1:n for j = i:n) + 
		      sum(A[j,i]*S1[j,i] for i =1:n for j = 1:i-1))

# tr(-D*ω*(Dω)^T)
	@NLexpression(system_id, AS12, sum(a[i]*S1[n+i,i] for i = 1:n))

# tr(A*S1)
	@NLexpression(system_id, T2, AS11 + AS12)


# tr((L_{l,:})^T*L_{l,:}*F(k))
	@NLexpression(system_id, AAF1, 
		      sum(A[l,i]*A[l,j]*Fk[j,i] for i = l:n for j = l:n) + 
		      2*sum(A[i,l]*A[l,j]*Fk[j,i] for i = 1:l-1 for j = l:n) + 
		      sum(A[i,l]*A[j,l]*Fk[j,i] for i = 1:l-1 for j = 1:l-1))

# tr((L_{l,:})^T*D_{l,:}*F(k)) = tr((D_{l,:})^T*L_{l,:}*F(k))
	@NLexpression(system_id, AAF2, 
		      sum(A[l,i]*a[l]*Fk[n+l,i] for i = l:n) + 
		      sum(A[i,l]*a[l]*Fk[n+l,i] for i = 1:l-1))

# tr((D_{l,:})^T*D_{l,:}*F(k))
	@NLexpression(system_id, AAF3, a[l]^2*Fk[n+l,n+l])

# tr((A_{l,:})^T*A_{l,:}*F(k))
	@NLexpression(system_id, T3, AAF1 + 2*AAF2 + AAF3)

	
# (f_l(k))^T*A_{l,:}
	@NLexpression(system_id, T4, sum(flk[i]*A[l,i] for i = l:n) + sum(flk[i]*A[i,l] for i = 1:l-1) + flk[n+l]*a[l])


	@NLexpression(system_id, γ2, γ^2)

	
	@NLobjective(system_id, Min, T1 + 2*T2 + .5*γ2 - 2*γ/sqrt(N)*sqrt(T3 + 2*T4 + glk))

	optimize!(system_id)

	t = time()-t0
	TTMM = readdlm("data/times.csv",',')
	TT = TTMM[1]
	MM = TTMM[2]
	avt = (TT*MM + t)/(MM+1)
	writedlm("data/times.csv",[avt,MM+1],',')

	@info "===================================================================================="
	@info "Full optimization took $(round(t,digits=3))'', avg time is $(round(avt,digits=3))''."
	@info "===================================================================================="

	return objective_value(system_id), value(γ)
end



"""
	Lmin_l1(x::Array{Float64,2}, Dx::Array{Float64,2}, xt::Array{Complex{Float64},2}, Dxt::Array{Complex{Float64,2}, k::Int64, A1h::Array{Float64,2}, a2h::Array{Float64,1}, b::Float64=0., μ::Float64=1e-1, bp::Float64=1e-1)

Minimizes the quadratic error in the estimation of the forced trajectory, for a fixed frequency (k) of the forcing. The optimization parameters are the dynamics matrix (A1), the damings (a2), and the forcing amplitude (γ). 

_INPUT_:
`x`: Time series of the phase angles.
`Dx`: Time series of the phase frequencies. 
`xt`: (Inverse) Fourier transform of x.
`Dxt`: (Inverse) Fourier transform of Dx.
`k`: Fixed index of the forcing frequency (ν = k/T).
`A1h`: Warm start for A1.
`a2h`: Warm starrt for a2.
`b`: Regularizaion parameter to avoid overfitting. 
`μ`: Initial value of the barrier parameter (in IPOPT).
`bp`: Initial value of the bound_push parameter (in IPOPT).

_OUTPUT_:
`objective`: Value of the optimized objective function.
`A1`: Best estimate of the dynamcis matrix.
`a2`: Best estimate of the dampings.
`γ`: Best estimate of the forcing amplitude at each possible location.
"""
function Lmin_l1_Aa(x::Array{Float64,2}, Dx::Array{Float64,2}, xt::Array{Complex{Float64},2}, Dxt::Array{Complex{Float64},2}, k::Int64, A::Array{Float64,2}, a::Array{Float64,1}, b::Tuple{Float64,Float64}=(0.,0.), μ::Float64=1e-1, bp::Float64=1e-1)
	t0 = time()
	
	nn,N = size(x)
	n = Int(nn/2)
	b1,b2 = b

# Definition of the needed parameters.
	S0 = (x*x')/N
	S1 = (x*Dx')/N

	xtk = xt[:,k+1] # THE FIRST COLUMN (of xt) CORRESPONDS TO f=0, WHICH WE DON'T WANT TO TREAT.
	
	Fk = real.(xtk*xtk')
	flk = [real.(Dxt[l,k+1]*xtk') for l = 1:n] # THE FIRST COLUMN (of Dxt) CORRESPONDS TO f=0, WHICH WE DON'T WANT TO TREAT.

	glk = [norm(Dxt[l,k+1])^2 for l = 1:n]# THE FIRST COLUMN (of Dxt) CORRESPONDS TO f=0, WHICH WE DON'T WANT TO TREAT.


# Definition of the optimization problem. 
	system_id = Model(optimizer_with_attributes(Ipopt.Optimizer, "mu_init" => μ, "bound_push" => bp))

	@variable(system_id, γ[l = 1:n])
	for l in 1:n
		@constraint(system_id, γ[l] >= 0.)
	end
	set_start_value.(γ,ones(n))

# tr(L^2*θ*θ^T)
	@NLexpression(system_id, AtAS01, 
		      sum(A[i,j]*A[j,m]*S0[m,i] for i = 1:n for j = i:n for m = j:n) + 
		      sum(A[i,j]*A[m,j]*S0[m,i] for i = 1:n for j = i:n for m = 1:j-1) + 
		      sum(A[j,i]*A[j,m]*S0[m,i] for i = 1:n for j = 1:i-1 for m = j:n) + 
		      sum(A[j,i]*A[m,j]*S0[m,i] for i = 1:n for j = 1:i-1 for m = 1:j-1))

# tr(-L*D*ω*θ^T) = tr(-D*L*θ*ω^T)
	@NLexpression(system_id, AtAS02, 
		      sum(A[i,j]*a[j]*S0[n+j,i] for i = 1:n for j = i:n) + 
		      sum(A[j,i]*a[j]*S0[n+j,i] for i = 1:n for j = 1:i-1))

# tr(D^2*ω*ω^T)
	@NLexpression(system_id, AtAS03, sum(a[i]^2*S0[n+i,n+i] for i = 1:n))
	
# tr(A^T*A*S0)
	@NLexpression(system_id, T1, AtAS01 + 2*AtAS02 + AtAS03)


# tr(-L*\th*(Dω)^T)
	@NLexpression(system_id, AS11, 
		      sum(A[i,j]*S1[j,i] for i = 1:n for j = i:n) + 
		      sum(A[j,i]*S1[j,i] for i =1:n for j = 1:i-1))

# tr(-D*ω*(Dω)^T)
	@NLexpression(system_id, AS12, sum(a[i]*S1[n+i,i] for i = 1:n))

# tr(A*S1)
	@NLexpression(system_id, T2, AS11 + AS12)


# tr((L_{l,:})^T*L_{l,:}*F(k))
	@NLexpression(system_id, AAF1[l = 1:n], 
		      sum(A[l,i]*A[l,j]*Fk[j,i] for i = l:n for j = l:n) + 
		      2*sum(A[i,l]*A[l,j]*Fk[j,i] for i = 1:l-1 for j = l:n) + 
		      sum(A[i,l]*A[j,l]*Fk[j,i] for i = 1:l-1 for j = 1:l-1))

# tr((L_{l,:})^T*D_{l,:}*F(k)) = tr((D_{l,:})^T*L_{l,:}*F(k))
	@NLexpression(system_id, AAF2[l = 1:n], 
		      sum(A[l,i]*a[l]*Fk[n+l,i] for i = l:n) + 
		      sum(A[i,l]*a[l]*Fk[n+l,i] for i = 1:l-1))

# tr((D_{l,:})^T*D_{l,:}*F(k))
	@NLexpression(system_id, AAF3[l = 1:n], a[l]^2*Fk[n+l,n+l])

# tr((A_{l,:})^T*A_{l,:}*F(k))
	@NLexpression(system_id, T3[l = 1:n], AAF1[l] + 2*AAF2[l] + AAF3[l])

	
# (f_l(k))^T*A_{l,:}
	@NLexpression(system_id, T4[l = 1:n], sum(flk[l][i]*A[l,i] for i = l:n) + sum(flk[l][i]*A[i,l] for i = 1:l-1) + flk[l][n+l]*a[l])


	@NLexpression(system_id, γ2, sum(γ[l]^2 for l = 1:n))

	
	@NLobjective(system_id, Min, T1 + 2*T2 + .5*γ2 - 2/sqrt(N)*sum(γ[l]*sqrt(T3[l] + 2*T4[l] + glk[l]) for l = 1:n) + b2*sum(γ[l] for l in 1:n))

	optimize!(system_id)

	t = time()-t0
	TTMM = readdlm("data/times.csv",',')
	TT = TTMM[1]
	MM = TTMM[2]
	avt = (TT*MM + t)/(MM+1)
	writedlm("data/times.csv",[avt,MM+1],',')

	@info "===================================================================================="
	@info "Full optimization took $(round(t,digits=3))'', avg time is $(round(avt,digits=3))''."
	@info "===================================================================================="
	
	return objective_value(system_id), value.(γ)
end





