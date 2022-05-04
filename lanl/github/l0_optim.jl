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
function run_l0(Xs::Array{Float64,2}, τ::Float64, ls::Array{Int64,1}, ks::Array{Int64,1}, is_laplacian::Bool, plot::Bool=false, b::Float64=0., μ::Float64=1e-1, bp::Float64=1e-1)
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

# Compute warm start
#	XXX,A1h,a2h = get_Ah_correl(Xs,τ) # Performs poorly for Laplcian dynamics, warm start at zero is better.
	A1h = zeros(n,n)
	a2h = ones(n)

# Run the optimizations
	Ls_l0 = zeros(length(ls),length(ks))
	L_l0 = 1000.
	A_l0 = zeros(n,nn)
	d_l0 = zeros(n)
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
			if is_laplacian
				Lt = Lmin_l0_lap(x,Dx,xt,Dxt,l,k,A1h,a2h,b,μ,bp)
			else
				Lt = Lmin_l0(x,Dx,xt,Dxt,l,k,A1h,a2h,b,μ,bp)
			end
			Ls_l0[cl,ck] = Lt[1]
			if Lt[1] < L_l0
				L_l0,A_l0,d_l0,γ_l0 = Lt
				l_l0 = l
				k_l0 = k
			end
		end
	end

	@info "Best l0: l = $(l_l0), ω = $(2*π*(k_l0)/(τ*N)), L = $(L_l0)."	

	if plot
		T = N*τ
		figure("l0, T = $(N*τ)")
		plot_l0(Ls_l0, L_l0, l_l0, k_l0, ks, T, n)
	end

	rm("data/times.csv")

	return Ls_l0, (L_l0,A_l0,d_l0,γ_l0,l_l0,k_l0)
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
function Lmin_l0(x::Array{Float64,2}, Dx::Array{Float64,2}, xt::Array{Complex{Float64},2}, Dxt::Array{Complex{Float64},2}, l::Int64, k::Int64, A1h::Array{Float64,2}, a2h::Array{Float64,1}, b::Float64=0., μ::Float64=1e-1, bp::Float64=1e-1)
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
	@variable(system_id, A1[i = 1:n, j = i:n])
	for i = 1:n-1
		for j = i+1:n
			set_start_value(A1[i,j],A1h[i,j])
		end
	end
	@variable(system_id, a2[i = 1:n])
	set_start_value.(a2,a2h)
	@variable(system_id, γ >= 0.)
	set_start_value(γ,1.)

# tr(L^2*θ*θ^T)
	@NLexpression(system_id, AtAS01, 
		      sum(A1[i,j]*A1[j,m]*S0[m,i] for i = 1:n for j = i:n for m = j:n) + 
		      sum(A1[i,j]*A1[m,j]*S0[m,i] for i = 1:n for j = i:n for m = 1:j-1) + 
		      sum(A1[j,i]*A1[j,m]*S0[m,i] for i = 1:n for j = 1:i-1 for m = j:n) + 
		      sum(A1[j,i]*A1[m,j]*S0[m,i] for i = 1:n for j = 1:i-1 for m = 1:j-1))

# tr(-L*D*ω*θ^T) = tr(-D*L*θ*ω^T)
	@NLexpression(system_id, AtAS02, 
		      sum(A1[i,j]*a2[j]*S0[n+j,i] for i = 1:n for j = i:n) + 
		      sum(A1[j,i]*a2[j]*S0[n+j,i] for i = 1:n for j = 1:i-1))

# tr(D^2*ω*ω^T)
	@NLexpression(system_id, AtAS03, sum(a2[i]^2*S0[n+i,n+i] for i = 1:n))
	
# tr(A^T*A*S0)
	@NLexpression(system_id, T1, AtAS01 + 2*AtAS02 + AtAS03)


# tr(-L*θ*(Dω)^T)
	@NLexpression(system_id, AS11, 
		      sum(A1[i,j]*S1[j,i] for i = 1:n for j = i:n) + 
		      sum(A1[j,i]*S1[j,i] for i =1:n for j = 1:i-1))

# tr(-D*ω*(Dω)^T)
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


	@NLexpression(system_id, γ2, γ^2)

	
	@NLobjective(system_id, Min, T1 + 2*T2 + .5*γ2 - 2*γ/sqrt(N)*sqrt(T3 + 2*T4 + glk) + b*sum(abs(A1[i,j]) for i = 1:n-1 for j = i+1:n))

	optimize!(system_id)

	mL = zeros(n,n)
	for i in 1:n
		for j in i+1:n
			mL[i,j] = value(A1[i,j])
			mL[j,i] = value(A1[i,j])
		end
		mL[i,i] = value(A1[i,i])
	end

	t = time()-t0
	TTMM = readdlm("data/times.csv",',')
	TT = TTMM[1]
	MM = TTMM[2]
	avt = (TT*MM + t)/(MM+1)
	writedlm("data/times.csv",[avt,MM+1],',')

	@info "===================================================================================="
	@info "Full optimization took $(round(t,digits=3))'', avg time is $(round(avt,digits=3))''."
	@info "===================================================================================="

	return objective_value(system_id), mL, value.(a2), value(γ)
end



"""
	plot_l0(Ls_l0::Array{Float64,2}, L_l0::Float64, l_l0::Int64, k_l0::Int64, ks::Array{Int64,1}, T::Union{Float64,Int64}, n::Int64)

Plots the value of the objective vs. the estimated (discretized) forcing frequencies by the l0 approach.

_INPUT_:
`Ls_l0`: Values of the objective for the various values of l (rows) and k (columns). 
`L_l0`: Minimal objective value.
`l_l0`: Corresponding value of l.
`k_l0`: Corresponding value of k.
`ks`: Values of k assessed.
`T`: Observation time.
`n`: Number of agents.
"""
function plot_l0(Ls_l0::Array{Float64,2}, L_l0::Float64, l_l0::Int64, k_l0::Int64, ks::Array{Int64,1}, T::Union{Float64,Int64}, n::Int64)
	ωs = 2*π*ks/T
	for l in 1:n
		PyPlot.plot(ωs,Ls_l0[l,:],"o",label="l = $l")
	end
	xlabel("ω")
	ylabel("Log-likelihood")
	title("Best l0: l = $(l_l0), ω = $(round(2*π*k_l0/T,sigdigits=5)), L = $(round(L_l0,sigdigits=5))")
end


"""
    get_Ah_correl(Xs::Array{Float64,2}, dt::Float64)

Estimates the dynamics matrix based on Lokhov18.

_INPUT_:
`Xs`: Time series of the phase angles (rows 1:n) and of the phase frequencies (n+1:2*n).
`dt`: Time step.

_OUTPUT_: 
`Ah`: Estimate of the full dynamics matrix.
`Lh`: Estimate of the Laplacian matrix normalized by the inertias (M^{-1}*L).
`dh`: Estimate of the damping over inertia ratios.
"""
function get_Ah_correl(Xs::Array{Float64,2}, dt::Float64)
	nn,T = size(Xs)
	n = Int(nn/2)

	Σ0 = Xs*Xs' ./ T
	Σ1 = Xs[:,2:T]*Xs[:,1:T-1]' ./ (T-1)

	Ah = Σ1*inv(Σ0)
	Lh = -Ah[n+1:2*n,1:n]/dt
	dh = (ones(n) - diag(Ah[n+1:2*n,n+1:2*n]))./dt

	return Ah,Lh,dh
end

"""
	Lmin_l0_lap(x::Array{Float64,2}, Dx::Array{Float64,2}, xt::Array{Complex{Float64},2}, Dxt::Array{Complex{Float64,2}, l::Int64, k::Int64, A1h::Array{Float64,2}, a1h::Array{Float64,1}, b::Float64=0., μ::Float64=1e-1, bp::Float64=1e-1)

Same as Lmin_l0, but assumes that the dynamics matrix (A1) is a Laplacian, with nonpositive off-diagonal terms.
"""

function Lmin_l0_lap(x::Array{Float64,2}, Dx::Array{Float64,2}, xt::Array{Complex{Float64},2}, Dxt::Array{Complex{Float64},2}, l::Int64, k::Int64, A1h::Array{Float64,2}, a2h::Array{Float64,1}, b::Float64=0., μ::Float64=1e-1, bp::Float64=1e-1)
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
	@variable(system_id, A1[i = 1:n-1, j = i+1:n])
	for i = 1:n-1
		for j = i+1:n
			@constraint(system_id, A1[i,j] <= 0.)
			set_start_value(A1[i,j],A1h[i,j])
		end
	end
	@variable(system_id, a2[i = 1:n])
	set_start_value.(a2,a2h)
	@variable(system_id, γ >= 0.)
	set_start_value(γ,1.)

# tr(L^2*θ*θ^T)
	@NLexpression(system_id, AtAS01, 
		      sum(A1[j,i]*A1[m,j]*(S0[m,i]-S0[j,i]) for i = 1:n for j = 1:i-1 for m = 1:j-1) + 
		      sum(A1[j,i]*A1[j,m]*(S0[m,i]-S0[j,i]) for i = 1:n for j = 1:i-1 for m = j+1:n) -
		      sum(A1[j,i]*A1[m,i]*(S0[m,i]-S0[i,i]) for i = 1:n for j = 1:i-1 for m = 1:i-1) - 
		      sum(A1[j,i]*A1[i,m]*(S0[m,i]-S0[i,i]) for i = 1:n for j = 1:i-1 for m = i+1:n) +
		      sum(A1[i,j]*A1[m,j]*(S0[m,i]-S0[j,i]) for i = 1:n for j = i+1:n for m = 1:j-1) + 
		      sum(A1[i,j]*A1[j,m]*(S0[m,i]-S0[j,i]) for i = 1:n for j = i+1:n for m = j+1:n) -
		      sum(A1[i,j]*A1[m,i]*(S0[m,i]-S0[i,i]) for i = 1:n for j = i+1:n for m = 1:i-1) -
		      sum(A1[i,j]*A1[i,m]*(S0[m,i]-S0[i,i]) for i = 1:n for j = i+1:n for m = i+1:n))

# tr(-L*D*ω*θ^T) = tr(-D*L*θ*ω^T)
	@NLexpression(system_id, AtAS02, 
		      sum(A1[j,i]*(a2[j]*S0[n+j,i]-a2[i]*S0[n+i,i]) for i = 1:n for j = 1:i-1) + 
		      sum(A1[i,j]*(a2[j]*S0[n+j,i]-a2[i]*S0[n+i,i]) for i = 1:n for j = i+1:n))

# tr(D^2*ω*ω^T)
	@NLexpression(system_id, AtAS03, sum(a2[i]^2*S0[n+i,n+i] for i = 1:n))
	
# tr(A^T*A*S0)
	@NLexpression(system_id, T1, AtAS01 + 2*AtAS02 + AtAS03)


# tr(-L*\th*(Dω)^T)
	@NLexpression(system_id, AS11, 
		      sum(A1[i,j]*(S1[j,i]-S1[i,i]) for i = 1:n for j = i+1:n) + 
		      sum(A1[j,i]*(S1[j,i]-S1[i,i]) for i = 1:n for j = 1:i-1))

# tr(-D*ω*(Dω)^T)
	@NLexpression(system_id, AS12, sum(a2[i]*S1[n+i,i] for i = 1:n))

# tr(A*S1)
	@NLexpression(system_id, T2, AS11 + AS12)


# tr((L_{l,:})^T*L_{l,:}*F(k))
	@NLexpression(system_id, AAF1, 
		      sum(A1[i,l]*A1[j,l]*(Fk[j,i]-Fk[l,i]-Fk[j,l]+Fk[l,l]) for i = 1:l-1 for j = 1:l-1) + 
		      sum(A1[i,l]*A1[l,j]*(Fk[j,i]-Fk[l,i]-Fk[j,l]+Fk[l,l]) for i = 1:l-1 for j = l+1:n) + 
		      sum(A1[l,i]*A1[j,l]*(Fk[j,i]-Fk[l,i]-Fk[j,l]+Fk[l,l]) for i = l+1:n for j = 1:l-1) + 
		      sum(A1[l,i]*A1[l,j]*(Fk[j,i]-Fk[l,i]-Fk[j,l]+Fk[l,l]) for i = l+1:n for j = l+1:n))

# tr((L_{l,:})^T*D_{l,:}*F(k)) = tr((D_{l,:})^T*L_{l,:}*F(k))
	@NLexpression(system_id, AAF2, 
		      sum(A1[l,i]*a2[l]*(Fk[n+l,i]-Fk[n+l,l]) for i = l+1:n) + 
		      sum(A1[i,l]*a2[l]*(Fk[n+l,i]-Fk[n+l,l]) for i = 1:l-1))

# tr((D_{l,:})^T*D_{l,:}*F(k))
	@NLexpression(system_id, AAF3, a2[l]^2*Fk[n+l,n+l])

# tr((A_{l,:})^T*A_{l,:}*F(k))
	@NLexpression(system_id, T3, AAF1 + 2*AAF2 + AAF3)

	
# (f_l(k))^T*A_{l,:}
	@NLexpression(system_id, T4, sum((flk[i]-flk[l])*A1[l,i] for i = l+1:n) + sum((flk[i]-flk[l])*A1[i,l] for i = 1:l-1) + flk[n+l]*a2[l])


	@NLexpression(system_id, γ2, γ^2)

	
	@NLobjective(system_id, Min, T1 + 2*T2 + .5*γ2 - 2*γ/sqrt(N)*sqrt(T3 + 2*T4 + glk) + b*sum(abs(A1[i,j]) for i = 1:n-1 for j = i+1:n))

	optimize!(system_id)

	mL = zeros(n,n)
	for i in 1:n-1
		for j in i+1:n
			mL[i,j] = value(A1[i,j])
			mL[j,i] = value(A1[i,j])
		end
	end
	for i in 1:n
		mL[i,i] = -sum(mL[i,j] for j = 1:n)
	end

	t = time()-t0
	TTMM = readdlm("data/times.csv",',')
	TT = TTMM[1]
	MM = TTMM[2]
	avt = (TT*MM + t)/(MM+1)
	writedlm("data/times.csv",[avt,MM+1],',')

	@info "===================================================================================="
	@info "Full optimization took $(round(t,digits=3))'', avg time is $(round(avt,digits=3))''."
	@info "===================================================================================="

	return objective_value(system_id), mL, value.(a2), value(γ)
end



function hms(T::Union{Int64,Float64})
	h = floor(Int64,T/3600)
	m = floor(Int64,(T - h*3600)/60)
	s = round(Int64,T - h*3600 - m*60)

	return h,m,s
end
