using Distributed, Dates

# Define the number of parallel threads.
n_thr = 3

if nworkers() < n_thr
	addprocs(n_thr - nworkers())
end

# If the directory "data" does not exists, then create it.
if !isdir("data/")
	mkdir("data/")
end

@everywhere include("final.jl")

#=
"""
	run_l0_par(id::String, Xs::Array{Float64,2}, τ::Float64, ls::Array{Int64,1}, ks::Array{Int64,1}, is_laplacian::Bool, b::Float64=0., μ::Float64=1e-1, bp::Float64=1e-1)

Parallelized verions of "run_l0".

The parameter "id" identifies the system to be identified.
"""
=#
@everywhere function run_l0_par(id::String, Xs::Array{Float64,2}, τ::Float64, ls::Array{Int64,1}, ks::Array{Int64,1}, is_laplacian::Bool, b::Float64=0., μ::Float64=1e-1, bp::Float64=1e-1)
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

	args = Array{Tuple{String,Array{Float64,2},Array{Float64,2},Array{Complex{Float64},2},Array{Complex{Float64},2},Int64,Int64,Array{Float64,2},Array{Float64,1},Float64,Float64,Float64},1}()

	for l in ls
		for k in ks
			push!(args,(id,x,Dx,xt,Dxt,l,k,A1h,a2h,b,μ,bp))
		end
	end

	if is_laplacian
		pmap(Lmin_l0_lap_par,args)
	else
		pmap(Lmin_l0_par,args)
	end

	@info "ID $id IS DONE !!!"
end

#=
"""
	run_l1_par(Xs::Array{Float64,2}, τ::Float64, ks::Array{Int64,1}, is_laplacian::Bool, b::Float64=0., μ::Float64=1e-1, bp::Float64=1e-1)

Parallelized version of "run_l1".

The parameter "id" identifies the system to be identified.
"""
=#
@everywhere function run_l1_par(id::String, Xs::Array{Float64,2}, τ::Float64, ks::Array{Int64,1}, is_laplacian::Bool, b::Float64=0., μ::Float64=1e-1, bp::Float64=1e-1)
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
#	XXX,A1h,a2h = get_Ah_correl(Xs,τ) # Performs poorly.
	A1h = zeros(n,n)
	a2h = zeros(n)

# Run the optimizations
	args = Array{Tuple{String,Array{Float64,2},Array{Float64,2},Array{Complex{Float64},2},Array{Complex{Float64},2},Int64,Array{Float64,2},Array{Float64,1},Float64,Float64,Float64},1}()

	for k in ks
		push!(args,(id,x,Dx,xt,Dxt,k,A1h,a2h,b,μ,bp))
	end

	if is_laplacian
		pmap(Lmin_l1_lap_par,args)
	else
		pmap(Lmin_l1_par,args)
	end
end

#=
"""
	Lmin_l0_par(id::String, x::Array{Float64,2}, Dx::Array{Float64,2}, xt::Array{Complex{Float64},2}, Dxt::Array{Complex{Float64,2}, l::Int64, k::Int64, A1h::Array{Float64,2}, a1h::Array{Float64,1}, b::Float64=0., μ::Float64=1e-1, bp::Float64=1e-1)

Parallelized version of "Lmin_l0", i.e., stores the data in files.
"""
=#
@everywhere function Lmin_l0_par(tups::Tuple{String,Array{Float64,2},Array{Float64,2},Array{Complex{Float64},2},Array{Complex{Float64},2},Int64,Int64,Array{Float64,2},Array{Float64,1},Float64,Float64,Float64})
	id,x,Dx,xt,Dxt,l,k,A1h,a2h,b,μ,bp = tups

	@info "===================================================================================="
	@info "Computing "*id*" with ℓ0 optimization: ℓ = $l, k = $k."
	@info "===================================================================================="
	
	t0 = time()
	
	nn,N = size(x)
	n = Int(nn/2)

# Definition of the needed parameters.
	S0 = (x*x')/N
	S1 = (x*Dx')/N

	xtk = xt[:,k+1]
	
	Fk = real.(xtk*xtk')
	flk = real.(Dxt[l,k+1]*xtk')
	glk = norm(Dxt[l,k+1])^2

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

	writedlm("data/"*id*"_l0_$(l).$(k)_obj.csv",objective_value(system_id),',')
#	writedlm("data/"*id*"_l0_$(l).$(k)_A1.csv",mL,',')
#	writedlm("data/"*id*"_l0_$(l).$(k)_a2.csv",value.(a2),',')
#	writedlm("data/"*id*"_l0_$(l).$(k)_g.csv",value(γ),',')
end


#=
"""
	Lmin_l1_par(id::String, x::Array{Float64,2}, Dx::Array{Float64,2}, xt::Array{Complex{Float64},2}, Dxt::Array{Complex{Float64,2}, k::Int64, A1h::Array{Float64,2}, a2h::Array{Float64,1}, b::Float64=0., μ::Float64=1e-1, bp::Float64=1e-1)

Parallelized version of Lmin_l1.
"""
=#
@everywhere function Lmin_l1_par(tups::Tuple{String,Array{Float64,2},Array{Float64,2},Array{Complex{Float64},2},Array{Complex{Float64},2},Int64,Array{Float64,2},Array{Float64,1},Float64,Float64,Float64})
	id,x,Dx,xt,Dxt,k,A1h,a2h,b,μ,bp = tups
	
	@info "===================================================================================="
	@info "Computing "*id*" with ℓ1 optimization: k = $k."
	@info "===================================================================================="
	
	t0 = time()
	
	nn,N = size(x)
	n = Int(nn/2)

# Definition of the needed parameters.
	S0 = (x*x')/N
	S1 = (x*Dx')/N

	xtk = xt[:,k+1]
	
	Fk = real.(xtk*xtk')
	flk = [real.(Dxt[l,k+1]*xtk') for l = 1:n]
	glk = [norm(Dxt[l,k+1])^2 for l = 1:n]

# Definition of the optimization problem. 
	system_id = Model(optimizer_with_attributes(Ipopt.Optimizer, "mu_init" => μ, "bound_push" => bp))

	@variable(system_id, A1[i = 1:n, j = i:n])
	for i in 1:n
		for j in i:n
			set_start_value(A1[i,j],A1h[i,j])
		end
	end
	@variable(system_id, a2[i = 1:n])
	set_start_value.(a2,a2h)
	@variable(system_id, γ[l = 1:n])
	for l in 1:n
		@constraint(system_id, γ[l] >= 0.)
	end
	set_start_value.(γ,ones(n))

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
	@NLexpression(system_id, AAF1[l = 1:n], 
		      sum(A1[l,i]*A1[l,j]*Fk[j,i] for i = l:n for j = l:n) + 
		      2*sum(A1[i,l]*A1[l,j]*Fk[j,i] for i = 1:l-1 for j = l:n) + 
		      sum(A1[i,l]*A1[j,l]*Fk[j,i] for i = 1:l-1 for j = 1:l-1))

# tr((L_{l,:})^T*D_{l,:}*F(k)) = tr((D_{l,:})^T*L_{l,:}*F(k))
	@NLexpression(system_id, AAF2[l = 1:n], 
		      sum(A1[l,i]*a2[l]*Fk[n+l,i] for i = l:n) + 
		      sum(A1[i,l]*a2[l]*Fk[n+l,i] for i = 1:l-1))

# tr((D_{l,:})^T*D_{l,:}*F(k))
	@NLexpression(system_id, AAF3[l = 1:n], a2[l]^2*Fk[n+l,n+l])

# tr((A_{l,:})^T*A_{l,:}*F(k))
	@NLexpression(system_id, T3[l = 1:n], AAF1[l] + 2*AAF2[l] + AAF3[l])

	
# (f_l(k))^T*A_{l,:}
	@NLexpression(system_id, T4[l = 1:n], sum(flk[l][i]*A1[l,i] for i = l:n) + sum(flk[l][i]*A1[i,l] for i = 1:l-1) + flk[l][n+l]*a2[l])


	@NLexpression(system_id, γ2, sum(γ[l]^2 for l = 1:n))

	
	@NLobjective(system_id, Min, T1 + 2*T2 + .5*γ2 - 2/sqrt(N)*sum(γ[l]*sqrt(T3[l] + 2*T4[l] + glk[l]) for l = 1:n) + b*sum(abs(A1[i,j]) for i = 1:n-1 for j = i+1:n))

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
	
	writedlm("data/"*id*"_l1_$(k)_obj.csv",objective_value(system_id),',')
#	writedlm("data/"*id*"_l1_$(k)_A1.csv",mL,',')
#	writedlm("data/"*id*"_l1_$(k)_a2.csv",value.(a2),',')
	writedlm("data/"*id*"_l1_$(k)_g.csv",value.(γ),',')
end

#=
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
=#
@everywhere function get_Ah_correl(Xs::Array{Float64,2}, dt::Float64)
	nn,T = size(Xs)
	n = Int(nn/2)

	Σ0 = Xs*Xs' ./ T
	Σ1 = Xs[:,2:T]*Xs[:,1:T-1]' ./ (T-1)

	Ah = Σ1*inv(Σ0)
	Lh = -Ah[n+1:2*n,1:n]/dt
	dh = (ones(n) - diag(Ah[n+1:2*n,n+1:2*n]))./dt

	return Ah,Lh,dh
end

#=
"""
	Lmin_l0_lap_par(id::String, x::Array{Float64,2}, Dx::Array{Float64,2}, xt::Array{Complex{Float64},2}, Dxt::Array{Complex{Float64,2}, l::Int64, k::Int64, A1h::Array{Float64,2}, a1h::Array{Float64,1}, b::Float64=0., μ::Float64=1e-1, bp::Float64=1e-1)

Same as Lmin_l0_par, but assumes that the dynamics matrix (A1) is a Laplacian, with nonpositive off-diagonal terms.
"""
=#
@everywhere function Lmin_l0_lap_par(tups::Tuple{String,Array{Float64,2},Array{Float64,2},Array{Complex{Float64},2},Array{Complex{Float64},2},Int64,Int64,Array{Float64,2},Array{Float64,1},Float64,Float64,Float64})
	id,x,Dx,xt,Dxt,l,k,A1h,a2h,b,μ,bp = tups
	
	t0 = time()
	
	nn,N = size(x)
	n = Int(nn/2)

# Definition of the needed parameters.
	S0 = (x*x')/N
	S1 = (x*Dx')/N

	xtk = xt[:,k+1]
	
	Fk = real.(xtk*xtk')
	flk = real.(Dxt[l,k+1]*xtk')
	glk = norm(Dxt[l,k+1])^2

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


# tr(-L*θ*(Dω)^T)
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

	writedlm("data/"*id*"_lap0_$(l).$(k)_obj.csv",objective_value(system_id),',')
#	writedlm("data/"*id*"_lap0_$(l).$(k)_A1.csv",mL,',')
#	writedlm("data/"*id*"_lap0_$(l).$(k)_a2.csv",value.(a2),',')
#	writedlm("data/"*id*"_lap0_$(l).$(k)_g.csv",value(γ),',')
end


#=
"""
	Lmin_l1_lap_par(id::String, x::Array{Float64,2}, Dx::Array{Float64,2}, xt::Array{Complex{Float64},2}, Dxt::Array{Complex{Float64,2}, k::Int64, A1h::Array{Float64,2}, a1h::Array{Float64,1}, b::Float64=0., μ::Float64=1e-1, bp::Float64=1e-1)

Same as Lmin_l1_par, but assumes that the dynamics matrix (A1) is a Laplacian, with nonpositive off-diagonal terms.
"""
=#
@everywhere function Lmin_l1_lap_par(tups::Tuple{String,Array{Float64,2},Array{Float64,2},Array{Complex{Float64},2},Array{Complex{Float64},2},Int64,Array{Float64,2},Array{Float64,1},Float64,Float64,Float64})
	id,x,Dx,xt,Dxt,k,A1h,a2h,b,μ,bp = tups
	
	t0 = time()
	
	nn,N = size(x)
	n = Int(nn/2)

# Definition of the needed parameters.
	S0 = (x*x')/N
	S1 = (x*Dx')/N

	xtk = xt[:,k+1]
	
	Fk = real.(xtk*xtk')
	flk = [real.(Dxt[l,k+1]*xtk') for l = 1:n]
	glk = [norm(Dxt[l,k+1])^2 for l = 1:n]

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
	@variable(system_id, γ[l = 1:n])
	for l in 1:n
		@constraint(system_id, γ[l] >= 0.)
	end
	set_start_value.(γ,ones(n))

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


# tr(-L*θ*(Dω)^T)
	@NLexpression(system_id, AS11, 
		      sum(A1[i,j]*(S1[j,i]-S1[i,i]) for i = 1:n for j = i+1:n) + 
		      sum(A1[j,i]*(S1[j,i]-S1[i,i]) for i = 1:n for j = 1:i-1))

# tr(-D*ω*(Dω)^T)
	@NLexpression(system_id, AS12, sum(a2[i]*S1[n+i,i] for i = 1:n))

# tr(A*S1)
	@NLexpression(system_id, T2, AS11 + AS12)


# tr((L_{l,:})^T*L_{l,:}*F(k))
	@NLexpression(system_id, AAF1[l = 1:n], 
		      sum(A1[i,l]*A1[j,l]*(Fk[j,i]-Fk[l,i]-Fk[j,l]+Fk[l,l]) for i = 1:l-1 for j = 1:l-1) + 
		      sum(A1[i,l]*A1[l,j]*(Fk[j,i]-Fk[l,i]-Fk[j,l]+Fk[l,l]) for i = 1:l-1 for j = l+1:n) + 
		      sum(A1[l,i]*A1[j,l]*(Fk[j,i]-Fk[l,i]-Fk[j,l]+Fk[l,l]) for i = l+1:n for j = 1:l-1) + 
		      sum(A1[l,i]*A1[l,j]*(Fk[j,i]-Fk[l,i]-Fk[j,l]+Fk[l,l]) for i = l+1:n for j = l+1:n))

# tr((L_{l,:})^T*D_{l,:}*F(k)) = tr((D_{l,:})^T*L_{l,:}*F(k))
	@NLexpression(system_id, AAF2[l = 1:n], 
		      sum(A1[l,i]*a2[l]*(Fk[n+l,i]-Fk[n+l,l]) for i = l+1:n) + 
		      sum(A1[i,l]*a2[l]*(Fk[n+l,i]-Fk[n+l,l]) for i = 1:l-1))

# tr((D_{l,:})^T*D_{l,:}*F(k))
	@NLexpression(system_id, AAF3[l = 1:n], a2[l]^2*Fk[n+l,n+l])

# tr((A_{l,:})^T*A_{l,:}*F(k))
	@NLexpression(system_id, T3[l = 1:n], AAF1[l] + 2*AAF2[l] + AAF3[l])

	
# (f_l(k))^T*A_{l,:}
	@NLexpression(system_id, T4[l = 1:n], sum((flk[l][i]-flk[l][l])*A1[l,i] for i = l+1:n) + sum((flk[l][i]-flk[l][l])*A1[i,l] for i = 1:l-1) + flk[l][n+l]*a2[l])


	@NLexpression(system_id, γ2, sum(γ[l]^2 for l = 1:n))

	
	@NLobjective(system_id, Min, T1 + 2*T2 + .5*γ2 - 2/sqrt(N)*sum(γ[l]*sqrt(T3[l] + 2*T4[l] + glk[l]) for l = 1:n) + b*sum(abs(A1[i,j]) for i = 1:n-1 for j = i+1:n))

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

	writedlm("data/"*id*"_lap2_$(k)_obj.csv",objective_value(system_id),',')
#	writedlm("data/"*id*"_lap2_$(k)_A1.csv",mL,',')
#	writedlm("data/"*id*"_lap2_$(k)_a2.csv",value.(a2),',')
	writedlm("data/"*id*"_lap2_$(k)_g.csv",value.(γ),',')
end

