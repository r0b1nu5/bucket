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
	run_l0_par(id::String, Xs::Array{Float64,2}, tau::Float64, ls::Tuple{Int64,Int64,Int64}, ks::Tuple{Int64,Int64,Int64}, is_laplacian::Bool, b::Float64=0., mu::Float64=1e-1, bp::Float64=1e-1)

Parallelized verions of "run_l0".

The parameter "id" identifies the system to be identified.
"""
=#
@everywhere function run_l0_par(id::String, Xs::Array{Float64,2}, tau::Float64, ls::Tuple{Int64,Int64,Int64}, ks::Tuple{Int64,Int64,Int64}, is_laplacian::Bool, b::Float64=0., mu::Float64=1e-1, bp::Float64=1e-1)
	nn,NN = size(Xs)
	n = Int(nn/2)
	N = NN-1

	lmin,lmax,dl = ls
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
#	XXX,A1h,a2h = get_Ah_correl(Xs,tau) # Performs poorly for Laplcian dynamics, warm start at zero is better.
	A1h = zeros(n,n)
	a2h = ones(n)

	args = Array{Tuple{String,Array{Float64,2},Array{Float64,2},Array{Complex{Float64},2},Array{Complex{Float64},2},Int64,Int64,Array{Float64,2},Array{Float64,1},Float64,Float64,Float64},1}()

	for l in lmin:dl:lmax
		for k in kmin:dk:kmax
			push!(args,(id,x,Dx,xt,Dxt,l,k,A1h,a2h,b,mu,bp))
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
	run_l2_par(Xs::Array{Float64,2}, tau::Float64, ks::Tuple{Int64,Int64,Int64}, is_laplacian::Bool, b::Float64=0., mu::Float64=1e-1, bp::Float64=1e-1)

Parallelized version of "run_l2".

The parameter "id" identifies the system to be identified.
"""
=#
@everywhere function run_l2_par(id::String, Xs::Array{Float64,2}, tau::Float64, ks::Tuple{Int64,Int64,Int64}, is_laplacian::Bool, b::Float64=0., mu::Float64=1e-1, bp::Float64=1e-1)
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
#	XXX,A1h,a2h = get_Ah_correl(Xs,tau) # Performs poorly.
	A1h = zeros(n,n)
	a2h = zeros(n)

# Run the optimizations
	args = Array{Tuple{String,Array{Float64,2},Array{Float64,2},Array{Complex{Float64},2},Array{Complex{Float64},2},Int64,Array{Float64,2},Array{Float64,2},Float64,Float64,Float64},1}()

	for k in kmin:dk:kmax
		push!(args,(id,x,Dx,xt,Dxt,k,A1h,a2h,b,mu,bp))
	end

	if is_laplacian
		pmap(Lmin_l2_lap_par,args)
	else
		pmap(Lmin_l2_par,args)
	end
end

#=
"""
	Lmin_l0_par(id::String, x::Array{Float64,2}, Dx::Array{Float64,2}, xt::Array{Complex{Float64},2}, Dxt::Array{Complex{Float64,2}, l::Int64, k::Int64, A1h::Array{Float64,2}, a1h::Array{Float64,1}, b::Float64=0., mu::Float64=1e-1, bp::Float64=1e-1)

Parallelized version of "Lmin_l0", i.e., stores the data in files.
"""
=#
@everywhere function Lmin_l0_par(tups::Tuple{String,Array{Float64,2},Array{Float64,2},Array{Complex{Float64},2},Array{Complex{Float64},2},Int64,Int64,Array{Float64,2},Array{Float64,1},Float64,Float64,Float64})
	id,x,Dx,xt,Dxt,l,k,A1h,a2h,b,mu,bp = tups
	
	t0 = time()
	
	nn,N = size(x)
	n = Int(nn/2)

# Definition of the needed parameters.
	S0 = (x*x')/N
	S1 = (x*Dx')/N

	xtk = xt[:,k+1]
	Dxtlk = Dxt[l,k+1]
	
	Fk = real.(xtk*xtk')
	flk = real.(Dxtlk*xtk')
	glk = norm(Dxtlk)^2

# Definition of the optimization problem.
	system_id = Model(optimizer_with_attributes(Ipopt.Optimizer, "mu_init" => mu, "bound_push" => bp))
	@variable(system_id, A1[i = 1:n, j = i:n])
	for i = 1:n-1
		for j = i+1:n
			set_start_value(A1[i,j],A1h[i,j])
		end
	end
	@variable(system_id, a2[i = 1:n])
	set_start_value.(a2,a2h)
	@variable(system_id, gamma >= 0.)
	set_start_value(gamma,1.)

# tr(L^2*\th*\th^T)
	@NLexpression(system_id, AtAS01, 
		      sum(A1[i,j]*A1[j,m]*S0[m,i] for i = 1:n for j = i:n for m = j:n) + 
		      sum(A1[i,j]*A1[m,j]*S0[m,i] for i = 1:n for j = i:n for m = 1:j-1) + 
		      sum(A1[j,i]*A1[j,m]*S0[m,i] for i = 1:n for j = 1:i-1 for m = j:n) + 
		      sum(A1[j,i]*A1[m,j]*S0[m,i] for i = 1:n for j = 1:i-1 for m = 1:j-1))

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

	
	@NLobjective(system_id, Min, T1 + 2*T2 + .5*g2 - 2*gamma/sqrt(N)*sqrt(T3 + 2*T4 + glk) + b*sum(abs(A1[i,j]) for i = 1:n-1 for j = i+1:n))

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

	writedlm("data/"*id*"_l0_$(l).$(k)_obj.csv",objective_value(system_id),',')
#	writedlm("data/"*id*"_l0_$(l).$(k)_A1.csv",mL,',')
#	writedlm("data/"*id*"_l0_$(l).$(k)_a2.csv",value.(a2),',')
#	writedlm("data/"*id*"_l0_$(l).$(k)_gamma.csv",value(gamma),',')
end

#=
"""
	Lmin_l0_lap_par(id::String, x::Array{Float64,2}, Dx::Array{Float64,2}, xt::Array{Complex{Float64},2}, Dxt::Array{Complex{Float64,2}, l::Int64, k::Int64, A1h::Array{Float64,2}, a1h::Array{Float64,1}, b::Float64=0., mu::Float64=1e-1, bp::Float64=1e-1)

Same as Lmin_l0_par, but assumes that the dynamics matrix (A1) is a Laplacian, with nonpositive off-diagonal terms.
"""
=#
@everywhere function Lmin_l0_lap_par(tups::Tuple{String,Array{Float64,2},Array{Float64,2},Array{Complex{Float64},2},Array{Complex{Float64},2},Int64,Int64,Array{Float64,2},Array{Float64,1},Float64,Float64,Float64})
	id,x,Dx,xt,Dxt,l,k,A1h,a2h,b,mu,bp = tups
	
	t0 = time()
	
	nn,N = size(x)
	n = Int(nn/2)

# Definition of the needed parameters.
	S0 = (x*x')/N
	S1 = (x*Dx')/N

	xtk = xt[:,k+1]
	Dxtlk = Dxt[l,k+1]
	
	Fk = real.(xtk*xtk')
	flk = real.(Dxtlk*xtk')
	glk = norm(Dxtlk)^2

# Definition of the optimization problem.
	system_id = Model(optimizer_with_attributes(Ipopt.Optimizer, "mu_init" => mu, "bound_push" => bp))
	@variable(system_id, A1[i = 1:n-1, j = i+1:n])
	for i = 1:n-1
		for j = i+1:n
			@constraint(system_id, A1[i,j] <= 0.)
			set_start_value(A1[i,j],A1h[i,j])
		end
	end
	@variable(system_id, a2[i = 1:n])
	set_start_value.(a2,a2h)
	@variable(system_id, gamma >= 0.)
	set_start_value(gamma,1.)

# tr(L^2*\th*\th^T)
	@NLexpression(system_id, AtAS01, 
		      sum(A1[j,i]*A1[m,j]*(S0[m,i]-S0[j,i]) for i = 1:n for j = 1:i-1 for m = 1:j-1) + 
		      sum(A1[j,i]*A1[j,m]*(S0[m,i]-S0[j,i]) for i = 1:n for j = 1:i-1 for m = j+1:n) -
		      sum(A1[j,i]*A1[m,i]*(S0[m,i]-S0[i,i]) for i = 1:n for j = 1:i-1 for m = 1:i-1) - 
		      sum(A1[j,i]*A1[i,m]*(S0[m,i]-S0[i,i]) for i = 1:n for j = 1:i-1 for m = i+1:n) +
		      sum(A1[i,j]*A1[m,j]*(S0[m,i]-S0[j,i]) for i = 1:n for j = i+1:n for m = 1:j-1) + 
		      sum(A1[i,j]*A1[j,m]*(S0[m,i]-S0[j,i]) for i = 1:n for j = i+1:n for m = j+1:n) -
		      sum(A1[i,j]*A1[m,i]*(S0[m,i]-S0[i,i]) for i = 1:n for j = i+1:n for m = 1:i-1) -
		      sum(A1[i,j]*A1[i,m]*(S0[m,i]-S0[i,i]) for i = 1:n for j = i+1:n for m = i+1:n))

# tr(-L*D*\om*\th^T) = tr(-D*L*\th*\om^T)
	@NLexpression(system_id, AtAS02, 
		      sum(A1[j,i]*(a2[j]*S0[n+j,i]-a2[i]*S0[n+i,i]) for i = 1:n for j = 1:i-1) + 
		      sum(A1[i,j]*(a2[j]*S0[n+j,i]-a2[i]*S0[n+i,i]) for i = 1:n for j = i+1:n))

# tr(D^2*\om*\om^T)
	@NLexpression(system_id, AtAS03, sum(a2[i]^2*S0[n+i,n+i] for i = 1:n))
	
# tr(A^T*A*S0)
	@NLexpression(system_id, T1, AtAS01 + 2*AtAS02 + AtAS03)


# tr(-L*\th*(D\om)^T)
	@NLexpression(system_id, AS11, 
		      sum(A1[i,j]*(S1[j,i]-S1[i,i]) for i = 1:n for j = i+1:n) + 
		      sum(A1[j,i]*(S1[j,i]-S1[i,i]) for i = 1:n for j = 1:i-1))

# tr(-D*\om*(D\om)^T)
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


	@NLexpression(system_id, g2, gamma^2)

	
	@NLobjective(system_id, Min, T1 + 2*T2 + .5*g2 - 2*gamma/sqrt(N)*sqrt(T3 + 2*T4 + glk) + b*sum(abs(A1[i,j]) for i = 1:n-1 for j = i+1:n))

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

	@info "Full optimization took $(time() - t0)''."

	writedlm("data/"*id*"_lap0_$(l).$(k)_obj.csv",objective_value(system_id),',')
#	writedlm("data/"*id*"_lap0_$(l).$(k)_A1.csv",mL,',')
#	writedlm("data/"*id*"_lap0_$(l).$(k)_a2.csv",value.(a2),',')
#	writedlm("data/"*id*"_lap0_$(l).$(k)_gamma.csv",value(gamma),',')
end

#=
"""
	Lmin_l2_par(id::String, x::Array{Float64,2}, Dx::Array{Float64,2}, xt::Array{Complex{Float64},2}, Dxt::Array{Complex{Float64,2}, k::Int64, A1h::Array{Float64,2}, a2h::Array{Float64,1}, b::Float64=0., mu::Float64=1e-1, bp::Float64=1e-1)

Parallelized version of Lmin_l2.
"""
=#
@everywhere function Lmin_l2_par(tups::Tuple{String,Array{Float64,2},Array{Float64,2},Array{Complex{Float64},2},Array{Complex{Float64},2},Int64,Array{Float64,2},Array{Float64,1},Float64,Float64,Float64})
	id,x,Dx,xt,Dxt,k,A1h,a2h,b,mu,bp = tups
	
	t0 = time()
	
	nn,N = size(x)
	n = Int(nn/2)

# Definition of the needed parameters.
	S0 = (x*x')/N
	S1 = (x*Dx')/N

	xtk = xt[:,k+1]
	
	Fk = real.(xtk*xtk')
	flk = [real.(Dxt[l,k]*xtk') for l = 1:n]
	glk = [norm(Dxt[l,k])^2 for l = 1:n]

# Definition of the optimization problem. 
	system_id = Model(optimizer_with_attributes(Ipopt.Optimizer, "mu_init" => mu, "bound_push" => bp))

	@variable(system_id, A1[i = 1:n, j = i:n])
	for i in 1:n
		for j in i:n
			set_start_value(A1[i,j],A1h[i,j])
		end
	end
	@variable(system_id, a2[i = 1:n])
	set_start_value.(a2,a2h)
	@variable(system_id, gamma[l = 1:n])
	for l in 1:n
		@constraint(system_id, gamma[l] >= 0.)
	end
	set_start_value.(gamma,ones(n))

# tr(L^2*\th*\th^T)
	@NLexpression(system_id, AtAS01, 
		      sum(A1[i,j]*A1[j,m]*S0[m,i] for i = 1:n for j = i:n for m = j:n) + 
		      sum(A1[i,j]*A1[m,j]*S0[m,i] for i = 1:n for j = i:n for m = 1:j-1) + 
		      sum(A1[j,i]*A1[j,m]*S0[m,i] for i = 1:n for j = 1:i-1 for m = j:n) + 
		      sum(A1[j,i]*A1[m,j]*S0[m,i] for i = 1:n for j = 1:i-1 for m = 1:j-1))

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


	@NLexpression(system_id, g2, sum(gamma[l]^2 for l = 1:n))

	
	@NLobjective(system_id, Min, T1 + 2*T2 + .5*g2 - 2/sqrt(N)*sum(gamma[l]*sqrt(T3[l] + 2*T4[l] + glk[l]) for l = 1:n) + b*sum(abs(A1[i,j]) for i = 1:n-1 for j = i+1:n))

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
	
	writedlm("data/"*id*"_l2_$(l).$(k)_obj.csv",objective_value(system_id),',')
#	writedlm("data/"*id*"_l2_$(l).$(k)_A1.csv",mL,',')
#	writedlm("data/"*id*"_l2_$(l).$(k)_a2.csv",value.(a2),',')
#	writedlm("data/"*id*"_l2_$(l).$(k)_gamma.csv",value.(gamma),',')
end

#=
"""
	Lmin_l2_lap_par(id::String, x::Array{Float64,2}, Dx::Array{Float64,2}, xt::Array{Complex{Float64},2}, Dxt::Array{Complex{Float64,2}, k::Int64, A1h::Array{Float64,2}, a1h::Array{Float64,1}, b::Float64=0., mu::Float64=1e-1, bp::Float64=1e-1)

Same as Lmin_l2_par, but assumes that the dynamics matrix (A1) is a Laplacian, with nonpositive off-diagonal terms.
"""
=#
@everywhere function Lmin_l2_lap_par(tups::Tuple{String,Array{Float64,2},Array{Float64,2},Array{Complex{Float64},2},Array{Complex{Float64},2},Int64,Array{Float64,2},Array{Float64,1},Float64,Float64,Float64})
	id,x,Dx,xt,Dxt,k,A1h,a2h,b,mu,bp = tups
	
	t0 = time()
	
	nn,N = size(x)
	n = Int(nn/2)

# Definition of the needed parameters.
	S0 = (x*x')/N
	S1 = (x*Dx')/N

	xtk = xt[:,k+1]
	
	Fk = real.(xtk*xtk')
	flk = [real.(Dxt[l,k]*xtk') for l = 1:n]
	glk = [norm(Dxt[l,k])^2 for l = 1:n]

# Definition of the optimization problem.
	system_id = Model(optimizer_with_attributes(Ipopt.Optimizer, "mu_init" => mu, "bound_push" => bp))
	@variable(system_id, A1[i = 1:n-1, j = i+1:n])
	for i = 1:n-1
		for j = i+1:n
			@constraint(system_id, A1[i,j] <= 0.)
			set_start_value(A1[i,j],A1h[i,j])
		end
	end
	@variable(system_id, a2[i = 1:n])
	set_start_value.(a2,a2h)
	@variable(system_id, gamma[l = 1:n])
	for l in 1:n
		@constraint(system_id, gamma[l] >= 0.)
	end
	set_start_value.(gamma,ones(n))

# tr(L^2*\th*\th^T)
	@NLexpression(system_id, AtAS01, 
		      sum(A1[j,i]*A1[m,j]*(S0[m,i]-S0[j,i]) for i = 1:n for j = 1:i-1 for m = 1:j-1) + 
		      sum(A1[j,i]*A1[j,m]*(S0[m,i]-S0[j,i]) for i = 1:n for j = 1:i-1 for m = j+1:n) -
		      sum(A1[j,i]*A1[m,i]*(S0[m,i]-S0[i,i]) for i = 1:n for j = 1:i-1 for m = 1:i-1) - 
		      sum(A1[j,i]*A1[i,m]*(S0[m,i]-S0[i,i]) for i = 1:n for j = 1:i-1 for m = i+1:n) +
		      sum(A1[i,j]*A1[m,j]*(S0[m,i]-S0[j,i]) for i = 1:n for j = i+1:n for m = 1:j-1) + 
		      sum(A1[i,j]*A1[j,m]*(S0[m,i]-S0[j,i]) for i = 1:n for j = i+1:n for m = j+1:n) -
		      sum(A1[i,j]*A1[m,i]*(S0[m,i]-S0[i,i]) for i = 1:n for j = i+1:n for m = 1:i-1) -
		      sum(A1[i,j]*A1[i,m]*(S0[m,i]-S0[i,i]) for i = 1:n for j = i+1:n for m = i+1:n))

# tr(-L*D*\om*\th^T) = tr(-D*L*\th*\om^T)
	@NLexpression(system_id, AtAS02, 
		      sum(A1[j,i]*(a2[j]*S0[n+j,i]-a2[i]*S0[n+i,i]) for i = 1:n for j = 1:i-1) + 
		      sum(A1[i,j]*(a2[j]*S0[n+j,i]-a2[i]*S0[n+i,i]) for i = 1:n for j = i+1:n))

# tr(D^2*\om*\om^T)
	@NLexpression(system_id, AtAS03, sum(a2[i]^2*S0[n+i,n+i] for i = 1:n))
	
# tr(A^T*A*S0)
	@NLexpression(system_id, T1, AtAS01 + 2*AtAS02 + AtAS03)


# tr(-L*\th*(D\om)^T)
	@NLexpression(system_id, AS11, 
		      sum(A1[i,j]*(S1[j,i]-S1[i,i]) for i = 1:n for j = i+1:n) + 
		      sum(A1[j,i]*(S1[j,i]-S1[i,i]) for i = 1:n for j = 1:i-1))

# tr(-D*\om*(D\om)^T)
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


	@NLexpression(system_id, g2, sum(gamma[l]^2 for l = 1:n))

	
	@NLobjective(system_id, Min, T1 + 2*T2 + .5*g2 - 2/sqrt(N)*sum(gamma[l]*sqrt(T3[l] + 2*T4[l] + glk[l]) for l = 1:n) + b*sum(abs(A1[i,j]) for i = 1:n-1 for j = i+1:n))

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

	@info "Full optimization took $(time() - t0)''."

	writedlm("data/"*id*"_lap2_$(l).$(k)_obj.csv",objective_value(system_id),',')
#	writedlm("data/"*id*"_lap2_$(l).$(k)_A1.csv",mL,',')
#	writedlm("data/"*id*"_lap2_$(l).$(k)_a2.csv",value.(a2),',')
#	writedlm("data/"*id*"_lap2_$(l).$(k)_gamma.csv",value.(gamma),',')
end

#=
"""
    get_Ah_correl(Xs::Array{Float64,2}, dt::Float64)

Estimates the dynamics matrix bases on Lokhov18.

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

	S0 = Xs*Xs' ./ T
	S1 = Xs[:,2:T]*Xs[:,1:T-1]' ./ (T-1)

	Ah = S1*inv(S0)
	Lh = -Ah[n+1:2*n,1:n]/dt
	dh = (ones(n) - diag(Ah[n+1:2*n,n+1:2*n]))./dt

	return Ah,Lh,dh
end

