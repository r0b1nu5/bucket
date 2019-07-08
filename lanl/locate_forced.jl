using JuMP, Ipopt, LinearAlgebra, DelimitedFiles, Dates

## Identification of the network Laplacian, using the covariance matrix method, proposed in Lokhov18

function system_identification_ipopt(X::Array{Float64,2}, dt::Float64)
	@info "$(now()) -- Start..."
	
	nn,T = size(X)
	n = Int(nn/2)
	
	Mi = diagm(0 => 1 ./ m)
	D = diagm(0 => d)
	l = .01

	S0 = zeros(2*n,2*n)
	S1 = zeros(2*n,2*n)
	
	for i in 1:T-1
		S0 += Xs[:,i]*Xs[:,i]' ./ (T-1)
		S1 += Xs[:,i+1]*Xs[:,i]' ./ (T-1)
	end
	
	Ah = S1*inv(S0)
	Adh = (Ah - diagm(0 => ones(2*n)))./dt
	
	@info "$(now()) -- Stop."
	
	return Adh
end


function locate_forcing(X::Array{Float64,2}, A::Array{Float64,2}, dt::Float64, nf::Int64 = 1)
	@info "$(now())"," -- Start..."
	
	nn,T = size(X)
	n = Int(nn/2)
	
	AX = [zeros(n,T); A[n+1:2*n,:]*X]
	
	locate_f = Model(with_optimizer(Ipopt.Optimizer, print_level=2, mumps_mem_percent=5))
	
## Variables
	@variable(locate_f, 0 <= c[i = 1:n] <= 5)
	@variable(locate_f, f[i = 1:n])
	@variable(locate_f, 0 <= phi[i = 1:n] <= 2*pi)
	
## Constraints
	@NLconstraint(locate_f, sum((c[i] > 0) for i = 1:n) == nf)
	
## Objective
	@NLobjective(system_id, Min, sum((X[n+i,t+1] - AX[n+i,t] - c[i]*cos(f[i]*t + phi[i]))^2 for i = 1:n for t = 1:T-1))
	
	optimize!(locate_f)
	
	@info "$(now())"," -- Stop."
	
	return (
		c = value.(c),
		f = value.(f),
		phi = value.(phi)
	)
end


function find_inertia_damping(





## Identification of the network Laplacian by an optimization formulation (Ipopt).
## !!!!!!!!!!!!!!!!!!! FOR TOO LARGE NETWORK/DATASET, IT FILLS THE RAM AND CRASHES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!¨

function system_identification_ipopt(X::Array{Float64,2}, m::Array{Float64,1}, d::Array{Float64,1}, dt::Float64)
	@info "$(now()) -- Start..."
	
	nn,T = size(X)
	n = Int(nn/2)
	
	Mi = diagm(0 => 1 ./ m)
	D = diagm(0 => d)
	l = .01

	system_id = Model(with_optimizer(Ipopt.Optimizer, print_level=2, mumps_mem_percent=5))
		
## Variables
	@variables(system_id, begin
		L[i = 1:n, j = 1:n]
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
#		c = value.(c),
#		f = value.(f),
#		phi = value.(phi)
	)
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
	




