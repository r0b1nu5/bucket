using LinearAlgebra, JuMP, Ipopt, Dates

## INPUT 
# X: time series of the dynamic variables. In our case, rows 1 to n are the angles and rows n+1 to 2*n are the velocities. Each column is a time step.
# m: Vector of inertias
# d: Vector of dampings
# dt: time step

## OUTPUT
# L: estimated Laplacian of the network.


function system_identification_ipopt(X::Array{Float64,2}, m::Array{Float64,1}, d::Array{Float64,1}, dt::Float64)
	
	T = size(X)[2]
	n = length(m)
	
	Mi = diagm(0 => 1 ./ m)
	D = diagm(0 => d)

	system_id = Model(with_optimizer(Ipopt.Optimizer))
@info "$(now()) -- Model is defined"

## Variables
	@variable(system_id, L[i = 1:n, j = 1:n])
@info "$(now()) -- Variables are defined"

# #=	
## Constraints
	for i in 1:n-1
		for j in i+1:n
			@constraint(system_id, L[i,j] <= 0)
			@constraint(system_id, L[i,j] == L[j,i])
		end
	end
@info "$(now()) -- Symmetry constraints defined"
	
	for i in 1:n
		@constraint(system_id, L[i,i] == -sum(L[i,j] for j = 1:i-1) - sum(L[i,j] for j = i+1:n))
	end
@info "$(now()) -- Zero row-sum constraints defined"
# =#
	
## Objective
# Defining the error at each node and each time step
#err = [X[n+i,t+1] - X[n+i,t] - dt * (sum(-L[i,k]*X[k,t] for k = 1:n)/m[i] - d[i]*X[n+i,t]/m[i]) for i in 1:n, for t in 1:T-1]
@info "$(now()) -- Intermediate expressions defined"

# Sum of squared errors
#	@objective(system_id, Min, sum(err[i,t]^2 for i in 1:n for t in 1:T-1)) 
	@objective(system_id, Min, sum((X[n+i,t+1] - X[n+i,t] - dt * (sum(-L[i,k]*X[k,t] for k = 1:n)/m[i] - d[i]*X[n+i,t]/m[i]))*(X[n+i,t+1] - X[n+i,t] - dt * (sum(-L[i,k]*X[k,t] for k = 1:n)/m[i] - d[i]*X[n+i,t]/m[i])) for i in 1:n for t in 1:T-1))

@info "$(now()) -- Objective function defined"

	optimize!(system_id)
@info "$(now()) -- Optimization completed"

	return value.(L)
end


function system_identification_ipopt_no_constraints(X::Array{Float64,2}, m::Array{Float64,1}, d::Array{Float64,1}, dt::Float64)
	
	T = size(X)[2]
	n = length(m)
	
	Mi = diagm(0 => 1 ./ m)
	D = diagm(0 => d)

	system_id = Model(with_optimizer(Ipopt.Optimizer))
@info "$(now()) -- Model is defined"

## Variables
	@variable(system_id, L[i = 1:n, j = 1:n])
@info "$(now()) -- Variables are defined"

	
## Objective
# Defining the error at each node and each time step
	@NLexpression(system_id, err[i = 1:n, t = 1:T-1], X[n+i,t+1] - X[n+i,t] - dt * (sum(-L[i,k]*X[k,t] for k = 1:n)/m[i] - d[i]*X[n+i,t]/m[i]))
@info "$(now()) -- Intermediate expressions defined"

# Sum of squared errors
	@NLobjective(system_id, Min, sum(err[i,t]^2 for i in 1:n for t in 1:T-1)) 
@info "$(now()) -- Objective function defined"

	optimize!(system_id)
@info "$(now()) -- Optimization completed"

	return value.(L)
end
