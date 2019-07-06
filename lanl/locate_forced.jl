using JuMP, Ipopt, LinearAlgebra, DelimitedFiles, Dates


function locate_forced(X::Array{Float64,2}, m::Array{Float64,1}, d::Array{Float64,1}, dt::Float64)
	@info "$(now()) -- Start..."
	
	nn,T = size(X)
	n = Int(nn/2)
	
	M = diagm(0 => m)
	D = diagm(0 => d)
	l = .01
		
	system_id = Model(with_optimizer(Ipopt.Optimizer, print_level=2))
	
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

		
	
	
	




