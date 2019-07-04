using JuMP, Ipopt

function locate_forced(X::Array{Float64,2})
	n,T = size(X)
	
	l = ones(n)
	
	system_id = Model(with_optimizer(Ipopt.Optimizer, print_level=0))
	
	@variables(system_id, begin
		A[i = 1:n, j = 1:n],
		c[i = 1:n], 
		f[i = 1:n], 
		phi[i = 1:n]
	end)
	
	@NLobjective(system_id, Min, sum(norm(X[t+1] - A*X[t] - c .* cos.(t * f + phi))^2 for t in 1:T-1) + sum(l .* c))
	
	for i in 1:n-1
		for j in i+1:n
			@constraint(system_id, A[i,j] == A[j,i])
		end
	end
	
	for i in 1:n
		@constraint(system_id, sum(A[i,:]) == 0)
	end
	
	optimize!(system_id)
	
	return (
		A = value.(A),
		c = value.(c),
		f = value.(f),
		phi = value.(phi)
	)
end

		
	
	
	




