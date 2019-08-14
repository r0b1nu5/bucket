

function find_cos(X::Array{Float64,1}, dt::Float64)
	T = length(X)

	find_afp = Model(with_optimizer(Ipopt.Optimizer))

	@variable(find_afp, a >= 0)
	@variable(find_afp, f >= 0)
	@variable(find_afp, 0 <= p <= 2pi)

	@NLexpression(find_afp, err[t = 1:T], X[t] - a*cos(2*pi*f*dt*t + p))

	@NLobjective(find_afp, Min, sum(err[t]^2 for t = 1:T))
	
	optimize!(find_afp)

	return (value(a), value(f), value(p))
end









