using JuMP, Ipopt

function find_cos(X::Array{Float64,1}, dt::Float64)
	T = length(X)

	find_afp = Model(with_optimizer(Ipopt.Optimizer, mu_init=1e-5))

	@variable(find_afp, a >= 0, start = (maximum(X) - minimum(X))/2)
	@variable(find_afp, f >= 0)
	@variable(find_afp, 0 <= p <= 2pi)
	@variable(find_afp, sh)

	@NLexpression(find_afp, err[t = 1:T], X[t] - a*cos(2*pi*f*dt*t + p) - sh)

	@NLobjective(find_afp, Min, sum(err[t]^2 for t = 1:T))
	
	optimize!(find_afp)

	return (value(a), value(f), value(p), value(sh))
end


function find_sin(X::Array{Float64,1}, dt::Float64)
	T = length(X)

	find_afp = Model(with_optimizer(Ipopt.Optimizer, mu_init=1e-5))

	@variable(find_afp, a >= 0, start = (maximum(X) - minimum(X))/2)
	@variable(find_afp, f >= 0)
	@variable(find_afp, 0 <= p <= 2pi)
	@variable(find_afp, sh)

	@NLexpression(find_afp, err[t = 1:T], X[t] - a*sin(2*pi*f*dt*t + p) - sh)

	@NLobjective(find_afp, Min, sum(err[t]^2 for t = 1:T))
	
	optimize!(find_afp)

	return (value(a), value(f), value(p), value(sh))
end







