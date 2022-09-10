using SpecialFunctions

function xi_th(x::Float64, eps::Float64, mu::Float64, sig::Float64)
	xi = ((-sig^2)*(exp(-(x+eps-mu)^2/(2*sig^2)) - 
			exp(-(x-eps-mu)^2/(2*sig^2)) + 
			exp(-(x+eps+mu)^2/(2*sig^2)) - 
			exp(-(x-eps+mu)^2/(2*sig^2))) - 
	     sqrt(2*pi)*sig*((x-mu)*(erf(x+eps-mu) - 
				     erf(x-eps-mu)) + 
			     (x+mu)*(erf(x+eps+mu) - erf(x-eps+mu)))/2)/(2*sqrt(2*pi)*sig)

	return xi
end

function nm(x::Float64, eps::Float64, mu::Float64, sig::Float64)
	e1 = erf((x-mu)/(sqrt(2)*sig))
	e2 = erf((x-eps-mu)/(sqrt(2)*sig))
	e3 = erf((x+mu)/(sqrt(2)*sig))
	e4 = erf((x-eps+mu)/(sqrt(2)*sig))

	return .5*(e1 - e2 + e3 - e4)
end

function np(x::Float64, eps::Float64, mu::Float64, sig::Float64)
	e1 = erf((x+eps-mu)/(sqrt(2)*sig))
	e2 = erf((x-mu)/(sqrt(2)*sig))
	e3 = erf((x+eps+mu)/(sqrt(2)*sig))
	e4 = erf((x+mu)/(sqrt(2)*sig))

	return .5*(e1 - e2 + e3 - e4)
end

function int_non(D::Float64, eps::Float64, mu::Float64, sig::Float64, dx::Float64=.01)
	nx = round(Int,1/dx)
	
	xs = LinRange(0,D,nx)
	nps = [np(x,eps,mu,sig) for x in xs]
	nms = [nm(x,eps,mu,sig) for x in xs]

	X = nps./nms.*(exp.(-(xs .- mu).^2/(2*sig^2)) + exp.(-(xs .+ mu).^2/(2*sig^2)))/(2*sqrt(2*pi)*sig)
	
	return sum(X)*dx
end

function int_nmn(D::Float64, eps::Float64, mu::Float64, sig::Float64, dx::Float64=.01)
	nx = round(Int,1/dx)
	
	xs = LinRange(-D,0,nx)
	nps = [np(x,eps,mu,sig) for x in xs]
	nms = [nm(x,eps,mu,sig) for x in xs]

	X = (nps - nms).*(exp.(-(xs .- mu).^2/(2*sig^2)) + exp.(-(xs .+ mu).^2/(2*sig^2)))/(2*sqrt(2*pi)*sig)
	
	return sum(X)*dx
end

function int_npn(D::Float64, eps::Float64, mu::Float64, sig::Float64, dx::Float64=.01)
	nx = round(Int,1/dx)
	
	xs = LinRange(-D,0,nx)
	nps = [np(x,eps,mu,sig) for x in xs]
	nms = [nm(x,eps,mu,sig) for x in xs]

	X = (nps + nms).*(exp.(-(xs .- mu).^2/(2*sig^2)) + exp.(-(xs .+ mu).^2/(2*sig^2)))/(2*sqrt(2*pi)*sig)
	
	return sum(X)*dx
end

