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

