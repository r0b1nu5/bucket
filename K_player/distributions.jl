using SpecialFunctions

# Normalization for power law
function C_pl(s::Float64, mi::Int64, ma::Int64)
	return 1/zeta(s,mi)
end

# Value of the power law distribution
function z_pl(s::Float64, mi::Int64, ma::Int64)
	C = C_pl(s,mi,ma)

	return C * ((mi:ma).^(-s))
end

# Normalization for power law with cutoff
function C_plc(a::Float64, l::Float64, mi::Int64, ma::Int64)
	return 1/(real(polylog(a,Complex(exp(-l)))) - sum((1:mi-1).^(-a).*exp.(-l*(1:mi-1))))
end

# Value of the power law with cutoff distribution
function z_plc(a::Float64, l::Float64, mi::Int64, ma::Int64)
	C = C_plc(a,l,mi,ma)

	return C * (mi:ma).^(-a) .* exp.(-l.*(mi:ma))
end

# Normalization for Yule-Simon
function C_ys(a::Float64, mi::Int64, ma::Int64)
	return 1/(1-(a-1)*sum(beta.(1:(mi-1),a)))
end

# Value of the Yule-Simon distribution
function z_ys(a::Float64, mi::Int64, ma::Int64)
	C = C_ys(a,mi,ma)

	return C * (a-1)*beta.(mi:ma,a)
end

# Normalization for exponential
function C_exp(b::Float64, mi::Int64, ma::Int64)
	return (1 - exp(-b))/exp(-b*mi)
end

# Value of the exponential distribution
function z_exp(s::Float64, mi::Int64, ma::Int64)
	C = C_exp(b,mi,ma)

	return C * exp.(-b*(mi:ma))
end




