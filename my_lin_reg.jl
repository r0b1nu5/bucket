function lin_reg(x::Array{Float64,1},y::Array{Float64,1})
	n = length(x)
	
	Sx = sum(x)
	Sy = sum(y)
	Sxx = sum(x.^2)
	Sxy = sum(x.*y)

	a = (Sxy - Sx*Sy/n)/(Sxx - Sx^2/n)
	b = Sy/n - a*Sx/n

	return a,b
end



