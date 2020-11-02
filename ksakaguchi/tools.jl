function winding(th::Array{Float64,1}, C::Array{Int64,1})
	th1 = th[C]
	th2 = th[[C[2:end];C[1]]]

	dth = th2 - th1

	q = sum(mod.(dth .+ pi,2pi) .- pi)/(2pi)

	return q
end

function jacobian(L::Array{Float64,2}, th::Array{Float64,1}, a::Float64)
	n = length(th)

	A = -L.*(1 .- diagm(0 => ones(n)))

	dth = th*ones(1,n) - ones(n)*th'

	J = A.*cos.(dth .- a)
end






