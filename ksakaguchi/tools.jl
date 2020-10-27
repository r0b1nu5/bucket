function winding(th::Array{Float64,1}, C::Array{Int64,1})
	th1 = th[C]
	th2 = th[[C[2:end];C[1]]]

	dth = th2 - th1

	q = sum(mod.(dth .+ pi,2pi) .- pi)/(2pi)

	return q
end





