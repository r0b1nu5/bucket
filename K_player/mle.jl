

function mle_pl_cutoff(x::Array{Float64,1})
	side = 2.99
	atemp = 3.
	ltemp = 3.
	
	L = zeros(100,100)
	
	while side > 1e-4
		as = LinRange(max(atemp-side,side/100),atemp+side,100)
		ls = LinRange(max(ltemp-side,side/100),ltemp+side,100)
		
		for i in 1:100
			for j in 1:100
				L[i,j] = -length(x)*log(real(polylog(as[i],Complex(exp(-ls[j]))))) - sum(as[i].*log.(x) + ls[j].*x)
			end
		end
		
		atemp = as[findmax(L)[2][1]]
		ltemp = ls[findmax(L)[2][2]]
		side /= 10
	end
	
	return atemp,ltemp
end

function mle_exp(x::Array{Float64,1},xmin::Float64)
	side = 4.99
	ltemp = 5.
	
	L = zeros(100)

	while side > 1e-4
		ls = LinRange(max(ltemp-side,side/100),ltemp+side,100)

		for i in 1:100
			L[i] = length(x)*(1-exp(-ls[i])) + ls[i]*xmin - ls[i]*sum(x)
		end

		ltemp = ls[findmax(L)[2]]
		side /= 10
	end

	return ltemp
end

function mle_yule(x::Array{Float64,1},xmin::Float64)
	side = 4.99
	atemp = 6.0
	n = length(x)

	L = zeros(100)

	while side > 1e-4
		as = LinRange(max(atemp-side,1+side/100),atemp+side,100)

		for i in 1:100
#			L[i] = n*log(as[i]-1) + n*log(gamma(xmin+as[i]-1)) + sum(log.(gamma.(x))) - sum(log.(gamma.(x.+as[i])))
			L[i] = n*log(as[i]-1) + sum(log.(beta.(x,as[i])))
		end

		atemp = as[findmax(L)[2]]
		side /= 10
	end

	return atemp
end

function mle_poisson(x::Array{Float64,1},xmin::Float64)
	side = 4.99
	mtemp = 5.0
	n = length(x)

	L = zeros(100)

	while side > 1e-4
		ms = LinRange(max(mtemp-side,0.0),mtemp+side,100)

		for i in 1:100
			L[i] = -n*log(exp(ms[i])-sum(ms[i].^(0:xmin-1)./(factorial.(0:xmin-1)))) .+ sum(log(ms[i]).*x)
		end

		mtemp = ms[findmax(L)[2]]
		side /= 10
	end

	return mtemp
end


	
	

function polylog(s::Float64,z::Complex{Float64})
	Li = gamma(1-s)/(2pi)^(1-s)*((im)^(1-s)*zeta(1-s,.5+log(-z)/(2pi*im)) + (im)^(s-1)*zeta(1-s,.5-log(-z)/(2pi*im)))
	return Li
end

