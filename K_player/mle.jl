# We fit the data with a Maximum-Likelihood Estimator (MLE)


function new_mle_pl(x::Array{Int64,2})
	mi = minimum(x[1,:])
	ma = maximum(x[1,:])
	n = sum(x[2,:])
	
	l = 1.01
	u = 10.
	s = 10.
	
	sum_data = sum(x[2,:].*log.(x[1,:]))
	
	error = 1000
	
	while error > 1e-5
		ss = LinRange(l,u,100)
		L = -n*log.([zeta(si,mi) for si in ss]) - ss.*sum_data
		id = findmax(L)[2]
		l = max(l,ss[max(id-1,1)])
		u = min(u,ss[min(id+1,100)])
		error = u - l
	end
	
	s = (l+u)/2
	
	return s
end


function new_mle_plc(x::Array{Int64,2})
	side = 2.99
	atemp = 3.
	ltemp = 3.
	
	L = zeros(100,100)
	
	while side > 1e-4
		as = LinRange(max(atemp-side,side/100),atemp+side,100)
		ls = LinRange(max(ltemp-side,side/100),ltemp+side,100)
		
		for i in 1:100
			for j in 1:100
				L[i,j] = -sum(x[2,:])*log(real(polylog(as[i],Complex(exp(-ls[j]))))) - sum(x[2,:].*(as[i].*log.(x[1,:]) + ls[j].*x[1,:]))
			end
		end
		
		atemp = as[findmax(L)[2][1]]
		ltemp = ls[findmax(L)[2][2]]
		side /= 10
	end
	
	return atemp,ltemp
end

function new_mle_yule(x::Array{Int64,2}, xmin::Int64)
	side = 4.99
	atemp = 6.
	n = sum(x[2,:])
	
	L = zeros(100)
	
	while side > 1e-4
		as = LinRange(max(atemp-side,1+side/100),atemp+side,100)
		
		for i in 1:100
			L[i] = -n*log(1-(as[i]-1)*sum(beta.(1:(xmin-1),as[i]))) + n*log(as[i]-1) + sum(x[2,:].*log.(beta.(x[1,:],as[i])))
		end
		
		atemp = as[findmax(L)[2]]
		side /= 10
	end
	
	return atemp
end


function new_mle_exp(x::Array{Int64,2}, xmin::Int64)
	side = 4.99
	btemp = 5.
	n = sum(x[2,:])

	L = zeros(100)

	while side > 1e-4
		bs = LinRange(max(btemp-side,side/100),btemp+side,100)
		
		for i in 1:100
			L[i] = n*(log(1-exp(-bs[i])) + bs[i]*xmin) - bs[i]*sum(x[2,:].*x[1,:])
		end

		btemp = bs[findmax(L)[2]]
		side /= 10
	end

	return btemp
end


function new_mle_poisson(x::Array{Int64,2}, xmin::Int64)
	side = 4.99
	mtemp = 5.0
	n = sum(x[2,:])

	L = zeros(100)

	while side > 1e-4
		ms = LinRange(max(mtemp-side,0.),mtemp+side,100)

		for i in 1:100
			L[i] = -n*log(exp(ms[i]) - sum((ms[i].^(0:xmin-1))./(factorial.(0:xmin-1)))) + ms[i]*sum(x[2,:].*x[1,:]) - sum(x[2,:].*log_factorial(Int.(x[1,:])))
		end

		mtemp = ms[findmax(L)[2]]
		side /= 10
	end

	return mtemp
end


###################### SIDE FUNCTIONS #################################
	

function polylog(s::Float64,z::Complex{Float64})
	Li = gamma(1-s)/(2pi)^(1-s)*((im)^(1-s)*zeta(1-s,.5+log(-z)/(2pi*im)) + (im)^(s-1)*zeta(1-s,.5-log(-z)/(2pi*im)))
	return Li
end

# Computes log(x!)
function log_factorial(x::Int64)
	return sum(log.(1:x))
end

function log_factorial(x::Array{Int64,1})
	lf = Array{Float64,1}()
	for xx in x
		push!(lf,log_factorial(xx))
	end
	return lf
end

##################### OLD MLE #####################################
#
function mle_pl(x::Array{Float64,1})
	mi = minimum(x)
	ma = maximum(x)
	n = length(x)
	
	l = 1.01
	u = 10.
	s = 10.
	
	sum_data = sum(log.(x))
	
	error = 1000
	
	while error > 1e-5
		ss = LinRange(l,u,100)
		L = -n*log.([zeta(si,mi) for si in ss]) - ss.*sum_data
		id = findmax(L)[2]
		l = max(l,ss[max(id-1,1)])
		u = min(u,ss[min(id+1,100)])
		error = u - l
	end
	
	s = (l+u)/2
	
	return s
end


function mle_plc(x::Array{Float64,1})
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


function mle_yule(x::Array{Float64,1},xmin::Float64)
	side = 4.99
	atemp = 6.0
	n = length(x)

	L = zeros(100)

	while side > 1e-4
		as = LinRange(max(atemp-side,1+side/100),atemp+side,100)

		for i in 1:100
#			L[i] = n*log(as[i]-1) + n*log(gamma(xmin+as[i]-1)) + sum(log.(gamma.(x))) - sum(log.(gamma.(x.+as[i])))
			L[i] = -n*log(1-(as[i]-1)*sum(beta.(1:(xmin-1),as[i]))) + n*log(as[i]-1) + sum(log.(beta.(x,as[i])))
		end

		atemp = as[findmax(L)[2]]
		side /= 10
	end

	return atemp
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

