# Computes the s coefficient of the Zipf's law with max likelihood for the data set x.

function ml_zipf(x::Array{Float64,1})
	mi = minimum(x)
	ma = maximum(x)
	
	l = .01
	u = 10.
	
	s = 10.
	
	f = Array{Float64,1}()
	for i in mi:ma
		push!(f,sum(x.==i)/length(x))
	end
	suf = sum(f)
	su = sum(f.*log(mi:ma))
	
	error = 1000
	
	while error > 1e-5
		ss = linspace(l,u,100)
		Hs = Array{Float64,1}()
		for s in ss
			push!(Hs,sum(1 ./((mi:ma).^s)))
		end
		La = -ss*su - suf*log(Hs)
		
		id = indmax(La)
		l = max(l,ss[max(id-1,1)])
		u = min(u,ss[min(id+1,100)])
		error = u - l
	end
	
	s = (l+u)/2

	return s
end

function ml_clauset(x::Array{Float64,1})
	mi = minimum(x)
	ma = maximum(x)
	n = length(x)
	
	l = .01
	u = 10.
	s = 10.
	
	sum_data = sum(log(x))
	
	error = 1000
	
	while error > 1e-5
		ss = linspace(l,u,100)
		zeta = Array{Float64,1}()
		for si in ss
			push!(zeta,sum((mi:ma).^(-si)))
		end
		
		L = -n*log(zeta) - ss.*sum_data
		
		id = indmax(L)
		l = max(l,ss[max(id-1,1)])
		u = min(u,ss[min(id+1,100)])
		error = u - l
	end
	
	s = (l+u)/2
	
	return s
end


