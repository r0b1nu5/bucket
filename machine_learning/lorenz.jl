function lorenz(x0::Array{Float64,1},n_iter::Int64=1000,h::Float64=.01,sig::Float64=10.,beta::Float64=8/3,rho::Float64=28.)
	xs = Array{Float64,2}(undef,3,0)
	xs = [xs x0]
	
	i = 0
	
	while i < n_iter
		i += 1
		if i%1000 == 0
			@info "$i/$n_iter"
		end
		
		xyz = xs[:,end]
		
		k1 = f_lorenz(xyz,sig,beta,rho)
		k2 = f_lorenz(xyz + h/2*k1,sig,beta,rho)
		k3 = f_lorenz(xyz + h/2*k2,sig,beta,rho)
		k4 = f_lorenz(xyz + h*k3,sig,beta,rho)
	
		dx = (k1+2*k2+2*k3+k4)/6
		
		xs = [xs xs[:,end]+h*dx]
	end
	
	return xs
end

function f_lorenz(xyz::Array{Float64,1},s::Float64=10.,b::Float64=8/3,r::Float64=28.)
	x,y,z = xyz

	return [s*(y - x), x*(r - z) - y, x*y - b*z]
end



