function lorentz(x0::Array{Float64,1},n_iter::Int64=1000,h::Float64=.01,sig::Float64=10.,beta::Float64=8/3,rho::Float64=28.)
	xs = Array{Float64,2}(undef,3,0)
	xs = [xs x0]
	
	i = 0
	
	while i < n_iter
		i += 1
		if i%1000 == 0
			@info "$i/$n_iter"
		end
		
		x,y,z = xs[:,end]
		
		k1 = [sig*(y-x), 
		      x*(rho-z)-y, 
		      x*y-beta*z]
		k2 = [sig*((y+h/2*k1[2])-(x+h/2*k1[1])), 
		      (x+h/2*k1[1])*(rho-(z+h/2*k1[3]))-(y+h/2*k1[2]), 
		      (x+h/2*k1[1])*(y+h/2*k1[2])-beta*(z+h/2*k1[3])]
		k3 = [sig*((y+h/2*k2[2])-(x+h/2*k2[1])), 
		      (x+h/2*k2[1])*(rho-(z+h/2*k2[3]))-(y+h/2*k2[2]), 
		      (x+h/2*k2[1])*(y+h/2*k2[2])-beta*(z+h/2*k2[3])]
		k4 = [sig*((y+h*k3[2])-(x+h*k3[1])), 
		      (x+h*k3[1])*(rho-(z+h*k3[3]))-(y+h*k3[2]), 
		      (x+h*k3[1])*(y+h*k3[2])-beta*(z+h*k3[3])]
		
		dx = (k1+2*k2+2*k3+k4)/6
		
		xs = [xs xs[:,end]+h*dx]
	end
	
	return xs
end





