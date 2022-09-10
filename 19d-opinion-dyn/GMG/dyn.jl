using LinearAlgebra

function dyn(L::Array{Float64,2}, x0::Array{Float64,2}, w0::Array{Float64,2}, xi::Array{Float64,2}, max_iter::Int64=10000, eps::Float64=1e-8, h::Float64=.1)
	n,p = size(x0)
	
	D = diagm(0 => diag(L))
	Ld = diagm(0 => 1 ./ diag(L))*L

	M = Ld + diagm(0 => ones(n))

	err = 1000.
	iter = 0

	x1 = copy(xi)
	x2 = copy(xi)

	xs = Array{Float64,3}(undef,n,p,0)
	xs = cat(xs,xi,dims=3)

	while err > eps && iter < max_iter
		iter += 1
		if iter%1000 == 0
			@info "$iter"
		end

		x1 = copy(x2)

		k1 = -M*x1 + (x0+w0)
		k2 = -M*(x1 + h/2*k1) + (x0+w0)
		k3 = -M*(x1 + h/2*k2) + (x0+w0)
		k4 = -M*(x1 + h*k3) + (x0+w0)

		dx = (k1+2*k2+2*k3+k4)/6

		x2 = x1 + h*dx

		xs = cat(xs,x2,dims=3)
		
		err = maximum(abs.(dx))
	end

	return xs
end







